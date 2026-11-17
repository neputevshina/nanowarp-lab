package main

import (
	"flag"
	"math"
	"math/cmplx"
	"os"

	"github.com/neputevshina/nanowarp-lab/common"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/dspio/wavio"
	"github.com/neputevshina/nanowarp/oscope"
	"github.com/neputevshina/nanowarp/pffft"
	"golang.org/x/exp/constraints"
)

var finputa = flag.String("a", "", "source WAV")
var fhrad = flag.Int("h", 0, "horizontal radius, if less than or 0 equal to -l")
var fvrad = flag.Int("v", 13, "vertical radius")
var folap = flag.Int("l", 5, "overlap")
var fnfft = flag.Int("nfft", 3600, "fft grain size")

var fvib = flag.Int("vib", 7, "vertical vibrato radius, 0 to disable")
var ftonal = flag.Bool("t", false, "output tonals")
var fhperc = flag.Int("ph", 50, "horizontal percentile")
var fvperc = flag.Int("pv", 50, "vertical percentile")

var fcif = flag.Bool("cif", false, "backtrack hpss mask by cif displacement")

var fwin = flag.String("w", `kaiser`, "use a `window`: hann, welch, vorbis, tri, kaiser")

func main() {
	flag.Parse()
	oscope.Enable = true

	fa, err := os.Open(*finputa)
	if err != nil {
		panic(err)
	}
	name := "./perc.wav"
	if *ftonal {
		name = "./harm.wav"
	}
	fo, err := os.Create(name)
	if err != nil {
		panic(err)
	}

	nfft := pffft.NearestNfft(*fnfft)
	hop := nfft / *folap
	nbins := nfft/2 + 1

	wavr, err := wavio.NewDecoder(fa)
	if err != nil {
		panic(err)
	}
	p := wavr.Properties()
	wavw, err := wavio.NewEncoder(fo, p.Samplerate, p.Nch, wavio.FormatFloat, 32)
	if err != nil {
		panic(err)
	}
	defer wavw.Close()

	gr := dspio.NewOfflineGrainReader(nfft, hop, wavr)
	fft := pffft.New(nfft)

	// HPSS works best with windows with the smallest main lobe width,
	// but a window's derivative must have compact support (zeros only before and after the window)
	// to be perfectly invertible after modifications.
	// Meaning, a rectangular or a Tukey window will produce artifacts.
	win := make([]float64, nfft)
	switch *fwin {
	case `hann`:
		hann(win)
	case `welch`:
		for w := range nfft {
			x := unmix(0, float64(nfft), float64(w))
			win[w] = 1 - (2*x-1)*(2*x-1) // Welch, lags
		}
	case `vorbis`:
		for w := range nfft {
			x := unmix(0, float64(nfft), float64(w))
			c := math.Cos(math.Pi * (x + 0.5))
			win[w] = math.Sin(math.Pi / 2 * c * c) // Vorbis
			// win[w] = math.Sqrt(1 - x*x) // Circular, doesn't work
		}
	case `tri`:
		for w := range nfft / 2 {
			win[w] = unmix(0, float64(nfft/2), float64(w))
		}
		for w := range nfft / 2 {
			win[nfft/2:][w] = unmix(float64(nfft/2), 0, float64(w))
		}
	case `kaiser`:
		win = kaiser(nfft, 3)
	default:
		panic(`unknown window type`)
	}

	if *fhrad == 0 {
		*fhrad = *folap
	}

	bank := make([]*mediator[float64, struct{}], nbins)
	for i := range bank {
		bank[i] = mediatorNew[float64, struct{}](*fhrad, *fhrad, float64(*fhperc)/100)
	}
	vert := mediatorNew[float64, struct{}](*fvrad, *fvrad, float64(*fvperc)/100)
	smolv := mediatorNew[float64, struct{}](*fvib, *fvib, 1)
	smolbank := make([]*mediator[float64, struct{}], nbins)
	for i := range smolbank {
		smolbank[i] = mediatorNew[float64, struct{}](*folap, *folap, 1)
	}
	magt := make([]float64, nbins)
	magt1 := make([]float64, nbins)
	first := true
	now := *fhrad
	vib := func(out []float64, in [][]complex128, ok bool) {
		if *fvib != 0 && ok {
			smolv.Reset(*fvib)
			for w := *fvib - 1; w >= 0; w-- { // Mirror the beginning of vibrato removal
				mid := .0
				for ch := range in {
					mid += cmplx.Abs(in[ch][w])
				}
				smolv.Insert(mid, struct{}{})
			}
		}
		for w := range nbins {
			mid := .0
			for ch := range in {
				mid += cmplx.Abs(in[ch][w])
			}
			if *fvib != 0 && ok {
				out[w], _ = smolv.Filt(mid, struct{}{})
			} else {
				out[w] = mid
			}
		}
	}
	roto := make2[bool](now*2+1, nbins)
	newroto := make([]int, nbins)
	common.StftHandle(gr, wavw, fft, [][]float64{win, windowDx(win), windowT(win)}, now, now, func(out [][]complex128, windowed [][][][]complex128, preanalyze bool) {
		frames := windowed[0]
		// Accumulate lookahead.
		if first {
			for f := range frames[now:] {
				vib(magt, frames[f], false)
				for w := range nbins {
					bank[w].Insert(magt[w], struct{}{}) // Insert two times imitating the mirroring.
					bank[w].Insert(magt[w], struct{}{}) // Order probably won't matter.
				}
			}
			first = false
		}
		if !preanalyze {
			vib(magt1, frames[now], false)
			vib(magt, frames[now], true)
			vert.Reset(*fvrad)
			for w := *fvrad - 1; w >= 0; w-- { // Mirror the beginning
				vert.Insert(magt1[w], struct{}{})
			}
			for w := range nbins {
				v, _ := vert.Filt(magt1[w], struct{}{})
				h, _ := bank[w].Take()
				roto[now][w] = v > h
				newroto[w] = int(boolfloat(roto[now][w]))
			}
			if *fcif {
				// Backpropagate the mask vertically by vertical bin shift (CIF).
				clear(newroto)
				for w := range nbins {
					cif := imag(windowed[1][now][0][w]/frames[now][0][w]) / math.Pi / 2
					y := max(0, min(int(math.Round(float64(w)+cif)), nbins-1))
					if roto[now][y] {
						newroto[w]++
					} else {
						newroto[w]--
					}
					if roto[now][w] {
						newroto[y]++
					} else {
						newroto[y]--
					}
				}
			}
			for w := range nbins {
				m := boolfloat(newroto[w] > 0)
				if *ftonal {
					m = 1 - m
				}
				for ch := range out {
					out[ch][w] = frames[now][ch][w] * complex(m, 0)
				}
			}

			roto[len(roto)-1], roto[0] = roto[0], roto[len(roto)-1]
			copy(roto, roto[1:])
		}
		for w := range nbins { // Vibrato filtering
			bank[w].Insert(magt[w], struct{}{})
		}
	})

	oscope.Dump(nil, ".")
}

const (
	right = 1 << iota
	down
	up
	topleft
	topright
	ridgemask = right | down | up | topleft | topright
)

func fill[T any](s []T, e T) {
	for i := range s {
		s[i] = e
	}
}

func boolfloat(b bool) float64 {
	if b {
		return 1
	}
	return 0
}

func make2[T any](j, i int) (v [][]T) {
	v = make([][]T, j)
	for j := range j {
		v[j] = make([]T, i)
	}
	return
}

func make3[T any](k, j, i int) (v [][][]T) {
	v = make([][][]T, k)
	for k := range k {
		v[k] = make([][]T, j)
		for j := range j {
			v[k][j] = make([]T, i)
		}
	}
	return
}

func rotate[T any](frames []T) {
	last := len(frames) - 1
	t := frames[0]
	copy(frames, frames[1:])
	frames[last] = t
}

func unmix[F constraints.Float](a, b, x F) F {
	return (x - a) / (b - a)
}

func mix[F constraints.Float](a, b, x F) F {
	return a*(1-x) + b*x
}

// hann is a Hann window function.
//
// Dual to itself (achieves COLA/perfect reconstruction) at any integer overlap greater than 2.
func hann(out []float64) {
	for i := range out {
		x := float64(i) / float64(len(out))
		out[i] = 0.5 * (1 - math.Cos(2*math.Pi*x))
	}
}

// windowT calculates the time ramp multiplied window
// function for time-frequency reassignment.
func windowT(w []float64) (out []float64) {
	out = make([]float64, len(w))
	n := float64(len(w))
	for i := range w {
		out[i] = w[i] * mix(-n/2, n/2+1, float64(i)/n)
	}
	return
}

// windowDx calculates the derivative of w using DFT.
func windowDx(w []float64) (out []float64) {
	out = make([]float64, len(w))
	f := pffft.New(len(w))
	s := f.Coefficients(nil, w)
	for i := range s {
		s[i] *= complex(0, float64(i))
		s[i] /= complex(float64(len(w)), 0)
	}
	s[len(w)/2] = 0
	f.Sequence(out, s)
	for i := range out {
		out[i] = -out[i] * math.Pi * 2
	}
	return
}

func grainwrap(nfft int, x float64) float64 {
	nf := float64(nfft)
	r := math.Mod(x+nf/2, nf)
	if r < 0 {
		r += nf
	}
	r -= nf / 2
	if r != r || math.IsInf(r, 0) {
		return 0
	}
	return r
}
func clamp[T constraints.Ordered](a, b, x T) T {
	return max(a, min(b, x))
}
