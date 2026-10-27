package main

import (
	"flag"
	"math"
	"math/cmplx"
	"os"
	"slices"

	"github.com/neputevshina/nanowarp-lab/common"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/dspio/wavio"
	"github.com/neputevshina/nanowarp/oscope"
	"github.com/neputevshina/nanowarp/pffft"
	"golang.org/x/exp/constraints"
)

var finputa = flag.String("a", "", "source WAV")
var fhrad = flag.Int("h", 6, "horizontal radius")
var fvrad = flag.Int("v", 24, "vertical radius")

var fvib = flag.Int("vib", 0, "vertical vibrato radius")
var ftonal = flag.Bool("t", false, "output tonals")
var fhperc = flag.Int("ph", 50, "horizontal percentile")
var fvperc = flag.Int("pv", 50, "vertical percentile")

var fcif = flag.Bool("cif", false, "backtrack hpss mask by cif displacement")
var fseam = flag.Bool("seam", false, "output seam carved clips")

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

	const nfft, hop = 3840, 640
	const nbins = nfft/2 + 1

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

	gr := dspio.NewGrainReader(nfft, hop, wavr)
	fft := pffft.New(nfft)

	// HPSS works best with windows with the smallest main lobe width,
	// but a window's derivative must have compact support (zeros only before and after the window)
	// to be perfectly invertible after modifications.
	// Meaning, a rectangular or a Tukey window will produce artifacts.
	win := make([]float64, nfft)
	hann(win)

	zeroh, zerov := *fhrad**fhperc/100 == 0, *fvrad**fvperc/100 == 0
	bank := make([]*mediator[float64, struct{}], nbins)
	for i := range bank {
		if zeroh {
			bank[i] = mediatorNew[float64, struct{}](*fhrad, *fhrad, 1)
		} else {
			bank[i] = mediatorNew[float64, struct{}](*fhrad, *fhrad, float64(*fhperc)/100)
		}
	}
	vert := mediatorNew[float64, struct{}](*fvrad, *fvrad, float64(*fvperc)/100)
	if zerov {
		vert = mediatorNew[float64, struct{}](*fvrad, *fvrad, 1)
	}
	smolv := mediatorNew[float64, struct{}](*fvib, *fvib, 1)
	magt := make([]float64, nbins)
	first := true
	now := *fhrad
	vib := func(out []float64, in [][]complex128) {
		if *fvib != 0 {
			for range *fvib {
				smolv.Insert(0, struct{}{})
			}
		}
		for w := range nbins {
			mid := .0
			for ch := range in {
				mid += cmplx.Abs(in[ch][w])
			}
			if *fvib != 0 {
				out[w], _ = smolv.Filt(mid, struct{}{})
			} else {
				out[w] = mid
			}
		}
	}
	roto := make2[bool](now*2+1, nbins)
	newroto := make([]int, nbins)
	var field, energy, energy2, ueoai, soutput [][]float64
	var chains [][]int
	if *fseam {
		field = make([][]float64, 0)
		energy = make([][]float64, 0)
		energy2 = make([][]float64, 0)
		ueoai = make([][]float64, 0)
		soutput = make([][]float64, nbins)
		chains = make([][]int, 0)
	}
	common.StftHandle(gr, wavw, fft, [][]float64{win, windowDx(win), windowT(win)}, now, now, func(out [][]complex128, windowed [][][][]complex128, preanalyze bool) {
		frames := windowed[0]
		// Accumulate lookahead.
		if first {
			for f := range frames[now:] {
				vib(magt, frames[f])
				for w := range nbins {
					bank[w].Insert(magt[w], struct{}{}) // Insert two times imitating the mirroring.
					bank[w].Insert(magt[w], struct{}{}) // Order probably won't matter.
				}
			}
			first = false
		}
		vib(magt, frames[len(frames)-1])
		if !preanalyze {
			vib(magt, frames[now])
			vert.Reset(*fvrad)
			for w := range *fvrad { // Mirror the beginning
				mid := magt[w]
				if zerov {
					mid = -mid
				}
				vert.Insert(mid, struct{}{})
			}
			ue := make([]float64, nbins)
			for w := range nbins {
				mid := magt[w]
				if zerov {
					mid = -mid
				}
				v, _ := vert.Filt(mid, struct{}{})
				if zerov {
					v = -v
				}
				h, _ := bank[w].Take()
				if zeroh {
					h = -h
				}
				// lgd := grainwrap(nfft, real(windowed[2][now][0][w]/frames[now][0][w]))
				// disp := int(math.Round(lgd / float64(hop)))
				// m := boolfloat(v > h || roto[min(now+disp, now)][w])
				// _ = disp
				// m := boolfloat(v > h)
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
				// cif := imag(windowed[1][now][0][w]/frames[now][0][w]) / math.Pi / 2
				// lgd := -grainwrap(nfft, -real(windowed[2][now][0][w]/frames[now][0][w])) / float64(hop)
				// // ue[w] = m * math.Log1p(magt[w])
				// gn := min(math.Hypot(lgd, cif), 4)
				// if gn != gn || math.IsInf(gn, 0) {
				// 	gn = 0
				// }
				// ue[w] = gn
				// ue[w] = float64(newroto[w])
			}
			if *fseam {
				field = append(field, slices.Clone(ue))
				energy = append(energy, slices.Clone(ue))
				energy2 = append(energy2, slices.Clone(ue))
				asd := slices.Clone(ue)
				clear(asd)
				ueoai = append(ueoai, asd)
			}

			// oscope.Oscope(ue, oscope.Name(`asdf`))
			roto[len(roto)-1], roto[0] = roto[0], roto[len(roto)-1]
			copy(roto, roto[1:])
		}
		// Vibrato filtering (doesn't work).
		for w := range nbins {
			mid := magt[w]
			if zeroh {
				mid = -mid
			}
			bank[w].Insert(mid, struct{}{})
		}
		// // Include the point if its LGD maps to the percussive frame.
		// // Use finite phase reassignment as a bodge.
		// c := float64(wavr.Properties().Samplerate*nfft) / (4 * math.Pi)
		// (We actually need to save the mask for this)
	})

	if *fseam {
		// for w := range field[0] {
		// 	if w == 0 {
		// 		continue
		// 	}
		// 	for t := range field {
		// 		if t == 0 || t == len(field)-1 {
		// 			continue
		// 		}
		// 		energy[t][w] += min(energy[t-1][w-1], energy[t][w-1], energy[t+1][w-1])
		// 	}
		// }
		// for w := nbins - 1; w >= 0; w-- {
		// 	if w == nbins-1 {
		// 		continue
		// 	}
		// 	for t := range field {
		// 		if t == 0 || t == len(field)-1 {
		// 			continue
		// 		}
		// 		energy2[t][w] += min(energy2[t-1][w+1], energy2[t][w+1], energy2[t+1][w+1])
		// 		energy[t][w] += energy2[t][w]
		// 	}
		// }

		const wdth = 12
		for t := range energy[wdth:] {
			mn, mni := energy[t][0], 0
			ecrop := energy[t : t+wdth]
			lc := make([]int, 0)
			for w := range nbins {
				if t+mni < 0 {
					lc = append(lc, 0)
				} else {
					lc = append(lc, t+mni)
				}
				w2 := 1
				mn = math.MaxFloat64
				for nmni := -w2; nmni <= w2; nmni++ {
					cu := ecrop[max(0, min(nmni, wdth-1))][w]
					if cu <= mn {
						mn = cu
						mni = nmni
					}
				}
			}
			chains = append(chains, lc)
		}

		repeats := make([][]int, nbins)
		for _, ch := range chains {
			for w := range nbins {
				repeats[w] = append(repeats[w], ch[w])
			}
		}
		for w := range repeats {
			// fmt.Println(repeats[w])
			for t := range repeats[w] {
				ueoai[repeats[w][t]][w]++
			}

		}
		for w := range nbins {
			slices.Reverse(repeats[w])
			oscope.Oscope(slices.Clone(repeats[w]), oscope.Name(`repeats`))
		}
		// sinput := make([][]float64, nbins)
		shift := make([]int, nbins)
		for w := range nbins {
			for t := range field {
				soutput[w] = append(soutput[w], field[t][w])
				for len(repeats[w]) > 0 && repeats[w][len(repeats[w])-1] <= t+shift[w] {
					shift[w]++
					soutput[w] = append(soutput[w], 0)
					// ueoai[t][w]++
					repeats[w] = repeats[w][:len(repeats[w])-1]
				}
			}
		}

		for _, f := range soutput {
			d := len(f) - 3000
			if d <= 0 {
				// b := len(f)
				f = append(f, make([]float64, -d)...)
				// println(b, len(f))
			}
			oscope.Oscope(f, oscope.Name(`soutput`))
		}
		// for _, f := range soutput {
		// 	m0 = min(m0, len(f))
		// }
		for _, f := range energy {
			oscope.Oscope(f, oscope.Name(`energy`))
		}
		for _, f := range ueoai {
			oscope.Oscope(f, oscope.Name(`ueoai`))
		}
	}

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
