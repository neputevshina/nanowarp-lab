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
var fhrad = flag.Int("h", 6, "horizontal radius")
var fvrad = flag.Int("v", 6, "vertical radius")

var fvib = flag.Int("vib", 0, "vertical vibrato radius")
var ftonal = flag.Bool("t", false, "output tonals")
var fhperc = flag.Int("ph", 50, "horizontal percentile")
var fvperc = flag.Int("pv", 50, "vertical percentile")

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
	common.StftHandle(gr, wavw, fft, win, now, now, func(out [][]complex128, frames [][][]complex128, preanalyze bool) {
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
				m := boolfloat(v > h)
				if *ftonal {
					m = 1 - m
				}
				for ch := range out {
					out[ch][w] = frames[now][ch][w] * complex(m, 0)
				}
			}
		}
		for w := range nbins {
			mid := magt[w]
			if zeroh {
				mid = -mid
			}
			bank[w].Insert(mid, struct{}{})
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
