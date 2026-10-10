package main

import (
	"flag"
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
var fhrad = flag.Int("h", 31, "horizontal radius")
var fvrad = flag.Int("v", 31, "vertical radius")
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

	// HPSS works best with a rectangular window.
	// Probably because of unit main lobe width.
	win := make([]float64, nfft)
	fill(win, 1)

	bank := make([]*mediator[float64, struct{}], nbins)
	for i := range bank {
		bank[i] = mediatorNew[float64, struct{}](*fhrad, *fhrad, float64(*fhperc)/100)
	}
	vert := mediatorNew[float64, struct{}](*fvrad, *fvrad, float64(*fvperc)/100)
	first := true
	mid := make([]float64, nbins)
	now := *fhrad
	common.StftHandle(gr, wavw, fft, win, now, now, func(out [][]complex128, frames [][][]complex128, preanalyze bool) {
		// Accumulate lookahead.
		if first {
			for f := range frames[now:] {
				clear(mid)
				for w := range nbins {
					for ch := range out {
						mid[w] += cmplx.Abs(frames[f][ch][w])
					}
					bank[w].Insert(mid[w], struct{}{}) // Insert two times imitating the mirroring.
					bank[w].Insert(mid[w], struct{}{}) // Order probably won't matter.
				}
			}
			first = false
		}
		clear(mid)
		for w := range nbins {
			for ch := range out {
				mid[w] += cmplx.Abs(frames[len(frames)-1][ch][w])
			}
			bank[w].Insert(mid[w], struct{}{})
		}
		if !preanalyze {
			clear(mid)
			for w := range nbins {
				for ch := range out {
					mid[w] += cmplx.Abs(frames[now][ch][w])
				}
				v, _ := vert.Filt(mid[w], struct{}{})
				h, _ := bank[w].Take()
				m := boolfloat(v > h)
				if *ftonal {
					m = 1 - m
				}
				for ch := range out {
					out[ch][w] = frames[now][ch][w] * complex(m, 0)
				}
			}
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
