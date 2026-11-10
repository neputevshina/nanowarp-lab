package main

import (
	"flag"
	"math"
	"os"

	"github.com/neputevshina/nanowarp-lab/common"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/dspio/wavio"
	"github.com/neputevshina/nanowarp/oscope"
	"github.com/neputevshina/nanowarp/pffft"
	"golang.org/x/exp/constraints"
)

var finputa = flag.String("a", "", "source WAV")
var folap = flag.Int("l", 6, "overlap")
var fnfft = flag.Int("nfft", 3840, "fft grain size")

func main() {
	flag.Parse()
	oscope.Enable = true

	fa, err := os.Open(*finputa)
	if err != nil {
		panic(err)
	}

	nfft := *fnfft
	hop := nfft / *folap
	nbins := nfft/2 + 1

	wavr, err := wavio.NewDecoder(fa)
	if err != nil {
		panic(err)
	}
	p := wavr.Properties()
	_ = nbins
	_ = p

	gr := dspio.NewGrainReader(nfft, hop, wavr)
	fft := pffft.New(nfft)

	// HPSS works best with windows with the smallest main lobe width,
	// but a window's derivative must have compact support (zeros only before and after the window)
	// to be perfectly invertible after modifications.
	// Meaning, a rectangular or a Tukey window will produce artifacts.
	win := make([]float64, nfft)
	hann(win)

	// var field, energy, energy2, ueoai, soutput [][]float64
	// var chains [][]int
	common.StftHandle(gr, nil, fft, [][]float64{win, windowDx(win), windowT(win)}, 0, 0, func(out [][]complex128, windowed [][][][]complex128, preanalyze bool) {
		// frames := windowed[0]
		// Accumulate lookahead.
		return
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
