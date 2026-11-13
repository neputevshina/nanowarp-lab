package main

import (
	"flag"
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"slices"

	"github.com/neputevshina/nanowarp-lab/common"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/dspio/wavio"
	"github.com/neputevshina/nanowarp/oscope"
	"github.com/neputevshina/nanowarp/pffft"
	"golang.org/x/exp/constraints"
)

var finputa = flag.String("a", "", "source WAV")
var folap = flag.Int("l", 6, "overlap")
var fnfft = flag.Int("nfft", 4096, "fft grain size")
var fnf = flag.Float64("nf", -90, "noise floor in `db`")
var ftrim = flag.Float64("trim", 0, "trim in `db`")
var fhpss = flag.Bool("hpss", true, "enable hpss")

func main() {
	flag.Parse()
	oscope.Enable = true

	if *finputa == `` {
		flag.PrintDefaults()
		return
	}

	fa, err := os.Open(*finputa)
	if err != nil {
		panic(err)
	}

	nfft := *fnfft
	hop := nfft / *folap
	nbins := nfft/2 + 1

	wavr, err := wavio.NewDecoder(fa)
	if err != nil {
		// Try to call ffmpeg, we've probably got an MP3.
		if err == wavio.ErrNotAWav {
			_ = fa.Close()

			_, err = exec.LookPath(`ffmpeg`)
			if err != nil {
				fmt.Fprintln(os.Stderr, `can process only WAV files without ffmpeg`)
				os.Exit(1)
			}

			s2, err := filepath.Abs(*finputa)
			if err != nil {
				panic(err)
			}

			ex := path.Join(os.TempDir(), path.Base(s2)) + `.wav`

			cmd := exec.Command(`ffmpeg`, `-v`, `warning`, `-y`, `-i`, s2, `-acodec`, `pcm_f32le`, ex)
			cmd.Stderr = os.Stderr
			cmd.Stdout = os.Stdout
			err = cmd.Run()
			if err != nil {
				if _, ok := err.(*exec.ExitError); ok {
					os.Exit(1)
				} else {
					panic(err)
				}
			}

			fa, err = os.Open(ex)
			if err != nil {
				panic(err)
			}
			wavr, err = wavio.NewDecoder(fa)
			if err != nil {
				panic(err)
			}
		}
	}
	p := wavr.Properties()

	gr := dspio.NewGrainReader(nfft, hop, wavr)
	pr, pw := dspio.GoPipe(p.Nch)
	gr2 := dspio.NewGrainReader(nfft, hop, dspio.LimitReader(pr, p.Samples))

	win := make([]float64, nfft)
	hann(win)

	go func() { // Process 1: HPSS tonal extraction.
		fft := pffft.New(nfft)
		bank := make([]*mediator[float64, struct{}], nbins)
		hrad, vrad := *folap, 23
		now := hrad
		for i := range bank {
			bank[i] = mediatorNew[float64, struct{}](hrad, hrad, 0.5)
		}
		vert := mediatorNew[float64, struct{}](vrad, vrad, 0.5)
		magt := make([]float64, nbins)
		first := true
		roto := make2[bool](now*2+1, nbins)
		newroto := make([]int, nbins)
		vib := func(out []float64, in [][]complex128) {
			for w := range nbins {
				mid := .0
				for ch := range in {
					mid += cmplx.Abs(in[ch][w])
				}
				out[w] = mid
			}
		}
		common.StftHandle(gr, pw, fft, [][]float64{win, windowDx(win), windowT(win)}, now, now, func(out [][]complex128, windowed [][][][]complex128, preanalyze bool) {
			frames := windowed[0]
			if first { // Accumulate lookahead.
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

				vert.Reset(vrad)
				for w := range vrad { // Mirror the beginning
					mid := magt[w]
					vert.Insert(mid, struct{}{})
				}

				for w := range nbins {
					mid := magt[w]
					v, _ := vert.Filt(mid, struct{}{})
					h, _ := bank[w].Take()
					roto[now][w] = v > h
					newroto[w] = int(boolfloat(roto[now][w]))
				}

				for w := range nbins {
					m := boolfloat(newroto[w] <= 0) // Extract tonals
					for ch := range out {
						out[ch][w] = frames[now][ch][w] * complex(m, 0)
					}
				}
				roto[len(roto)-1], roto[0] = roto[0], roto[len(roto)-1]
				copy(roto, roto[1:])
			}
			for w := range nbins {
				mid := magt[w]
				bank[w].Insert(mid, struct{}{})
			}
		})
	}()

	// Process 2: reassigned chromagram accumulation.
	fft := pffft.New(nfft)
	tex := make2[float64](600, 200)
	mx := 0.0
	cut := int(hztobin(2500, nfft, float64(p.Samplerate)))
	_, c0 := math.Modf(math.Log2(hztobin(32.703, nfft, float64(p.Samplerate))))
	common.StftHandle(gr2, nil, fft, [][]float64{win, windowDx(win), windowT(win)}, 0, 0, func(out [][]complex128, windowed [][][][]complex128, preanalyze bool) {

		for ch := range p.Nch {
			for w := range cut {
				abs := cmplx.Abs(windowed[0][0][ch][w])
				if abs == 0 {
					continue
				}
				cif := imag(windowed[1][0][ch][w]/windowed[0][0][ch][w]) / math.Pi / 2
				_, x := math.Modf(math.Log2(float64(w)+cif) - c0)
				if x != bitsafe(x) {
					continue
				}
				if x < 0 {
					x++
				}
				y := clamp(0, 1, unmix(*fnf, *ftrim, atodb(abs)-60))
				px := int(math.Round(mix(0, float64(len(tex)-1), x)))
				py := int(math.Round(mix(float64(len(tex[0])-1), 0, y)))
				if py == len(tex[0])-1 {
					continue
				}
				tex[px][py] += 1 / math.Sqrt(float64(w+1))
				mx = max(mx, tex[px][py])
			}
		}
	})

	// Output.
	d := -1.
	for n := range tex {
		a, _ := math.Modf(unmix(0, float64(len(tex)), float64(n)) * 12)
		if a != d {
			tex[n][len(tex[0])-1] = mx
		}
		d = a
	}
	for _, l := range tex {
		oscope.Oscope(slices.Clone(l), oscope.Name(path.Base(*finputa)))
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

// dbtoa converts a value in decibels full scale (dBFS) to amplitude.
func dbtoa(db float64) float64 {
	return math.Pow(10, db/20)
}

// atodb converts an amplidude value to decibels full scale (dBFS).
func atodb(a float64) float64 {
	return 20 * math.Log10(a)
}

func bitsafe(v float64) float64 {
	if v != v || math.IsInf(v, 0) {
		return 0
	}
	return v
}

// hztobin converts a frequency value in hertz to bin number of a DFT
// of size nfft on a signal with known sample rate.
func hztobin(hz float64, nfft int, samplerate float64) float64 {
	return hz * float64(nfft) / samplerate
}
