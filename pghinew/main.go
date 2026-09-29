package main

import (
	"flag"
	"math/bits"
	"os"
	"slices"

	"github.com/neputevshina/nanowarp-lab/common"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/dspio/wavio"
	"github.com/neputevshina/nanowarp/oscope"
	"github.com/neputevshina/nanowarp/pffft"
	"gonum.org/v1/gonum/cmplxs"
)

var finputa = flag.String("a", "", "source a WAV")

func main() {
	flag.Parse()
	oscope.Enable = true

	fa, err := os.Open(*finputa)
	if err != nil {
		panic(err)
	}
	fo, err := os.Create("./output.wav")
	if err != nil {
		panic(err)
	}

	nfft, hop := 3840, 640
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

	gr := dspio.NewGrainReader(nfft, hop, wavr)
	fft := pffft.New(nfft)

	const (
		trlen  = 12
		height = 5
		infl   = 3
	)

	heap := make(hp, nbins)
	arm := make([]bool, nbins)
	mag := make([]float64, nbins)
	pmags := make2[float64](p.Nch, nbins)
	ridges := make2[uint](p.Nch, nbins)
	trace := make3[float64](p.Nch, trlen+1, nbins)
	traceaccum := make2[float64](p.Nch, nbins)
	// armget := func(t, w int) bool {
	// 	if t != 0 {
	// 		return false
	// 	}
	// 	if w < 0 || w >= len(arm) {
	// 		return false
	// 	}
	// 	return arm[w]
	// }
	// armset := func(w int, v bool) {
	// 	if w < 0 || w >= len(arm) {
	// 		return
	// 	}
	// 	arm[w] = v
	// }
	rl := bits.OnesCount(ridgemask)
	common.StftHandle(gr, wavw, fft, kaiser(nfft, 3), trlen, 0, func(out [][]complex128, frames [][][]complex128, preanalyze bool) {
		for ch := range out {
			heap = heap[:nbins]
			cmplxs.Abs(mag, frames[len(frames)-1][ch])
			for w := range mag {
				heap[w] = heaptriple{mag: pmags[ch][w], w: w, t: -1}
				ridges[ch][w] <<= rl
			}
			fill(arm, true)
			heapInit(&heap)

			// hood := [][]uint{
			// 	{up, 0},
			// 	{0, right},
			// 	{down, 0},
			// }
			for len(heap) > 0 {
				h := heapPop(&heap)
				w := h.w
				switch h.t {
				case -1:
					if arm[w] {
						ridges[ch][w] |= right << rl
						arm[w] = false
						heapPush(&heap, heaptriple{mag[w], w, 0})
					}
				case 0:
					if w >= 1 && arm[w-1] {
						ridges[ch][w] |= down
						arm[w-1] = false
						heapPush(&heap, heaptriple{mag[w-1], w - 1, 0})
					}
					if w < nbins-1 && arm[w+1] {
						ridges[ch][w] |= up
						arm[w+1] = false
						heapPush(&heap, heaptriple{mag[w+1], w + 1, 0})
					}
				}
				// for y, t := range hood {
				// 	for x := range t {
				// 		m := h.w + y - 1
				// 		f := h.t + x
				// 		if hood[y][x] > 0 && armget(f, m) {
				// 			ridges[ch][h.w] |= hood[y][x] << (rl * x)
				// 			armset(m, false)
				// 			heapPush(&heap, heaptriple{mag: mag[m], w: m, t: 0})
				// 		}
				// 	}
				// }
			}
			trace[ch][len(trace)-1] = trackridges(trace[ch][len(trace)-1], traceaccum[ch], ridges[ch], height, infl)

			// for w := range mag {
			// 	trace[ch][len(trace)-1][w] = boolfloat(traceaccum[w] >= 1)
			// }
			rotate(trace[ch])
			copy(pmags[ch], mag)
		}
		oscope.Oscope(slices.Clone(traceaccum[0]), oscope.Name(`accum`))
		oscope.Oscope(slices.Clone(trace[0][0]), oscope.Name(`trace`))

		if !preanalyze {
			for ch := range out {
				for w := range mag {
					out[ch][w] = frames[0][ch][w] * complex(boolfloat(trace[1][ch][w] >= trlen), 0)
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

func trackridges(out, trace []float64, ridges []uint, HighRidgeHeight, InfluenceRadius int) []float64 {
	rl := bits.OnesCount(ridgemask)
	for w, v := range ridges {
		p := boolfloat(bits.OnesCount(v&(ridgemask<<rl)) >= 2)
		trace[w] = trace[w]*p + p
	}

	// Propagate vertically.
	l := -1
	for i := range trace {
		if l < 0 && trace[i] != 0 {
			l = i
		}
		if l >= 0 && trace[i] == 0 {
			v := slices.Max(trace[l:i])
			// Reset the track on a PGHI-detected transient.
			if i-l >= HighRidgeHeight {
				v = 0
			}
			fill(trace[l:i], v)
		}
		if trace[i] == 0 {
			l = -1
		}
	}
	if l > 0 {
		fill(trace[l:], slices.Max(trace[l:]))
	}
	// Propagate each trace to its native (per PGHI directions) region of influence,
	// limited by InfluenceRadius hyperparameter.
	clear(out)
	for w, v := range trace {
		if v == 0 {
			continue
		}
		for e := w - 1; e >= 0 && w-e <= InfluenceRadius; e-- {
			if trace[e] == 0 && ridges[e]&(down<<rl) > 0 {
				out[e] = trace[w]
			} else {
				break
			}
		}
		for e := w + 1; e < len(trace) && e-w <= InfluenceRadius; e++ {
			if trace[e] == 0 && ridges[e]&(up<<rl) > 0 {
				out[e] = trace[w]
			} else {
				break
			}
		}
	}
	// Add original traces to the output.
	for w, v := range trace {
		if v == 0 {
			continue
		}
		out[w] = v
	}
	return out
}

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
