package main

import (
	"flag"
	"math"
	"math/bits"
	"os"
	"slices"

	"github.com/neputevshina/nanowarp-lab/common"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/dspio/wavio"
	"github.com/neputevshina/nanowarp/oscope"
	"github.com/neputevshina/nanowarp/pffft"
	"golang.org/x/exp/constraints"
	"gonum.org/v1/gonum/cmplxs"
)

var finputa = flag.String("a", "", "source a WAV")
var ftonal = flag.Bool("t", false, "output tonals")

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
		trlen  = 6
		height = 1
		infl   = 5
	)

	heap := make(hp, nbins)
	arm := make([]bool, nbins)
	mag := make([]float64, nbins)
	pmags := make2[float64](p.Nch, nbins)
	ridges := make2[uint](p.Nch, nbins)
	trace := make3[float64](p.Nch, trlen, nbins)
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
			rotate(trace[ch])

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

			copy(trace[ch][trlen-1], trace[ch][trlen-2])
			trackridges(trace[ch][trlen-1], ridges[ch], height)

			backprop(trace[ch])
			copy(pmags[ch], mag)
		}

		if !preanalyze {
			for ch := range out {
				fatten(traceaccum[ch], trace[ch][1], ridges[ch], infl)
				for w := range mag {
					c := traceaccum[ch][w]
					what := c < trlen
					if *ftonal {
						what = c >= trlen
					}
					out[ch][w] = frames[0][ch][w] * complex(boolfloat(what), 0)
				}
			}
		}

		// oscope.Oscope(slices.Clone(traceaccum[0]), oscope.Name(`accum`))
		// oscope.Oscope(slices.Clone(trace[0][0]), oscope.Name(`trace`))
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

func trackridges(trace []float64, ridges []uint, octheight float64) {
	rl := bits.OnesCount(ridgemask)
	for w, v := range ridges {
		p := boolfloat(bits.OnesCount(v&(ridgemask<<rl)) >= 2)
		trace[w] = trace[w]*p + p
	}

	// Propagate vertically.
	l := -1
	height := func(oct float64, i int) int {
		log := math.Log2
		type f = float64
		v := oct * (log(f(i)) - log(1))
		return int(max(3, v))
	}
	for i := range trace {
		if l < 0 && trace[i] != 0 {
			l = i
		}
		if l >= 0 && trace[i] == 0 {
			v := slices.Max(trace[l:i])
			// Reset the track on a PGHI-detected transient.
			if i-l >= height(octheight, i) {
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
}

func fatten(out, trace []float64, ridges []uint, InfluenceRadius int) {
	rl := bits.OnesCount(ridgemask)
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
}

func backprop(traces [][]float64) {
	for t := len(traces) - 2; t >= 0; t-- {
		for w := range traces[t] {
			if traces[t][w] > 0 {
				traces[t][w] = max(traces[t][w], traces[t+1][w])
			}
		}
		for w := 1; w < len(traces[t]); w++ {
			if traces[t][w] > 0 {
				traces[t][w] = max(traces[t][w-1], traces[t][w])
			}
		}
		for w := len(traces[t]) - 1; w <= 0; w-- {
			if traces[t][w] > 0 {
				traces[t][w] = max(traces[t][w], traces[t][w+1])
			}
		}
	}
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

func unmix[F constraints.Float](a, b, x F) F {
	return (x - a) / (b - a)
}

func mix[F constraints.Float](a, b, x F) F {
	return a*(1-x) + b*x
}
