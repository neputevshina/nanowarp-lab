package main

import (
	"flag"
	"math/bits"
	"os"
	"slices"

	"github.com/neputevshina/nanowarp-lab/common"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/dspio/wavio"
	"github.com/neputevshina/nanowarp/pffft"
	"gonum.org/v1/gonum/cmplxs"
)

var finputa = flag.String("a", "", "source a WAV")

func main() {
	flag.Parse()

	fa, err := os.Open(*finputa)
	if err != nil {
		panic(err)
	}

	nfft, hop := 3600, 600
	nbins := nfft/2 + 1

	wavr, err := wavio.NewDecoder(fa)
	p := wavr.Properties()
	wavw, err := wavio.NewEncoder(fa, p.Samplerate, p.Nch, p.Format, 32)

	gr := dspio.NewGrainReader(nfft, hop, wavr)
	fft := pffft.New(nfft)

	heap := make(hp, nbins)
	arm := make([]bool, nbins)
	mag := make([]float64, nbins)
	ridges := make2[uint64](p.Nch, nbins)
	trace := make3[float64](p.Nch, 7, nbins)
	common.StftHandle(gr, wavw, fft, kaiser(nfft, 3), 7, 0, func(out [][]complex128, frames [][][]complex128, preanalyze bool) {
		for ch := range out {
			cmplxs.Abs(mag, frames[len(frames)-1][ch])
			for w := range mag {
				heap[w] = heaptriple{mag: mag[w], w: w, t: -1}
				ridges[ch][w] <<= bits.OnesCount(ridgemask)
			}
			fill(arm, true)
			heap = heap[:nbins]
			heapInit(&heap)

			hood := [][]uint64{
				{up, 0},
				{0, right},
				{down, 0},
			}
			for len(heap) > 0 {
				h := heapPop(&heap)
				for y, t := range hood {
					for x := range t {
						if hood[y][x] > 0 && arm[h.w+y-1] {
							ridges[ch][h.w] |= hood[y][x]
							arm[h.w+y-1] = false
						}
					}
				}
			}

			rotate(trace[ch])
		}
	})

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
	for w, v := range ridges {
		p := boolfloat(bits.OnesCount(v&(ridgemask<<3)) >= 2)
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
			if trace[e] == 0 && ridges[e]&(down<<3) > 0 {
				out[e] = trace[w]
			} else {
				break
			}
		}
		for e := w + 1; e < len(trace) && e-w <= InfluenceRadius; e++ {
			if trace[e] == 0 && ridges[e]&(up<<3) > 0 {
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
	v = make([][][]T, j)
	for k := range k {
		v[k] = make([][]T, i)
		for j := range j {
			v[k][j] = make([]T, i)
		}
	}
	return
}

func rotate[T any](frames []T) {
	last := len(frames) - 1
	frames[last], frames[0] = frames[0], frames[last]
	copy(frames, frames[1:])
}
