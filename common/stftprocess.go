package common

import (
	"io"
	"slices"

	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/pffft"
	"gonum.org/v1/gonum/floats"
)

// StftHandle automatically operates a time-uniform short-time Fourier transform pipleine.
//
// It reads grains from the signal, applies lookahead of ‘future’ number of frames,
// applies windowing, gain normalization, performs FFT on the signal and gives it for consumer to process it.
// The transform can be oversampled/downsampled if rd.N() ≠ fft.Len().
// Synthesis windowing is applied with automatically generated dual of the given window function.
// The size of frames in consumer is future+past+1, where the current frame (with the same number as ‘out’)
// is given by frames[past].
//
// It does not try to write to wr if it is nil. Changing gr.Hop while this routine is working will
// produce broken synthesis result.
func StftHandle(gr *dspio.GrainReader, wr dspio.SignalWriter,
	fft pffft.Fourier, window []float64, future int, past int,
	consumer func(out [][]complex128, frames [][][]complex128, preanalyze bool)) error {
	nch := gr.NchRead()
	if nch != wr.NchWrite() {
		panic(`StftHandle: different number of channels between reader and writer`)
	}

	dual, norm := windowDualUniform(window, gr.Hop)
	gw := dspio.NewGrainWriter(gr.N(), gr.Hop, wr)
	grain := make2[float64](gr.NchRead(), gr.N())
	makebins := func() [][]complex128 { return make2[complex128](gr.NchRead(), fft.Len()/2+1) }
	out := makebins()
	frames := make([][][]complex128, future+past+1)
	for n := range frames {
		frames[n] = makebins()
	}
	last := future + past
	countdown := future
	countup := -1

	nfft := fft.Len()
	for {
		if countup >= 0 {
			countup++
		} else {
			_, rerr := gr.SignalRead(nil, grain)
			if rerr == io.EOF {
				countup = 0 // Prepare to close.
			} else if rerr != nil {
				return rerr
			}
		}
		frames[last], frames[0] = frames[0], frames[last]
		copy(frames, frames[1:])
		for ch := range nch {
			floats.Mul(grain[ch], window)
			fft.Coefficients(frames[last][ch], grain[ch])
		}
		consumer(out, frames, countdown > 0)
		if countdown > 0 { // Accumulate lookahead.
			countdown--
			continue
		}
		for ch := range nch {
			fft.Sequence(grain[ch], out[ch])
			floats.Scale(1/(float64(nfft)*float64(gr.Hop)*norm), grain[ch])
			floats.Mul(grain[ch], dual)
		}
		_, werr := gw.SignalWrite(nil, grain)
		if werr != nil {
			return werr
		}
		if countup == future {
			return nil
		}
	}
}

func windowDualUniform(w []float64, hop int) (dual []float64, gain float64) {
	dual = make([]float64, len(w))
	buf := make([]float64, 4*len(w))
	w2 := slices.Clone(w)
	floats.Mul(w2, w)
	half := len(w) / 2
	for i := half; i < len(buf)-half; i += hop {
		floats.Add(buf[i-half:i+half], w2)
	}
	copy(dual, w)
	floats.Mul(dual, buf[len(w):2*len(w)])
	cl := slices.Clone(w)
	floats.Mul(cl, dual)
	gain = floats.Sum(cl) / float64(len(w))
	return
}

func make2[T any](j, i int) (v [][]T) {
	v = make([][]T, j)
	for j := range j {
		v[j] = make([]T, i)
	}
	return
}
