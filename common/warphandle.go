package common

import (
	"io"

	"github.com/neputevshina/nanowarp"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/pffft"
	"gonum.org/v1/gonum/floats"
)

// WarpHandle automatically operates a non-time-uniform short-time Fourier transform pipleine.
//
// Same as [StftHandle], but all windows are applied and producer can push a frame (which will be inverted
// with dual of windows[0]) anytime. Pushing frames does not work while preanalyze is true.
//
// frames of producer is addressed like [window][frame][channel][bin]complex128.
//
// It does not try to write to wr if it is nil. Changing gr.Hop while this routine is working will
// produce broken synthesis result.
func WarpHandle(gr dspio.GrainSeeker, wr dspio.SignalWriter, phasor *nanowarp.Phasor,
	fft pffft.Fourier, grainsize, synthhop int, windows [][]float64, future int, past int,
	producer func(push func(out [][]complex128), frames [][][][]complex128, stretch float64, preanalyze bool)) error {
	nch := gr.NchRead()
	if nch != wr.NchWrite() {
		panic(`WarpHandle: different number of channels between reader and writer`)
	}

	dual, norm := windowDualUniform(windows[0], synthhop)
	gw := dspio.NewGrainWriter(grainsize, synthhop, wr)
	grain := make2[float64](gr.NchRead(), grainsize)

	makebins := func() [][]complex128 { return make2[complex128](gr.NchRead(), fft.Len()/2+1) }
	frames := make([][][][]complex128, len(windows))
	for w := range frames {
		for n := range frames[w] {
			frames[w][n] = makebins()
		}
	}
	last := future + past
	countdown := future
	countup := -1

	var (
		lead [][]float64
		rerr error
	)

	nfft := fft.Len()
	for j := -grainsize / 2; ; j += synthhop {
		i := int(phasor.ReverseSample(float64(j)))
		c := 1 / phasor.Dy(float64(j)) // Stretch, inverse of scan speed, which Dy is.

		if countup >= 0 {
			countup++
		} else {
			lead, rerr = gr.GrainSeek(nil, int64(i-grainsize/2), grainsize)
			if rerr == io.EOF {
				countup = 0 // Prepare to close.
			} else if rerr != nil {
				return rerr
			}
		}
		for wn, fr := range frames {
			fr[last], fr[0] = fr[0], fr[last]
			copy(fr, fr[1:])
			for ch := range nch {
				floats.Mul(lead[ch], windows[wn])
				fft.Coefficients(frames[wn][last][ch], lead[ch])
			}
		}
		var werr error
		producer(func(out [][]complex128) {
			if werr != nil || countdown > 0 {
				return
			}
			for ch := range nch {
				fft.Sequence(grain[ch], out[ch])
				floats.Scale(1/(float64(nfft)*float64(nfft)/float64(synthhop)*norm), grain[ch])
				floats.Mul(grain[ch], dual)
			}
			_, werr = gw.SignalWrite(nil, grain)
		}, frames, c, countdown > 0)
		if werr != nil {
			return werr
		}
		if countdown > 0 { // Accumulate lookahead.
			countdown--
			continue
		}
		if countup == future {
			return nil
		}
	}
}
