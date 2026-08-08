package main

import (
	"io"
	"math"
	"math/cmplx"
	"slices"

	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/oscope"
	"gonum.org/v1/gonum/dsp/fourier"
	"gonum.org/v1/gonum/floats"
)

type warper struct {
	nfft  int     // DFT size, a power of 2
	nbuf  int     // Effective window size, nbuf<nfft
	hop   int     // Window output hop size
	lah   int     // Non-causal PGHI lookahead in frames (hop sizes)
	nbins int     // nfft/2+1, Number of DFT bins
	olap  int     // nbuf/hop, Window ovelap
	osamp float64 // nfft/nbuf, Zero-padding ratio

	fft         *fourier.FFT
	arm         [][]bool // PGHI done mask
	norm, wgain float64  // Global normalization factor and grain-only normalization factor
	heap        hp       // PGHI heap

	a wbufs
}

type wbufs struct {
	Mid, S                     []float64      // Scratch buffers
	Ph, M, P, F                []float64      `size:"nbins"` // Current phase
	Past, Future               []float64      // Phase accumulators
	Fadv, Tadv                 []float64      `size:"nbins"`
	W, Wr, Wd, Wt, Wdt         []float64      // Window functions
	X, Y, Xd, Xt, L, R, Lo, Ro []complex128   // Complex spectra
	C, Co                      [][]complex128 // Channels

	Cs                    [][][]complex128
	Phs, Fadvs, Tadvs, Ms [][]float64 `size:"lah"`
}

func (w *warper) process(ain dspio.SignalReader) error {
	agr := dspio.NewGrainReader(w.nbuf, w.hop, ain)

	agrain := make2(ain.NchRead(), w.nbuf)
	var err error
	for {
		_, err = agr.SignalRead(nil, agrain)
		if err != nil {
			break
		}
		w.advance(agrain)
	}
	if err == io.EOF {
		err = nil
	}
	return nil
}

func warperNew(nbuf, osamp, olap, nch int) (n *warper) {
	// FIXME Only 2x oversampling works, no more, no less.
	nfft := nextpow2(nbuf * osamp)
	n = &warper{
		nfft:  nfft,
		nbins: nfft/2 + 1,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		olap:  olap,
		osamp: float64(osamp),
	}
	a := &n.a

	makeslices(a, n.nbins, nfft, nch, n.lah)

	n.arm = make([][]bool, n.lah)
	for i := range n.arm {
		n.arm[i] = make([]bool, n.nbins)
	}

	s := func(w []float64) []float64 {
		// return w[nfft/2-nbuf/2 : nfft/2+nbuf/2]
		return w[:nbuf]
	}
	blackmanHarris(s(a.W))

	// FIXME Destination is the first argument by convention.
	windowDx(s(a.W), s(a.Wd))
	windowT(s(a.W), s(a.Wt))
	windowT(s(a.Wd), s(a.Wdt))

	copy(s(a.Wr), s(a.W))
	slices.Reverse(s(a.Wr))
	n.wgain = windowGain(n.a.W)
	n.norm = float64(nfft) * float64(n.olap) * n.osamp * n.wgain

	n.fft = fourier.NewFFT(nfft)
	n.heap = make(hp, n.lah*n.nbins) // 2 for future and past.

	return
}

func (n *warper) enfft(x []complex128, w, grain []float64) {
	a := &n.a
	clear(a.S)
	copy(a.S, grain)
	if w != nil {
		mul(a.S, w)
	}
	n.fft.Coefficients(x, a.S)
}

func (n *warper) defft(out []float64, x []complex128, w bool) {
	a := &n.a
	n.fft.Sequence(a.S, x)
	floats.Scale(1./n.norm, a.S)
	if w {
		mul(a.S, a.Wr)
	}
	copy(out, a.S)
}

func (n *warper) resetPast(present [][]float64) {
	a := &n.a
	clear(a.Mid)
	for ch := range present {
		floats.Add(a.Mid, present[ch])
		n.enfft(a.C[ch], a.W, present[ch])
	}
	n.enfft(a.X, a.W, a.Mid)
	for w := range a.Ph {
		a.P[w], a.Past[w] = cmplx.Polar(a.X[w])
	}
}

func (n *warper) bash(present [][]float64, C [][]complex128, M, Ph []float64, Mid []float64) {
	a := &n.a
	clear(Mid)
	for ch := range present {
		floats.Add(Mid, present[ch])
		n.enfft(C[ch], a.W, present[ch])
	}
	n.enfft(a.X, a.W, Mid)
	for w := range n.nbins {
		M[w], Ph[w] = cmplx.Polar(a.X[w])
	}
}

func (n *warper) analyze(present [][]float64, C [][]complex128, Fadv, Tadv, M, Mid []float64) {
	a := &n.a

	clear(Mid)
	for ch := range present {
		floats.Add(Mid, present[ch])
		n.enfft(C[ch], a.W, present[ch])
	}

	n.enfft(a.X, a.W, Mid)
	n.enfft(a.Xd, a.Wd, Mid)
	n.enfft(a.Xt, a.Wt, Mid)
	n.enfft(a.Y, a.Wdt, Mid)

	for w := range a.X {
		Fadv[w] = princarg(getfadv(a.X, a.Xt, 2./n.osamp)(w))
		Tadv[w] = tdx(a.X, a.Xd, w, float64(n.olap)*n.osamp)
	}

	for w := range a.X {
		m := mag(a.X[w])
		p := a.X[w] / complex(m, 0)
		if m < 1e-6 {
			p = complex(1, 0)
		}
		M[w] = m
		a.Y[w] = p
	}
	for ch := range present {
		for w := range a.X {
			C[ch][w] /= a.Y[w]
		}
	}
}

func (n *warper) integrate(Fadv, Tadv, M [][]float64, Ph [][]float64, arm [][]bool) {
	n.heap = n.heap[:0]
	clear(n.heap)

	// Frames on sides of the framebuffer are the past frame after the
	// previous transient, and the future frame, which is the first
	// frame of a transient.
	//
	// They are ground truths for the current step of integration, so they are
	// added to the heap and not armed for phase reconstruction.
	for t := range Ph {
		if t > 0 && t < n.lah {
			fill(arm[t], true)
		} else {
			fill(arm[t], false)
			for w := range n.nbins {
				n.heap = append(n.heap, heaptriple{M[t][w], w, t})
			}
		}
	}
	heapInit(&n.heap)

	// Do PGHI.
	for len(n.heap) > 0 {
		h := heapPop(&n.heap)
		w, t := h.w, h.t
		if t >= 1 && arm[t-1][w] {
			Ph[t-1][w] = Ph[t][w] - Tadv[t-1][w]
			arm[t-1][w] = false
			heapPush(&n.heap, heaptriple{M[t-1][w], w, t - 1})
		}
		if t < len(Ph)-1 && arm[t+1][w] {
			Ph[t+1][w] = Ph[t][w] + Tadv[t+1][w]
			arm[t+1][w] = false
			heapPush(&n.heap, heaptriple{M[t+1][w], w, t + 1})
		}
		if w >= 1 && arm[t][w-1] {
			Ph[t][w-1] = Ph[t][w] - Fadv[t][w-1]
			arm[t][w-1] = false
			heapPush(&n.heap, heaptriple{M[t][w-1], w - 1, t})
		}
		if w < n.nbins-1 && arm[t][w+1] {
			Ph[t][w+1] = Ph[t][w] + Fadv[t][w+1]
			arm[t][w+1] = false
			heapPush(&n.heap, heaptriple{M[t][w+1], w + 1, t})
		}
	}
}

func (n *warper) advance(presenta [][]float64) {
	a := &n.a

	n.analyze(presenta, a.C, a.Fadv, a.Tadv, a.M, a.Mid)

	oscope.Oscope(slices.Clone(a.Fadv), oscope.Name(`y`))
	oscope.Oscope(slices.Clone(a.Tadv), oscope.Name(`x`))

	// n.integrate([][]float64{nil, a.Fadv}, [][]float64{nil, a.Tadv}, [][]float64{a.P, a.M}, [][]float64{a.Past, a.Ph}, n.arm)

	// for w := range a.Ph {
	// 	a.Past[w] = princarg(a.Ph[w])
	// }
	// copy(a.P, a.M)

}

func getfadv(x, xt []complex128, stretch float64) func(w int) float64 {
	return func(j int) float64 {
		if mag(x[j]) == 0 {
			return 0
		}
		// NOTE Try len(x)-1 instead. Sounds worse on my $4 speakers.
		// FIXME Works ONLY with nbuf=4096, nfft=8192 (oversampling 2).
		return -real(xt[j]/x[j])/float64(len(x))*math.Pi*stretch - math.Pi/2
	}
}

func tdx(x, xd []complex128, j int, olap float64) float64 {
	if mag(x[j]) < 1e-6 {
		return 0
	}
	e := imag(xd[j]/x[j]) / (olap / 2)
	if abs(e) > olap/2 {
		e = math.Copysign(olap/2, e)
	}
	return e
}

func gettadv(x, xd []complex128, olap float64) func(w int) float64 {
	return func(j int) float64 {
		if mag(x[j]) < 1e-6 {
			return 0
		}
		return (math.Pi*float64(j) + imag(xd[j]/x[j])) / (olap / 2)
	}
}
