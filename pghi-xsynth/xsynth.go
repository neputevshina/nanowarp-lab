package main

import (
	"io"
	"math"
	"math/cmplx"
	"slices"

	"github.com/neputevshina/nanowarp/dspio"
	"gonum.org/v1/gonum/cmplxs"
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
	S, Mid, M, P, F            []float64 // Scratch buffers
	Ph                         []float64 `size:"nbins"` // Current phase
	Past, Future               []float64 // Phase accumulators
	Fadv, Tadv                 []float64
	W, Wr, Wd, Wt, Wdt         []float64      // Window functions
	X, Y, Xd, Xt, L, R, Lo, Ro []complex128   // Complex spectra
	C, Co                      [][]complex128 // Channels

	Cs                    [][][]complex128
	Phs, Fadvs, Tadvs, Ms [][]float64 `size:"lah"`
}

func (w *warper) process(ain, bin dspio.SignalReader, out dspio.SignalWriter, pm, mm *float64) error {
	agr := dspio.NewGrainReader(w.nbuf, w.hop, ain)
	bgr := dspio.NewGrainReader(w.nbuf, w.hop, bin)
	gw := dspio.NewGrainWriter(w.nbuf, w.hop, out)

	agrain := make2(ain.NchRead(), w.nbuf)
	bgrain := make2(ain.NchRead(), w.nbuf)
	outgrain := make2(ain.NchRead(), w.nbuf)
	var err error
	for {
		_, err = agr.SignalRead(nil, agrain)
		_, err = bgr.SignalRead(err, bgrain)
		if err != nil {
			break
		}
		w.advance(agrain, bgrain, outgrain, *pm, *mm)
		_, err = gw.SignalWrite(nil, outgrain)
		if err != nil {
			break
		}
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
		lah:   olap*120 + 1,
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

	// waveform.Dump(nil, a.W)
	// waveform.Dump(nil, a.Wt)
	// waveform.Dump(nil, a.Wd)

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

func (n *warper) bypassGrain(present, output [][]float64) {
	a := &n.a
	for ch := range present {
		// This is the only known way to correctly scale gains.
		n.enfft(a.Co[ch], a.W, present[ch])
		n.defft(output[ch], a.Co[ch], true)
		// And this works with -40 dB difference.
		// copy(output[ch], present[ch])
		// mul(output[ch], a.W)
		// mul(output[ch], a.Wr)
		// scale(output[ch], n.wgain*float64(n.nbuf)/float64(n.hop))
		// Probably there is something wrong in gonum FFT implementation.
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

func (n *warper) analyze(present [][]float64, C [][]complex128, Fadv, Tadv, M, Mid []float64, speedup float64) {
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

	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
	for w := range a.X {
		// TODO Probably it will be more numerically stable to limit the phase accuum to
		// 0..1 and scale back to -π..π range at the poltocar conversion.
		// fadv must return 0..1 accordingly, simply defer π multiplication till the end.
		Fadv[w] = princarg(getfadv(a.X, a.Xt, 2./n.osamp/speedup)(w))
		Tadv[w] = gettadv(a.X, a.Xd, float64(n.olap)*n.osamp)(w)
	}

	for w := range a.X {
		// Encode stereo phase differences and stretch mid only, keep original magnitudes.
		// NB: Phase difference in polar coordinates is complex division in cartesian.
		//     Phase sum is conversely a multiply.
		//     Hypot and multiplication are always cheaper than Atan2 and Sincos.
		//
		// See “Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a
		// real-time realization for a commercial music product.”
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

func (n *warper) crossAnalyze(phasemix, magnitudemix float64, presenta [][]float64, presentb [][]float64, C [][]complex128, Fadv, Tadv, M, Mid []float64) {
	a := &n.a

	clear(Mid)
	for ch := range presenta {
		floats.Add(Mid, presenta[ch])
		n.enfft(C[ch], a.W, presenta[ch])
	}

	n.enfft(a.X, a.W, Mid)
	n.enfft(a.Xd, a.Wd, Mid)
	n.enfft(a.Xt, a.Wt, Mid)

	// Stretch coefficient in cross-synthesis is 1, we're not changing the duration of sound.
	for w := range a.X {
		Fadv[w] = getfadv(a.X, a.Xt, 2./n.osamp)(w)
		Tadv[w] = gettadv(a.X, a.Xd, float64(n.olap)*n.osamp)(w)
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

	// Repeating with presentb:
	clear(Mid)
	for ch := range presentb {
		floats.Add(Mid, presentb[ch])
		n.enfft(C[ch], a.W, presentb[ch])
	}

	n.enfft(a.X, a.W, Mid)
	n.enfft(a.Xd, a.Wd, Mid)
	n.enfft(a.Xt, a.Wt, Mid)

	// Mix the phase.
	for w := range a.X {
		Fadv[w] = princarg(mix(Fadv[w], getfadv(a.X, a.Xt, 2./n.osamp)(w), phasemix))
		Tadv[w] = mix(Tadv[w], gettadv(a.X, a.Xd, float64(n.olap)*n.osamp)(w), phasemix)
	}

	for w := range a.X {
		m := mag(a.X[w])
		p := a.X[w] / complex(m, 0)
		if m < 1e-6 {
			p = complex(1, 0)
		}
		M[w] = mix(M[w], m, magnitudemix)
		a.Y[w] = p
	}
	for ch := range presenta {
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

func (n *warper) synthesize(output [][]float64, C [][]complex128, Ph []float64) {
	a := &n.a
	for w := range Ph {
		// Add stereo phase differences back through complex multiplication.
		a.Y[w] = cmplx.Rect(1, Ph[w])
	}
	for ch := range output {
		cmplxs.MulTo(a.X, C[ch], a.Y)
		n.defft(output[ch], a.X, true)
	}
}

func (n *warper) advance(presenta, presentb, output [][]float64, pm, mm float64) {
	a := &n.a

	n.crossAnalyze(pm, mm, presenta, presentb, a.C, a.Fadv, a.Tadv, a.M, a.Mid)

	n.integrate([][]float64{nil, a.Fadv}, [][]float64{nil, a.Tadv}, [][]float64{a.P, a.M}, [][]float64{a.Past, a.Ph}, n.arm)

	for w := range a.Ph {
		a.Past[w] = princarg(a.Ph[w])
	}
	copy(a.P, a.M)

	n.synthesize(output, a.C, a.Ph)
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

func gettadv(x, xd []complex128, olap float64) func(w int) float64 {
	return func(j int) float64 {
		if mag(x[j]) < 1e-6 {
			return 0
		}
		return (math.Pi*float64(j) + imag(xd[j]/x[j])) / (olap / 2)
	}
}
