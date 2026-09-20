package main

import (
	"io"
	"math"
	"math/cmplx"
	"slices"

	"github.com/neputevshina/geom"
	"github.com/neputevshina/nanowarp/dspio"
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

	points []geom.Point
	j      float64

	a wbufs
}

type wbufs struct {
	Mid, S, T                  []float64      // Scratch buffers
	Ph, M, P, F                []float64      `size:"nbins"` // Current phase
	Past, Future               []float64      // Phase accumulators
	Px, Py, Cone               []float64      `size:"nbins"`
	W, Wr, Wd, Wt, Wdt, Wcone  []float64      // Window functions
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

	j := int(w.j)
	{
		cnv := make2(j, w.nbins)
		for _, p := range w.points {
			e := &cnv[int(math.Round(clamp(0, w.j-1, p.Y)))][int(math.Round(clamp(0, float64(w.nbins-1), p.X)))]
			*e = min(*e+1, 8)
		}
		for _, sl := range cnv {
			slices.Reverse(sl)
			Oscope(sl, Name(`densbin`))
		}
	}
	if err == io.EOF {
		err = nil
	}
	return nil
}

func (n *warper) advance(presenta [][]float64) {
	a := &n.a

	n.j++
	n.analyze(presenta, a.C, a.Px, a.Py, a.M, a.Mid, a.Cone)

	for i := range n.nbins {
		n.points = append(n.points, geom.Pt(float64(i)+a.Px[i], n.j))
		// n.points = append(n.points, geom.Pt(a.Fadv[i], n.j))
		n.points = append(n.points, geom.Pt(float64(i), n.j+a.Py[i]))
	}

	// crop := func(e float64) float64 { return math.Copysign(bitsafe(math.Sqrt(abs(e))), e) }
	ue := slices.Clone(a.Px)
	uo := slices.Clone(a.Px)
	for i := range a.Px {
		a.Px[i] = clamp(-4, 4, a.Px[i])
		a.Py[i] = clamp(-4, 4, a.Py[i])
		ue[i], uo[i] = cmplx.Polar(complex(a.Px[i], a.Py[i]))
	}
	slices.Reverse((a.Px))
	slices.Reverse((a.Py))
	slices.Reverse((ue))
	slices.Reverse((uo))
	slices.Reverse(a.Cone)
	Oscope(slices.Clone(a.Px), Name(`x`))
	Oscope(slices.Clone(a.Py), Name(`y`))
	Oscope(slices.Clone(ue), Name(`mag`))
	Oscope(slices.Clone(uo), Name(`phase`))
	Oscope(slices.Clone(a.Cone), Name(`cone`))
	for i := range a.Px {
		a.M[i] = math.Log(a.M[i] + 1)
	}
	slices.Reverse(a.M)
	Oscope(slices.Clone(a.M), Name(`justmag`))
	// oscope.Oscope(slices.Clone(a.Tadv), oscope.Name(`x`))

	// n.integrate([][]float64{nil, a.Fadv}, [][]float64{nil, a.Tadv}, [][]float64{a.P, a.M}, [][]float64{a.Past, a.Ph}, n.arm)

	// for w := range a.Ph {
	// 	a.Past[w] = princarg(a.Ph[w])
	// }
	// copy(a.P, a.M)

}

func (n *warper) analyze(present [][]float64, C [][]complex128, Px, Py, M, Mid, Cone []float64) {
	a := &n.a

	clear(Mid)
	for ch := range present {
		floats.Add(Mid[:len(present[ch])], present[ch])
		n.enfft(C[ch], a.W, present[ch])
	}

	n.coneShapeKernel(Cone, Mid)

	n.enfft(a.X, a.W, Mid)
	n.enfft(a.Xd, a.Wd, Mid)
	n.enfft(a.Xt, a.Wt, Mid)

	for w := range a.X {
		// Default TFR:
		// Px[w] = cif(a.X, a.Xd, w) * n.osamp
		// Py[w] = -lgd(a.X, a.Xt, w) / float64(n.hop)

		// Amplidude reassignment, see section 3 of
		// “Hainsworth, Stephen, and Malcolm Macleod. Time frequency reassignment:
		// A review and analysis. University of Cambridge, Department of Engineering, 2003.”
		Px[w] = magvert(a.X, a.Xt, w) / float64(n.hop)
		Py[w] = maghor(a.X, a.Xd, w) * n.osamp
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

func (n *warper) coneShapeKernel(out, Mid []float64) {
	// https://en.wikipedia.org/wiki/Cone-shape_distribution_function
	a := &n.a

	n.enfft(a.X, nil, Mid)

	for w := range a.X {
		a.X[w] *= cmplx.Conj(a.X[w])
	}

	n.defft(a.T, a.X, nil, false)
	// slices.Reverse(a.T)
	n.enfft(a.X, a.Wcone, a.T)

	for w := range a.X {
		out[w] = math.Log1p(max(0, 4*real(a.X[w])))
	}
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
	copy(a.Wcone, s(a.W))
	a.Wcone[nbuf/2-1] *= 0.5

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

func (n *warper) defft(out []float64, x []complex128, w []float64, noscale bool) {
	a := &n.a
	n.fft.Sequence(a.S, x)
	if !noscale {
		floats.Scale(1./n.norm, a.S)
	}
	if w != nil {
		mul(a.S, w)
	}
	copy(out, a.S)
}

func lgd(x, xt []complex128, j int) float64 {
	if mag(x[j]) == 0 {
		return 0
	}
	e := -real(xt[j] / x[j])
	return e
}

func cif(x, xd []complex128, j int) float64 {
	if mag(x[j]) < 1e-6 {
		return 0
	}
	e := imag(xd[j]/x[j]) / math.Pi / 2
	return e
}

func magvert(x, xt []complex128, j int) float64 {
	if mag(x[j]) == 0 {
		return 0
	}
	e := imag(xt[j] / x[j])
	return e
}

func maghor(x, xd []complex128, j int) float64 {
	if mag(x[j]) < 1e-6 {
		return 0
	}
	e := real(xd[j]/x[j]) / math.Pi / 2
	return e
}
