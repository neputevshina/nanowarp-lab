package main

import (
	"fmt"
	"math"
	"os"
	"slices"

	"gonum.org/v1/gonum/dsp/fourier"
)

type detector struct {
	nfft        int
	nbuf        int
	nbins       int
	hop         int
	norm, wgain float64
	corr        float64
	thresh      float64
	fs          int

	fft  *fourier.FFT
	arm  []bool
	heap hp
	img  [][]float64

	a dbufs
}
type dbufs struct {
	S         []float64
	W, Wf, Wr []float64
	X, Y      []complex128
}

func detectorNew(nfft, fs int) (n *detector) {
	corr := math.Ceil(float64(fs) / 48000)
	nbuf := nfft * int(corr)
	nbins := nfft/2 + 1
	olap := 16

	n = &detector{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		corr:  corr,
		fs:    fs,
	}
	makeslices(&n.a, nbins, nfft)

	// Asymmetric window requires applying reversed copy of itself on synthesis stage.
	niemitalo(n.a.Wf)
	copy(n.a.Wr, n.a.Wf)
	slices.Reverse(n.a.Wr)

	n.wgain = windowGain(n.a.Wf)
	n.norm = float64(nfft) * float64(olap) * n.wgain
	n.fft = fourier.NewFFT(nfft)
	n.a.W = n.a.Wf
	n.arm = make([]bool, nbins)

	return
}

func (n *detector) process2(lin, a, b, c, d, e []float64) {
	fmt.Fprintln(os.Stderr, `(*detector).process`)

	t := make([]float64, n.nfft)
	for i := 0; i < len(lin); i += n.hop {
		pgrain := lin[min(len(lin)-n.hop, i):min(len(lin)-n.hop, i+n.nbuf-n.hop)]
		grain := lin[i:min(len(lin), i+n.nbuf)]
		so, div, ur, dr, sk := n.advance(pgrain, grain)

		add := func(v float64, o []float64) {
			fill(t, v)
			mul(t, n.a.Wr)
			add(o[i:min(len(lin), i+n.nbuf)], t)
		}

		add(so, a)
		add(div, b)
		add(ur, c)
		add(dr, d)
		add(sk, e)
	}

	return
}

func (n *detector) advance(pingrain, ingrain []float64) (nsource, ndiverge, nupright, ndownright, nsink float64) {
	a := &n.a
	enfft := func(x []complex128, w, grain []float64) {
		clear(a.S)
		copy(a.S, grain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	enfft(a.X, a.W, ingrain)
	enfft(a.Y, a.W, pingrain)

	n.heap = make(hp, n.nbins)
	clear(n.arm)

	for j := range a.X {
		n.arm[j] = true
		n.heap[j] = heaptriple{mag(a.Y[j]), j, -1}
	}
	heapInit(&n.heap)

	field := make([]uint, n.nbins)
	const (
		right = 1 << iota
		up
		down
	)

	// This is just PGHI without phase accumulation.
	for len(n.heap) > 0 {
		h := heapPop(&n.heap).(heaptriple)
		w := h.w
		switch h.t {
		case -1:
			if n.arm[w] {
				field[w] |= right
				n.arm[w] = false
				heapPush(&n.heap, heaptriple{mag(a.X[w]), w, 0})
			}
		case 0:
			if w > 1 && n.arm[w-1] {
				field[w] |= down
				n.arm[w-1] = false
				heapPush(&n.heap, heaptriple{mag(a.X[w-1]), w - 1, 0})
			}
			if w < n.nbins-1 && n.arm[w+1] {
				field[w] |= up
				n.arm[w+1] = false
				heapPush(&n.heap, heaptriple{mag(a.X[w+1]), w + 1, 0})
			}
		}
	}
	for i, v := range field {
		unit := min(1, 1/float64(i))
		switch v {
		case 0:
			nsink += unit
		case up | right:
			nupright += unit
		case down | right:
			ndownright += unit
		case up | down:
			ndiverge += unit
		case up | down | right:
			nsource += unit
		}
	}
	return
}
