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

func (n *detector) process2(lin, a, b, c []float64) {
	fmt.Fprintln(os.Stderr, `(*detector).process`)

	t := make([]float64, n.nfft)
	for i := 0; i < len(lin); i += n.hop {
		ax, bx, cx := n.advance(lin[min(len(lin)-n.hop, i):min(len(lin)-n.hop, i+n.nbuf-n.hop)], lin[i:min(len(lin), i+n.nbuf)])

		fill(t, ax)
		mul(t, n.a.Wr)
		add(a[i:min(len(lin), i+n.nbuf)], t)

		fill(t, bx)
		mul(t, n.a.Wr)
		add(b[i:min(len(lin), i+n.nbuf)], t)

		fill(t, cx)
		mul(t, n.a.Wr)
		add(c[i:min(len(lin), i+n.nbuf)], t)
	}

	return
}

func (n *detector) advance(pingrain, ingrain []float64) (up, down, right float64) {
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

	// This is just PGHI without phase accumulation.
	for len(n.heap) > 0 {
		h := heapPop(&n.heap).(heaptriple)
		w := h.w
		switch h.t {
		case -1:
			if n.arm[w] {
				right += mag(a.Y[w])
				n.arm[w] = false
				heapPush(&n.heap, heaptriple{mag(a.X[w]), w, 0})
			}
		case 0:
			if w > 1 && n.arm[w-1] {
				up += mag(a.X[w-1])
				n.arm[w-1] = false
				heapPush(&n.heap, heaptriple{mag(a.X[w-1]), w - 1, 0})
			}
			if w < n.nbins-1 && n.arm[w+1] {
				down += mag(a.X[w+1])
				n.arm[w+1] = false
				heapPush(&n.heap, heaptriple{mag(a.X[w+1]), w + 1, 0})
			}
		}
	}

	return
}
