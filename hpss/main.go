package main

import (
	"flag"
	"math"
	"math/cmplx"
	"os"

	"github.com/neputevshina/nanowarp-lab/common"
	"github.com/neputevshina/nanowarp-lab/common/window"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/dspio/wavio"
	"github.com/neputevshina/nanowarp/oscope"
	"github.com/neputevshina/nanowarp/pffft"
	"golang.org/x/exp/constraints"
)

var finputa = flag.String("a", "", "source WAV")
var fhrad = flag.Int("h", 0, "horizontal radius, if less than or 0 equal to -l")
var fvrad = flag.Int("v", 13, "vertical radius")
var folap = flag.Int("l", 5, "overlap")
var fnfft = flag.Int("nfft", 3600, "fft grain size")

var fvib = flag.Int("vib", 6, "vertical vibrato radius, 0 to disable")
var ftonal = flag.Bool("t", false, "output tonals")
var fhperc = flag.Int("ph", 50, "horizontal percentile")
var fvperc = flag.Int("pv", 50, "vertical percentile")

var fcif = flag.Bool("cif", true, "backtrack hpss mask by cif displacement")
var ffat = flag.Int("fat", 0, "dilate vertical activations horizontally, 0 disable, -1 → olap")

var fwin = flag.String("w", `kaiser`, "use a `window`: hann, welch, vorbis, tri, kaiser")

func main() {
	flag.Parse()
	oscope.Enable = true

	if *ffat == -1 {
		*ffat = *folap
	}

	name := "./perc.wav"
	if *ftonal {
		name = "./harm.wav"
	}
	fo, err := os.Create(name)
	if err != nil {
		panic(err)
	}

	nfft := pffft.NearestNfft(*fnfft)
	hop := nfft / *folap
	nbins := nfft/2 + 1

	wavr := common.INeedAWav(*finputa)
	p := wavr.Properties()
	wavw, err := wavio.NewEncoder(fo, p.Samplerate, p.Nch, wavio.FormatFloat, 32)
	if err != nil {
		panic(err)
	}
	defer wavw.Close()

	gr := dspio.NewOfflineGrainReader(nfft, hop, wavr)
	fft := pffft.New(nfft)

	// HPSS works best with windows with the smallest main lobe width,
	// but a window's derivative must have compact support (zeros only before and after the window)
	// to be perfectly invertible after modifications.
	// Meaning, a rectangular or a Tukey window will produce artifacts.
	win := make([]float64, nfft)
	switch *fwin {
	case `hann`:
		window.Hann(win)
	case `welch`:
		window.Welch(win)
	case `vorbis`:
		window.Vorbis(win)
	case `tri`:
		window.Tri(win)
	case `kaiser`:
		win = window.Kaiser(nfft, 3)
	default:
		panic(`unknown window type`)
	}

	if *fhrad == 0 {
		*fhrad = *folap
	}

	bank := newhoriz(nbins, *fhrad, float64(*fhperc)/100)
	vert := newvert(*fvrad, float64(*fvperc)/100)
	smolbank := newhoriz(nbins, *folap, 1)
	smolv := newvert(*fvib, 1)

	mid := make([]float64, nbins)
	a := make([]float64, nbins)
	b := make([]float64, nbins)
	c := make([]float64, nbins)
	d := make([]float64, nbins)

	first := true
	now := max(*fhrad, *ffat)

	mono := func(out []float64, in [][]complex128) {
		for w := range nbins {
			out[w] = .0
			for ch := range in {
				out[w] += cmplx.Abs(in[ch][w])
			}
		}
	}

	roto := make2[bool](now*2+1, nbins)
	newroto := make([]int, nbins)
	wins := [][]float64{win, window.Dx(win), window.T(win)}
	common.StftHandle(gr, wavw, fft, wins, now, now, func(out [][]complex128, windowed [][][][]complex128, preanalyze bool) {
		frames := windowed[0]

		// Accumulate lookahead.
		if first {
			for f := range frames[now:] {
				mono(mid, frames[f])
				if *fvib != 0 {
					smolv.Filt(mid)
				}
				bank.Insert(mid)
			}
			for f := range frames[now:] {
				mono(mid, frames[f])
				smolbank.Insert(mid)
			}
			first = false
		}
		if preanalyze {
			return
		}

		mono(mid, frames[now])
		if *ffat > 0 {
			smolbank.Insert(mid)
			smolbank.Take(d)
		} else {
			copy(d, mid)
		}
		vert.Filtout(d, a)

		bank.Take(b)
		if *fvib != 0 {
			smolv.Filtout(mid, c)
		}
		bank.Insert(c)

		for w := range nbins {
			roto[now][w] = a[w] > b[w]
		}

		if *fcif {
			// Backpropagate the mask vertically by vertical bin shift (CIF).
			clear(newroto)
			for w := range nbins {
				cif := imag(windowed[1][now][0][w]/frames[now][0][w]) / math.Pi / 2
				y := max(0, min(int(math.Round(float64(w)+cif)), nbins-1))
				if roto[now][y] {
					newroto[w]++
				} else {
					newroto[w]--
				}
				if roto[now][w] {
					newroto[y]++
				} else {
					newroto[y]--
				}
			}
		}
		for w := range nbins {
			m := boolfloat(newroto[w] > 0)
			if *ftonal {
				m = 1 - m
			}
			for ch := range out {
				out[ch][w] = frames[now][ch][w] * complex(m, 0)
			}
		}

		roto[len(roto)-1], roto[0] = roto[0], roto[len(roto)-1]
		copy(roto, roto[1:])

	})

	oscope.Dump(nil, ".")
}

type horiz struct {
	a []*mediator[float64, struct{}]
}

func newhoriz(nbins, n int, pc float64) *horiz {
	bank := make([]*mediator[float64, struct{}], nbins)
	for i := range bank {
		bank[i] = mediatorNew[float64, struct{}](n, n, pc)
	}
	return &horiz{a: bank}
}

func (b *horiz) Insert(f []float64) {
	for w := range len(f) {
		b.a[w].Insert(f[w], struct{}{})
	}
}

func (b *horiz) Take(f []float64) {
	for w := range len(f) {
		f[w], _ = b.a[w].Take()
	}
}

func (b *horiz) Filt(f []float64) {
	b.Insert(f)
	b.Take(f)
}

func (b *horiz) Filtout(in, out []float64) {
	b.Insert(in)
	b.Take(out)
}

type vert struct {
	a *mediator[float64, struct{}]
}

func newvert(n int, pc float64) *vert {
	return &vert{a: mediatorNew[float64, struct{}](n, n, pc)}
}

func (b *vert) Filt(f []float64) {
	b.a.Reset(b.a.N)
	for w := b.a.N - 1; w >= 0; w-- { // Mirror the beginning
		b.a.Insert(f[w], struct{}{})
	}
	for w := range f {
		f[w], _ = b.a.Filt(f[w], struct{}{})
	}
}

func (b *vert) Filtout(in, out []float64) {
	b.a.Reset(b.a.N)
	for w := b.a.N - 1; w >= 0; w-- { // Mirror the beginning
		b.a.Insert(in[w], struct{}{})
	}
	for w := range in {
		out[w], _ = b.a.Filt(in[w], struct{}{})
	}
}

const (
	right = 1 << iota
	down
	up
	topleft
	topright
	ridgemask = right | down | up | topleft | topright
)

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

func grainwrap(nfft int, x float64) float64 {
	nf := float64(nfft)
	r := math.Mod(x+nf/2, nf)
	if r < 0 {
		r += nf
	}
	r -= nf / 2
	if r != r || math.IsInf(r, 0) {
		return 0
	}
	return r
}
func clamp[T constraints.Ordered](a, b, x T) T {
	return max(a, min(b, x))
}
