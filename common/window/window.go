package window

import (
	"math"
	"slices"

	"golang.org/x/exp/constraints"
	"gonum.org/v1/gonum/dsp/fourier"
)

// Niemitalo is an asymmetric windowing function, designed by Olli Niemitalo.
//
// See https://dsp.stackexchange.com/questions/2337/fft-with-asymmetric-windowing
func Niemitalo(out []float64) {
	nfft := float64(len(out))
	clear(out)
	sin, cos := math.Sin, math.Cos
	for i := nfft / 4; i < nfft*7/8; i++ {
		x := 2 * math.Pi * ((i+0.5)/nfft - 1.75)
		out[int(i)] = 2.57392230162633461887 - 1.58661480271141974718*cos(x) + 3.80257516644523141380*sin(x) -
			1.93437090055110760822*cos(2*x) - 3.27163999159752183488*sin(2*x) + 3.26617449847621266201*cos(3*x) -
			0.30335261753524439543*sin(3*x) - 0.92126091064427817479*cos(4*x) + 2.33100177294084742741*sin(4*x) -
			1.19953922321306438725*cos(5*x) - 1.25098147932225423062*sin(5*x) + 0.99132076607048635886*cos(6*x) -
			0.34506787787355830410*sin(6*x) - 0.04028033685700077582*cos(7*x) + 0.55461815542612269425*sin(7*x) -
			0.21882110175036428856*cos(8*x) - 0.10756484378756643594*sin(8*x) + 0.06025986430527170007*cos(9*x) -
			0.05777077835678736534*sin(9*x) + 0.00920984524892982936*cos(10*x) + 0.01501989089735343216*sin(10*x)
	}
	for i := 0; i < int(nfft)/8; i++ {
		nfft := int(nfft)
		out[nfft-1-i] = (1 - out[nfft*3/4-1-i]*out[nfft*3/4+i]) / out[nfft/2+i]
	}
	copy(out, out[int(nfft)*2/8:])
	clear(out[int(nfft)*6/8:])
}

// Hann is a Hann window function.
//
// Dual to itself (achieves COLA/perfect reconstruction) at any integer overlap greater than 2.
func Hann(out []float64) {
	for i := range out {
		x := float64(i) / float64(len(out))
		out[i] = 0.5 * (1 - math.Cos(2*math.Pi*x))
	}
}

// Tri is a triangular window function.
func Tri(win []float64) {
	nfft := len(win)
	for w := range nfft / 2 {
		win[w] = unmix(0, float64(nfft/2), float64(w))
	}
	for w := range nfft / 2 {
		win[nfft/2:][w] = unmix(float64(nfft/2), 0, float64(w))
	}
}

// Vorbis is a window function used in Vorbis codec.
func Vorbis(win []float64) {
	nfft := len(win)
	for w := range nfft {
		x := unmix(0, float64(nfft), float64(w))
		c := math.Cos(math.Pi * (x + 0.5))
		win[w] = math.Sin(math.Pi / 2 * c * c)
	}
}

// Welch is a Welch (parabola) window function.
func Welch(win []float64) {
	nfft := len(win)
	for w := range nfft {
		x := unmix(0, float64(nfft), float64(w))
		win[w] = 1 - (2*x-1)*(2*x-1)
	}
}

// BlackmanHarrisClassic is, surprisingly, a three-term Nutall minimum sidelobe window function
// also known as “three-term Blackman-Harris window”[1].
//
// Dual to itself (achieves COLA/perfect reconstruction) at overlap 6.
//
// Fun fact: frequency response in the article “Window functon” in English Wikipedia is for
// four-term -92 dB Blackman-Harris. And the picture of it was created and uploaded by Olli Niemitalo,
// the creator of asymmetric window function defined below. Correct parameters for the pictured window
// can be found in German version of this article. This window still sounds better than its four-term version.
//
// [1]: Doerry, Armin W. Catalog of window taper functions for sidelobe control. No. SAND2017-4042.
// Sandia National Laboratories (SNL-NM), Albuquerque, NM (United States), 2017.
// https://www.osti.gov/servlets/purl/1365510
func BlackmanHarrisClassic(out []float64) {
	// See [fadv].
	for i := range out {
		x := float64(i) / float64(len(out)-1)
		out[i] = .4243801 - .4973406*math.Cos(2*math.Pi*x) + .0782793*math.Cos(4*math.Pi*x)
	}
}

// BlackmanHarris92dB is a four-term -92 dB Blackman-Harris window function.
// See blackmanHarrisClassic.
func BlackmanHarris92dB(out []float64) {
	for i := range out {
		x := float64(i) / float64(len(out))
		out[i] = .35875 - .48829*math.Cos(2*math.Pi*x) + .14128*math.Cos(4*math.Pi*x) - .01168*math.Cos(6*math.Pi*x)
	}
}

// AvciNacaroglu is an Avci-Nacaroglu window function [1].
// It is related to DPSS window like Kaiser, but doesn't require
// modified Bessel function to calculate.
// Duality to itself is achieved by tuning the parameter a.
// For example, COLA at overlap 4 can approximately be achieved with a = 1.78.
//
// [1]: Doerry, Armin W. Catalog of window taper functions for sidelobe control. No. SAND2017-4042.
// Sandia National Laboratories (SNL-NM), Albuquerque, NM (United States), 2017.
// https://www.osti.gov/servlets/purl/1365510
func AvciNacaroglu(out []float64, a float64) {
	for i := range out {
		x := float64(i) / float64(len(out))
		out[i] = math.Exp(math.Pi*a*math.Sqrt(1-(4*(x-0.5)*(x-0.5))) - 1)
	}
	m := slices.Max(out)
	for i := range out {
		out[i] /= m
	}
}

// T calculates the time ramp multiplied window
// function for time-frequency reassignment.
func T(w []float64) (out []float64) {
	out = make([]float64, len(w))
	n := float64(len(w))
	for i := range w {
		out[i] = w[i] * mix(-n/2, n/2+1, float64(i)/n)
	}
	return
}

// Dx calculates the derivative of w using DFT.
func Dx(w []float64) (out []float64) {
	out = make([]float64, len(w))
	f := fourier.NewFFT(len(w))
	s := f.Coefficients(nil, w)
	for i := range s {
		s[i] *= complex(0, float64(i))
		s[i] /= complex(float64(len(w)), 0)
	}
	s[len(w)/2] = 0
	f.Sequence(out, s)
	for i := range out {
		out[i] = -out[i] * math.Pi * 2
	}
	return
}

// Kaiser is a Kaiser window function.
//
// See https://en.wikipedia.org/wiki/Kaiser_window
func Kaiser(nfft int, α float64) (w []float64) {
	w = make([]float64, nfft)
	α *= math.Pi
	q := i0(α)
	for i := range w {
		x := float64(i) / float64(len(w))
		w[i] = i0(α*math.Sqrt(1-(2*x-1)*(2*x-1))) / q
	}
	return w
}

// Modified Bessel function of order zero.
//
// Required to calculate Kaiser window.
//
// Based on i0.c and chbevl.c from https://netlib.org/cephes/.
// Original copyright:
//
//	Cephes Math Library Release 2.8:  June, 2000
//	Copyright 1984, 1987, 2000 by Stephen L. Moshier

// i0 is a zeroth-order modified Bessel function of the first kind.
func i0(x float64) float64 {
	x = math.Abs(x)
	if x <= 8 {
		y := x/2 - 2
		return math.Exp(x) * chbevl(y, i0A)
	}
	return math.Exp(x) * chbevl(32./x-2, i0B) / math.Sqrt(x)
}

// chbevl evaluates Chebyshev series.
func chbevl(x float64, arr []float64) float64 {
	b0, b1, b2 := arr[0], 0., 0.
	for i := 1; i < len(arr); i++ {
		b2 = b1
		b1 = b0
		b0 = x*b1 - b2 + arr[i]
	}
	return 0.5 * (b0 - b2)
}

// Chebyshev coefficients for exp(-x) I0(x)
// in the interval [0,8].
//
//	lim(x->0){ exp(-x) I0(x) } = 1.
var i0A = []float64{
	-4.41534164647933937950e-18,
	3.33079451882223809783e-17,
	-2.43127984654795469359e-16,
	1.71539128555513303061e-15,
	-1.16853328779934516808e-14,
	7.67618549860493561688e-14,
	-4.85644678311192946090e-13,
	2.95505266312963983461e-12,
	-1.72682629144155570723e-11,
	9.67580903537323691224e-11,
	-5.18979560163526290666e-10,
	2.65982372468238665035e-9,
	-1.30002500998624804212e-8,
	6.04699502254191894932e-8,
	-2.67079385394061173391e-7,
	1.11738753912010371815e-6,
	-4.41673835845875056359e-6,
	1.64484480707288970893e-5,
	-5.75419501008210370398e-5,
	1.88502885095841655729e-4,
	-5.76375574538582365885e-4,
	1.63947561694133579842e-3,
	-4.32430999505057594430e-3,
	1.05464603945949983183e-2,
	-2.37374148058994688156e-2,
	4.93052842396707084878e-2,
	-9.49010970480476444210e-2,
	1.71620901522208775349e-1,
	-3.04682672343198398683e-1,
	6.76795274409476084995e-1,
}

// Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
// in the inverted interval [8,infinity].
//
//	lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
var i0B = []float64{
	-7.23318048787475395456e-18,
	-4.83050448594418207126e-18,
	4.46562142029675999901e-17,
	3.46122286769746109310e-17,
	-2.82762398051658348494e-16,
	-3.42548561967721913462e-16,
	1.77256013305652638360e-15,
	3.81168066935262242075e-15,
	-9.55484669882830764870e-15,
	-4.15056934728722208663e-14,
	1.54008621752140982691e-14,
	3.85277838274214270114e-13,
	7.18012445138366623367e-13,
	-1.79417853150680611778e-12,
	-1.32158118404477131188e-11,
	-3.14991652796324136454e-11,
	1.18891471078464383424e-11,
	4.94060238822496958910e-10,
	3.39623202570838634515e-9,
	2.26666899049817806459e-8,
	2.04891858946906374183e-7,
	2.89137052083475648297e-6,
	6.88975834691682398426e-5,
	3.36911647825569408990e-3,
	8.04490411014108831608e-1,
}

// mix is a linear interpolation.
func mix[F constraints.Float](a, b, x F) F {
	return a*(1-x) + b*x
}

// unmix is a linear extrapolation.
func unmix[F constraints.Float](a, b, x F) F {
	if a == b {
		return a
	}
	return (x - a) / (b - a)
}
