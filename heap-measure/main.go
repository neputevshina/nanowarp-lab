package main

import (
	"flag"
	"io"
	"os"

	"github.com/neputevshina/nanowarp/wav"
)

var finput = flag.String("i", "", "input WAV (or anything else, if ffmpeg is present) `path`")

func main() {
	flag.Parse()

	file, err := os.Open(*finput)
	if err != nil {
		panic(err)
	}

	wavrd := wav.NewReader(file)

	mid := []float64{}
	side := []float64{}
	for {
		samples, err := wavrd.ReadSamples()
		if err == io.EOF {
			break
		}

		for _, sample := range samples {
			l, r := wavrd.FloatValue(sample, 0), wavrd.FloatValue(sample, 1)
			mid = append(mid, l)
			side = append(side, r)
		}
	}

	sources := make([]float64, int(float64(len(mid))))
	diverges := make([]float64, int(float64(len(mid))))
	urs := make([]float64, int(float64(len(mid))))
	drs := make([]float64, int(float64(len(mid))))
	sinks := make([]float64, int(float64(len(mid))))

	f, err := wavrd.Format()
	if err != nil {
		panic(err)
	}
	fs := int(f.SampleRate)
	d := detectorNew(2048, fs)
	d.process2(mid, sources, diverges, urs, drs, sinks)

	for i := range sinks {
		sinks[i] = urs[i] + drs[i]
	}

	dump("sources.wav", sources, fs)
	dump("twos.wav", sinks, fs)
	dump("urs.wav", urs, fs)
	dump("drs.wav", drs, fs)
}
