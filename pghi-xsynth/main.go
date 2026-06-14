package main

import (
	"flag"
	"io"
	"os"

	"github.com/neputevshina/nanowarp/wav"
)

var finputa = flag.String("a", "", "source a WAV (or anything else, if ffmpeg is present) `path`")
var finputb = flag.String("b", "", "source b WAV (or anything else, if ffmpeg is present) `path`")
var pm = flag.String("p", "", "phase mix")
var mm = flag.String("m", "", "magnitude mix")

func main() {
	flag.Parse()

	file, err := os.Open(*finputa)
	if err != nil {
		panic(err)
	}

	wavrd := wav.NewReader(file)

	mid := []float64{}
	side := []float64{}
	for {
		_, samples, err := wavrd.ReadSamples(nil)
		if err == io.EOF {
			break
		}

		for _, sample := range samples {
			l, r := wavrd.FloatValue(sample, 0), wavrd.FloatValue(sample, 1)
			mid = append(mid, l)
			side = append(side, r)
		}
	}

	ups := make([]float64, int(float64(len(mid))))
	downs := make([]float64, int(float64(len(mid))))
	rights := make([]float64, int(float64(len(mid))))

	f, err := wavrd.Format()
	if err != nil {
		panic(err)
	}
	fs := int(f.SampleRate)
	// d := warperNew(4096, 2, 4, 2)
	// d.process1000(mid, ups, downs, rights)

	vert := make([]float64, int(float64(len(mid))))
	sub := make([]float64, int(float64(len(mid))))
	div := make([]float64, int(float64(len(mid))))
	for i := range ups {
		vert[i] = bitsafe(ups[i] + downs[i])
		sub[i] = bitsafe(ups[i] - downs[i])
		div[i] = bitsafe(rights[i] / vert[i])
	}

	dump("ups.wav", ups, fs)
	dump("downs.wav", downs, fs)
	dump("rights.wav", rights, fs)
	dump("vert.wav", vert, fs)
	dump("div.wav", div, fs)
	dump("sub.wav", sub, fs)

}
