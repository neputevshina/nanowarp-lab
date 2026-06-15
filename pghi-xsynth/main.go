package main

import (
	"flag"
	"os"
	"strconv"
)

var finputa = flag.String("a", "", "source a WAV (or anything else, if ffmpeg is present) `path`")
var finputb = flag.String("b", "", "source b WAV (or anything else, if ffmpeg is present) `path`")
var pm = flag.String("p", "", "phase mix")
var mm = flag.String("m", "", "magnitude mix")

func main() {
	flag.Parse()

	fa, err := os.Open(*finputa)
	if err != nil {
		panic(err)
	}
	fb, err := os.Open(*finputb)
	if err != nil {
		panic(err)
	}

	// Center-dilated novelty curve
	of, err := os.Create(`xsynth.wav`)
	defer of.Close()

	wsa, err := NewWavSignalReader(err, fa)
	wsb, err := NewWavSignalReader(err, fb)
	wsw, err := NewWavSignalWriter(err, of, func() int { return int(min(wsa.Size, wsb.Size)) },
		delay(2), func() int { return int(max(wsa.SampleRate, wsb.SampleRate)) })
	if err != nil {
		panic(err)
	}
	dt := warperNew(512, 2, 4, 2)
	pm, err := strconv.ParseFloat(*pm, 64)
	if err != nil {
		panic(err)
	}
	mm, err := strconv.ParseFloat(*mm, 64)
	if err != nil {
		panic(err)
	}
	err = dt.process(wsa, wsb, wsw, &pm, &mm)

}
