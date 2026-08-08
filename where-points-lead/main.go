package main

import (
	"flag"
	"os"

	"github.com/neputevshina/nanowarp/dspio/wavio"
	"github.com/neputevshina/nanowarp/oscope"
)

var finputa = flag.String("a", "", "source a WAV")

func main() {
	oscope.Enable = true
	flag.Parse()

	fa, err := os.Open(*finputa)
	if err != nil {
		panic(err)
	}

	wsa, err := wavio.NewDecoder(fa)

	dt := warperNew(4096, 1, 4, 2)
	err = dt.process(wsa)
	if err != nil {
		panic(err)
	}

	err = oscope.Dump(nil, `disp`)
	if err != nil {
		panic(err)
	}
}
