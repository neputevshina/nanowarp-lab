package main

import (
	"flag"
	"os"

	"github.com/neputevshina/nanowarp/dspio/wavio"
)

var finputa = flag.String("a", "", "source a WAV")

func main() {
	Enable = true
	flag.Parse()

	fa, err := os.Open(*finputa)
	if err != nil {
		panic(err)
	}

	wsa, err := wavio.NewDecoder(fa)

	dt := warperNew(4096, 1, 8, 2)
	err = dt.process(wsa)
	if err != nil {
		panic(err)
	}

	err = Dump(nil, `disp`)
	if err != nil {
		panic(err)
	}
}
