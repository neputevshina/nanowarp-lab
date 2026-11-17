package common

import (
	"fmt"
	"os"
	"os/exec"
	"path"
	"path/filepath"

	"github.com/neputevshina/nanowarp/dspio/wavio"
)

func INeedAWav(fp string) *wavio.Decoder {
	fa, err := os.Open(fp)
	if err != nil {
		panic(err)
	}

	wavr, err := wavio.NewDecoder(fa)
	if err != nil {
		// Try to call ffmpeg, we've probably got an MP3.
		if err == wavio.ErrNotAWav {
			_ = fa.Close()

			_, err = exec.LookPath(`ffmpeg`)
			if err != nil {
				fmt.Fprintln(os.Stderr, `can process only WAV files without ffmpeg`)
				os.Exit(1)
			}

			s2, err := filepath.Abs(fp)
			if err != nil {
				panic(err)
			}

			ex := path.Join(os.TempDir(), path.Base(s2)) + `.wav`

			cmd := exec.Command(`ffmpeg`, `-v`, `warning`, `-y`, `-i`, s2, `-acodec`, `pcm_f32le`, ex)
			cmd.Stderr = os.Stderr
			cmd.Stdout = os.Stdout
			err = cmd.Run()
			if err != nil {
				if _, ok := err.(*exec.ExitError); ok {
					os.Exit(1)
				} else {
					panic(err)
				}
			}

			fa, err = os.Open(ex)
			if err != nil {
				panic(err)
			}
			wavr, err = wavio.NewDecoder(fa)
			if err != nil {
				panic(err)
			}
		}
	}
	return wavr
}
