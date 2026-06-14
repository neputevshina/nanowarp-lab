module github.com/neputevshina/nanowarp-lab/pghi-xsynth

go 1.26.2

replace github.com/neputevshina/nanowarp/wav => ../../nanowarp/wav
replace github.com/neputevshina/nanowarp => ../../nanowarp

require (
	github.com/neputevshina/nanowarp/wav v0.0.0-20260612171917-dad73e50ae4f
	golang.org/x/exp v0.0.0-20260410095643-746e56fc9e2f
	gonum.org/v1/gonum v0.17.0
)

require (
	github.com/neputevshina/nanowarp v0.0.0-20260612171917-dad73e50ae4f // indirect
	github.com/youpy/go-riff v0.1.0 // indirect
	github.com/zaf/g711 v1.4.0 // indirect
	golang.org/x/sys v0.43.0 // indirect
)
