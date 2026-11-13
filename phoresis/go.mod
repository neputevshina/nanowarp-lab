module github.com/neputevshina/nanowarp-lab/phoresis

go 1.27.2

replace github.com/neputevshina/nanowarp-lab/common => ../common

require (
	github.com/neputevshina/nanowarp v0.6.17
	github.com/neputevshina/nanowarp-lab/common v0.0.0-20260929121432-2bd6e659f03c
	golang.org/x/exp v0.0.0-20251002181428-27f1f14c8bb9
	gonum.org/v1/gonum v0.17.0
)

require (
	golang.org/x/sys v0.47.0 // indirect
	golang.org/x/term v0.45.0 // indirect
)
