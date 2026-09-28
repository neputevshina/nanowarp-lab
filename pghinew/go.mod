module github.com/neputevshina/nanowarp-lab/pghinew

go 1.26.2

require github.com/neputevshina/nanowarp v0.6.13

require gonum.org/v1/gonum v0.17.0 // indirect

replace github.com/neputevshina/nanowarp-lab/common => ../common

require (
	github.com/neputevshina/nanowarp-lab/common v0.0.0-20260928105720-cdd690ff1d27
	golang.org/x/exp v0.0.0-20251002181428-27f1f14c8bb9 // indirect
	golang.org/x/sys v0.47.0 // indirect
)
