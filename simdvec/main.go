package main

import (
	"fmt"
	"simd"
	"unsafe"
)

func add(a, b []simd.Float32s) {
	for i := range a {
		a[i] = a[i].Add(b[i])
	}
}

func main() {
	s1 := []float32{0, 1, 2, 3, 4, 5, 6, 7}
	s2 := []float32{7, 6, 5, 4, 3, 2, 1, 0}
	vs := simd.VectorBitSize() / 32
	vs1 := make([]simd.Float32s, vs)
	for i := range len(s1) / vs {
		vs1[i] = simd.LoadFloat32s(s1[i*vs:])
	}
	vs2 := make([]simd.Float32s, vs)
	for i := range len(s2) / vs {
		vs2[i] = simd.LoadFloat32s(s2[i*vs:])
	}
	add(vs1, vs2)
	for i := range len(s1) / vs {
		vs1[i].Store(s1[i*vs:])
	}
	fmt.Println(s1, vs, simd.Emulated())
	a := uintptr(unsafe.Pointer(unsafe.SliceData(vs1)))
	b := uintptr(unsafe.Pointer(unsafe.SliceData(vs2)))
	fmt.Println(a, b, a&0xf, b&0xf)
}
