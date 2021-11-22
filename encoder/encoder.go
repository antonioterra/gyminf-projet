package encoder

import (
	"fmt"
	"math"
	"math/cmplx"

	"kazat.ch/lbcrypto/cMat"
	"kazat.ch/lbcrypto/polyMod"
)

type Encoder struct {
	M               float64
	mod             int
	scale           complex128
	ENC, DEC, basis cMat.CMat
}

func NewEncoer(M float64, scale complex128, mod int) Encoder {

	N := int(math.Floor(M / 2))
	exposant := complex(0, 2.0*math.Pi/M)
	xi := cmplx.Exp(exposant)
	v := make([]complex128, N*N)

	for i := 0; i < N; i++ {
		root := cmplx.Pow(xi, complex(float64(2*i)+1, 0))
		for j := 0; j < N; j++ {
			v[N*i+j] = cmplx.Pow(root, complex(float64(j), 0))
		}
	}

	DEC := cMat.NewCMat(N, N, v)
	R, S := DEC.Invert()
	ENC := cMat.Merge(&R, &S)

	var truc []complex128
	data := make([]complex128, N)
	for i := 0; i < N; i++ {
		if i > 0 {
			data[i-1] = 0
		}
		data[i] = 1
		G := cMat.NewCMat(N, 1, data)
		truc = append(truc, G.Mult(&DEC, &G).GetData()...)
	}

	basis := cMat.NewCMat(N, N, truc)

	res := Encoder{
		M:     M,
		mod:   mod,
		ENC:   ENC,
		DEC:   DEC,
		scale: scale,
		basis: basis,
	}
	return res
}

func (e *Encoder) Encode(v *cMat.CMat) polyMod.PolyInt {

	originaldata := v.GetData()
	l := len(originaldata)
	data := make([]complex128, 2*l)

	for i, val := range originaldata {

		data[i] = val
		data[i+l] = cmplx.Conj(originaldata[l-1-i])
	}

	newv := cMat.NewCMat(2*l, 1, data)
	fmt.Println("encoding :")
	newv.PPrint()
	newv.Scale(e.scale)

	newv = *cMat.ProjectOnRaws(&newv, &e.basis)
	newv.Mult(&e.ENC, &newv)

	return ToPol(&newv, e.mod)
}

func (e *Encoder) Decode(pol polyMod.PolyInt) cMat.CMat {

	v := ToMat(&pol)
	v.Mult(&e.DEC, &v)

	v.Scale(1 / e.scale)

	return v
}

func ToPol(m *cMat.CMat, mod int) polyMod.PolyInt {
	srcData := m.GetData()
	dstData := make([]int, len(srcData))
	for k, coef := range srcData {
		dstData[k] = int(math.Round(real(coef)))
	}

	return polyMod.NewPolyInt(dstData, mod)
}

func ToMat(pol *polyMod.PolyInt) cMat.CMat {
	srcData := pol.Coefs
	dstData := make([]complex128, len(srcData))
	for k, coef := range srcData {
		dstData[k] = complex(float64(coef), 0)
	}

	return cMat.NewCMat(len(srcData), 1, dstData)
}
