package encoder

import (
	"math"
	"math/cmplx"

	"kazat.ch/lbcrypto/cMat"
	"kazat.ch/lbcrypto/poly"
	"kazat.ch/lbcrypto/ring"
)

type Encoder struct {
	N        int
	mod      int
	scale    complex128
	ENC      cMat.CMat
	DEC      cMat.CMat
	Basis    cMat.CMat
	CycloPol poly.Poly
}

func (e *Encoder) GetMod() int {
	return e.mod
}
func NewEncoder(N int, scale complex128, mod int) Encoder {
	if N%2 != 0 {
		panic("Error : parameter N must be even")
	}
	exposant := complex(0, math.Pi/float64(N))
	xi := cmplx.Exp(exposant)
	v := make([]complex128, int(N*N))

	// Création de la matrice de Vandermonde et de son inverse
	// i.e. isomorphisme sigma et sigma inverse
	for i := 0; i < N; i++ {
		root := cmplx.Pow(xi, complex(float64(2*i)+1, 0))
		for j := 0; j < N; j++ {
			v[N*i+j] = cmplx.Pow(root, complex(float64(j), 0))
		}
	}
	DEC := cMat.NewCMat(N, N, v)
	R, S := DEC.Invert()
	ENC := cMat.Merge(&R, &S)

	// Image de la base canonique de l'ensemble des polynomes
	// par sigma (pour projection des vecteurs)
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

	// Création du Mème polynôme cyclotomique pour M = 2^p
	coefs := make([]float64, N+1)
	coefs[0], coefs[N] = 1, 1
	cycloPol := poly.NewPoly(coefs)

	res := Encoder{
		N:        N,
		mod:      mod,
		ENC:      ENC,
		DEC:      DEC,
		scale:    scale,
		Basis:    basis,
		CycloPol: cycloPol,
	}
	return res
}

// Encodes complex vectors of length e.N/2
func (e *Encoder) Encode(v *cMat.CMat) poly.Poly {

	N := e.N
	originaldata := v.GetData()
	if len(originaldata) != N/2 {
		panic("Error : message dimension does not match this encoder")
	}

	data := make([]complex128, N)

	for i, val := range originaldata {
		data[i] = val
		data[i+N/2] = cmplx.Conj(originaldata[N/2-1-i])
	}

	newv := cMat.NewCMat(N, 1, data)

	newv.Scale(e.scale)
	newv = *cMat.ProjectOnRaws(&newv, &e.Basis)
	newv.Mult(&e.ENC, &newv)
	pol := e.ToPol(&newv)

	return pol
}

func (e *Encoder) Decode(pol poly.Poly) cMat.CMat {

	v := e.ToMat(pol)
	v.Mult(&e.DEC, &v)
	v.Scale(1 / e.scale)
	data := v.GetData()[0 : e.N/2]
	res := cMat.NewCMat(e.N/2, 1, data)
	return res
}

// returns the ring element corresponding to the vector m
// ring parameters are deduced from enc
func (enc *Encoder) ToPol(m *cMat.CMat) poly.Poly {
	ring := ring.NewRing(enc.mod, enc.N)
	srcData := m.GetData()
	dstData := make([]float64, len(srcData))
	for k, coef := range srcData {
		dstData[k] = math.Round(real(coef))
	}
	//fmt.Println("dstData : ", dstData)
	pol := poly.NewPoly(dstData)
	//fmt.Println("pol before to ring :", pol)
	ring.ToRing(&pol)
	//fmt.Println("pol after to ring :", pol)

	return pol
}

func round(f float64) {
	panic("unimplemented")
}

func (enc *Encoder) ToMat(pol poly.Poly) cMat.CMat {
	srcData := pol.Coefs
	dstData := make([]complex128, len(srcData))
	for k, coef := range srcData {
		dstData[k] = complex(float64(coef), 0)
	}

	return cMat.NewCMat(len(srcData), 1, dstData)
}

func (e *Encoder) GetScale() complex128 {
	return e.scale
}

//*************************
