package encoder

import (
	"fmt"
	"math"
	"math/big"
	"math/cmplx"

	"kazat.ch/lbcrypto/cMat"
	"kazat.ch/lbcrypto/poly"
)

type Encoder struct {
	N        int
	scale    complex128
	ENC      cMat.CMat
	DEC      cMat.CMat
	Basis    cMat.CMat
	CycloPol poly.Poly
}

type PT struct {
	Pol   poly.Poly
	Scale complex128
}

func NewPT(pol poly.Poly, scale complex128) PT {
	PT := PT{Pol: pol, Scale: scale}
	return PT
}

func NewEncoder(N int, scale complex128) Encoder {
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
	coefs := make([]*big.Int, N+1)
	coefs[0], coefs[N] = big.NewInt(1), big.NewInt(1)
	cycloPol := poly.NewPoly(coefs)

	res := Encoder{
		N:        N,
		ENC:      ENC,
		DEC:      DEC,
		scale:    scale,
		Basis:    basis,
		CycloPol: cycloPol,
	}
	return res
}

// Encodes complex vectors of length e.N/2
func (enc *Encoder) Encode(v *cMat.CMat) PT {

	N := enc.N
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

	newv.Scale(enc.scale)
	newv = *cMat.ProjectOnRaws(&newv, &enc.Basis)
	newv.Mult(&enc.ENC, &newv)
	pol := enc.ToPol(&newv)

	res := PT{Pol: pol, Scale: enc.scale}

	return res
}

func (enc *Encoder) Decode(pt PT) cMat.CMat {

	v := enc.ToMat(pt.Pol)
	v.Mult(&enc.DEC, &v)
	v.Scale(1.0 / pt.Scale)
	data := v.GetData()[0 : enc.N/2]
	res := cMat.NewCMat(enc.N/2, 1, data)
	return res
}

func (enc *Encoder) DecodeDebug(pt PT) cMat.CMat {

	v := enc.ToMat(pt.Pol)
	fmt.Println("mat of pt:", v)
	v.Mult(&enc.DEC, &v)
	v.Scale(1.0 / pt.Scale)
	//fmt.Println("scale :", 1./pt.Scale)
	data := v.GetData()[0 : enc.N/2]
	res := cMat.NewCMat(enc.N/2, 1, data)
	return res
}

// returns the poly.Poly corresponding to the vector m
func (enc *Encoder) ToPol(m *cMat.CMat) poly.Poly {
	srcData := m.GetData()
	dstData := make([]*big.Int, len(srcData))
	for k, coef := range srcData {
		dstData[k] = big.NewInt(int64(math.Round(real(coef))))
	}
	pol := poly.NewPoly(dstData)

	return pol
}

//converts pol to a cMat.CMat
func (enc *Encoder) ToMat(pol poly.Poly) cMat.CMat {
	srcData := pol.Coefs
	dstData := make([]complex128, len(srcData))
	for k, coef := range srcData {
		temp := new(big.Float).SetInt(coef)
		coefFloat64, _ := temp.Float64()
		dstData[k] = complex(coefFloat64, 0)

	}

	return cMat.NewCMat(len(srcData), 1, dstData)
}

func (e *Encoder) GetScale() complex128 {
	return e.scale
}

//*************************
//returns the enc-compatible polynomial corresponding to the vector (k, ..., k) at scale scale
func (enc *Encoder) ConstToPT(k float64) PT {
	n := enc.N

	data := make([]complex128, n/2)
	for i := 0; i < n/2; i++ {
		data[i] = complex(k, 0)
	}

	v := cMat.NewCMat(n/2, 1, data)
	pt := enc.Encode(&v)
	return pt
}
