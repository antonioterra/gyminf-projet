package encoder

import (
	"math"
	"math/big"
	"math/cmplx"

	"kazat.ch/lbcrypto/cMat"
	"kazat.ch/lbcrypto/poly"
)

// Contitent les paramètres du schéma nécessaires à l'exécution des procédures d'encodag et décodage
// Sert essentiellement à éviter l'utilisation de variables globales.
type Encoder struct {
	N     int        // double de la dimension des vecteurs à encoder
	scale complex128 // mise à l'échelle effectuée durant l'encodage
	ENC   cMat.CMat  // matrice d'encodage
	DEC   cMat.CMat  // matrice de décodage
	Basis cMat.CMat  // base pour l'espace sur lequel on projète durant l'encodage
}

// Décrit un plaintext avec les infos nécessaires au décodage dans Scale
type PT struct {
	Pol   poly.Poly
	Scale complex128
}

// revoie le plaintext formé du polynôme <pol> et associé à l'échelle <scale>
func NewPT(pol poly.Poly, scale complex128) PT {
	PT := PT{Pol: pol, Scale: scale}
	return PT
}

// renvoie un encoder complet sur la base du double de la dimension des vecteurs à encoder <N>
// et de l'échelle de base faite durant l'encodage <scale>
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

	res := Encoder{
		N:     N,
		ENC:   ENC,
		DEC:   DEC,
		scale: scale,
		Basis: basis,
	}
	return res
}

// Renvoie un plaintext correspondant à l'encodage d'un vecteur <v> à l'aide du reciever <enc>
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

// Renvoie un vecteur correspondant au décodage d'un plaintext <pt> à l'aide du reciever <enc>
func (enc *Encoder) Decode(pt PT) cMat.CMat {

	v := enc.ToMat(pt.Pol)
	v.Mult(&enc.DEC, &v)
	v.Scale(1.0 / pt.Scale)
	data := v.GetData()[0 : enc.N/2]
	res := cMat.NewCMat(enc.N/2, 1, data)
	return res
}

// retourne l'objet de type poly.Poly dont les coefficiants sont donnés par le vecteur <m>
func (enc *Encoder) ToPol(m *cMat.CMat) poly.Poly {
	srcData := m.GetData()
	dstData := make([]*big.Int, len(srcData))
	for k, coef := range srcData {
		dstData[k] = big.NewInt(int64(math.Round(real(coef))))
	}
	pol := poly.NewPoly(dstData)

	return pol
}

// retourne l'objet de type cMat dont les coefficiants sont donnés par le polynôme <pol>
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

//*************************
//returns the <enc>-compatible polynomial corresponding to the vector (k, ..., k) at scale <scale>
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
