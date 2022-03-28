package random

import (
	"crypto/rand"
	"math"
	"math/big"
	mathrand "math/rand"

	"kazat.ch/lbcrypto/poly"
)

//Returns a pol of deg < N with coefs in [0, Q[
func RandomPol(N int, Q *big.Int) poly.Poly {
	coefs := make([]*big.Int, N)
	for i := 0; i < N; i++ {
		coefs[i], _ = rand.Int(rand.Reader, Q)
	}
	return poly.NewPoly(coefs)
}

//Hamming weight. Returns a slice of lenght N in {-1, 0, 1}^N with h non zero coordinates
func Hwt(N, h int) []*big.Int {

	res := make([]*big.Int, N)

	//remplissage avec des 0 puis avec -1 ou 1
	for i := 0; i < N; i++ {
		if i < N-h {
			res[i] = big.NewInt(0)
		} else {
			rand := mathrand.Intn(2)
			if rand == 0 {
				res[i] = big.NewInt(1)
			} else {
				res[i] = big.NewInt(-1)
			}
		}
	}

	//shake shake shake
	mathrand.Shuffle(N, func(i, j int) { res[i], res[j] = res[j], res[i] })
	return res
}

//returns n samples from ZO(r) with 0 < r < 1
func ZO(n int, r float64) []*big.Int {
	res := make([]*big.Int, n)

	//remplissage avec des -1, 0, 1 suivant la distribution ZO(r)
	for i := 0; i < n; i++ {
		rand := mathrand.Float64()
		if rand < r/2 {
			res[i] = big.NewInt(-1)
		} else if rand < r {
			res[i] = big.NewInt(1)
		} else {
			res[i] = big.NewInt(0)
		}
	}
	return res
}

//returns n samples from DG(sig^2) in a []*big.Int
func DG(n int, s2 float64, max *big.Int) []*big.Int {
	res := make([]*big.Int, n)
	var rand float64
	var randBigInt *big.Int
	min := new(big.Int).Neg(max)

	goodSamples := 0
	for goodSamples < n {

		rand = mathrand.NormFloat64() * math.Sqrt(s2)
		randBigInt = big.NewInt(int64(math.Round(rand)))
		if min.Cmp(randBigInt) == -1 && randBigInt.Cmp(max) == -1 {
			res[goodSamples] = randBigInt
		}
		goodSamples++
	}
	return res
}
