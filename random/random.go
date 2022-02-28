package random

import (
	"math"
	"math/rand"

	"kazat.ch/lbcrypto/poly"
)

//Returns a pol of deg < N with coefs in [0, Q[
func RandomPol(N, Q int) poly.Poly {
	coefs := make([]float64, N)
	for i := 0; i < N; i++ {
		coefs[i] = float64(rand.Intn(Q))
	}
	return poly.NewPoly(coefs)
}

//Hamming weight. Returns a slice of lenght N in {-1, 0, 1}^N with h non zero coordinates
func Hwt(N, h int) []float64 {

	res := make([]float64, N)

	//remplissage avec des 0 puis avec -1 ou 1
	for i := 0; i < N; i++ {
		if i < N-h {
			res[i] = 0
		} else {
			rand := rand.Intn(2) + 1
			res[i] = math.Pow(-1, float64(rand))
		}
	}

	//shake shake shake
	rand.Shuffle(N, func(i, j int) { res[i], res[j] = res[j], res[i] })
	return res
}

//returns n samples from ZO(r) with 0 < r < 1
func ZO(n int, r float64) []float64 {
	res := make([]float64, n)

	//remplissage avec des -1, 0, 1 suivant la distribution ZO(r)
	for i := 0; i < n; i++ {
		rand := rand.Float64()
		if rand < r/2 {
			res[i] = -1
		} else if rand < r {
			res[i] = 1
		} else {
			res[i] = 0
		}
	}
	return res
}

// TO BE IMPLEMENTED. FOR NOW, CALLS ZO(sig^2)
//returns n samples from DG(sig^2)
// TO BE IMPLEMENTED. FOR NOW, CALLS ZO(sig^2)

func DG(n int, s2 float64) []float64 {
	r := s2
	return ZO(n, r)
}
