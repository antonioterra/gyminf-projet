package ring

import (
	"kazat.ch/lbcrypto/poly"
)

// RING ELEMENT
type Ring struct {
	CMod int
	PMod poly.Poly
}

// returns a ring with coef modulus mod and pol modulus X^N + 1
func NewRing(mod, N int) Ring {
	pol := poly.ZeroPoly(N)
	pol.Coefs[0] = 1
	pol.Coefs[N] = 1
	return Ring{CMod: mod, PMod: pol}
}

// takes the mod of the coefs of pol (in place)
func (r *Ring) TakeCoefMod(pol *poly.Poly) {
	coefs := pol.Coefs
	for i, c := range coefs {
		x := int(c) % r.CMod
		if x < 0 {
			x = x + r.CMod
		}
		if x > r.CMod/2 {
			x = x - r.CMod
		}
		coefs[i] = float64(x)
	}
}

// changes pol to [pol % PMod]_q (in place)
func (r *Ring) ToRing(pol *poly.Poly) {

	_, reste := poly.PolyDiv(*pol, r.PMod)
	//r.TakeCoefMod(&reste)
	pol.Coefs = reste.Coefs
}

// returns the ring-sum of ring elements u and v
// u and v MUST already be ring elements
func (r *Ring) Add(u, v poly.Poly) poly.Poly {
	res := poly.PolyAdd(u, v)
	//r.TakeCoefMod(&res)
	return res
}

// returns the ring-product of ring elements u and v
// u and v MUST already be ring elements
func (r *Ring) Mult(u, v poly.Poly) poly.Poly {
	res := poly.Mult(u, v)
	//fmt.Println("before :", res)
	r.ToRing(&res)
	//fmt.Println("after :", res)

	res.Deflate()
	return res
}

// changes pol to k*pol as ring element (in place)
func (r *Ring) Scale(k float64, pol *poly.Poly) {
	coefs := pol.Coefs
	for i, c := range coefs {
		coefs[i] = k * c
	}
	r.ToRing(pol)
}
