package ckks

import (
	"kazat.ch/lbcrypto/poly"
	"kazat.ch/lbcrypto/random"
	"kazat.ch/lbcrypto/ring"
)

type CKKS struct {
	N     int        // for the X^N + 1
	Q     int        // initial coef modulus
	H     int        // for the HWT distribution
	Scale complex128 //original scaling factor for the encoder
}

func (ckks *CKKS) SKeyGen() [2]poly.Poly {

	c1 := poly.NewPoly([]float64{1})
	c2 := poly.NewPoly(random.Hwt(ckks.N, ckks.H))

	sk := [2]poly.Poly{c1, c2}
	return sk
}

func (ckks *CKKS) PKeyGen(sk [2]poly.Poly) [2]poly.Poly {

	rg := ring.NewRing(ckks.Q, ckks.N)
	a := random.RandomPol(ckks.N, ckks.Q)
	//fmt.Println("sk :", sk[1])

	//fmt.Println("a :", a)

	b := rg.Mult(a, sk[1])
	//fmt.Println("a*s :", b)
	rg.Scale(-1, &b)
	//fmt.Println("-a*s :", b)

	e := poly.NewPoly(random.ZO(ckks.N, 0.5))
	//fmt.Println("e :", e)
	b = rg.Add(b, e)
	//fmt.Println("b :", b)
	pk := [2]poly.Poly{b, a}

	return pk
}

func (ckks *CKKS) Encrypt(pt poly.Poly, pk [2]poly.Poly) [2]poly.Poly {
	rg := ring.NewRing(ckks.Q, ckks.N)

	//fmt.Println("pt :", pt)

	v := poly.NewPoly(random.ZO(ckks.N, 0.5))
	//rg.TakeCoefMod(&v)
	e0 := poly.NewPoly(random.ZO(ckks.N, 0.5))
	//rg.TakeCoefMod(&e0)
	e1 := poly.NewPoly(random.ZO(ckks.N, 0.5))
	//rg.TakeCoefMod(&e1)
	//fmt.Println("ENCRYPTION")
	//fmt.Println("PK :", pk)
	//fmt.Println("v: ", v)
	//fmt.Println("b =", pk[0])

	res0 := rg.Mult(pk[0], v)
	//fmt.Println("v*b =", res0)

	res0 = rg.Add(res0, pt)
	//fmt.Println("v*b + pt =", res0)

	res0 = rg.Add(res0, e0)
	//fmt.Println("e0 =", e0)

	//fmt.Println("v*b + pt + e0 =", res0)

	res1 := rg.Mult(pk[1], v)
	res1 = rg.Add(res1, e1)

	rg.TakeCoefMod(&res0)
	rg.TakeCoefMod(&res1)

	return [2]poly.Poly{res0, res1}
}

func (ckks *CKKS) Decrypt(c [2]poly.Poly, sk [2]poly.Poly) poly.Poly {

	rg := ring.NewRing(ckks.Q, ckks.N)
	//fmt.Println("c1 :", c[1])
	//fmt.Println("sk1 :", sk[1])
	pta := rg.Mult(c[1], sk[1])
	//fmt.Println("pta = c1 * sk1 :", pta)
	//fmt.Println("c0 :", c[0])
	ptb := rg.Add(pta, c[0])
	//fmt.Println("pta + c0 :", ptb)

	//fmt.Println(ptb)
	rg.TakeCoefMod(&ptb)

	return ptb
}
