package ckks

import (
	"math"

	"kazat.ch/lbcrypto/poly"
	"kazat.ch/lbcrypto/random"
)

type CKKS struct {
	N     int // for the X^N + 1
	Q     int // initial coef modulus
	P     int // for modulus in evaluation key
	H     int // for the HWT distribution
	Cyclo poly.Poly
	//Scale complex128 //original scaling factor for the encoder
}

type PT struct {
	Pol   poly.Poly
	Scale complex128
}

type CT struct {
	A     poly.Poly
	B     poly.Poly
	Mod   int
	Scale complex128
}

func NewCKKS(N, Q, P, H int) CKKS {
	cyclo := poly.ZeroPoly(N)
	cyclo.Coefs[0] = 1
	cyclo.Coefs[N] = 1
	CKKS := CKKS{N: N, Q: Q, P: P, H: H, Cyclo: cyclo}
	return CKKS
}

func NewPT(pol poly.Poly, scale complex128) PT {
	PT := PT{Pol: pol, Scale: scale}
	return PT
}

func NewCT(a, b poly.Poly, mod int, scale complex128) CT {
	CT := CT{A: a, B: b, Mod: mod, Scale: scale}
	return CT
}

func (ckks *CKKS) SKeyGen() [2]poly.Poly {

	c1 := poly.NewPoly([]float64{1})
	c2 := poly.NewPoly(random.Hwt(ckks.N, ckks.H))

	sk := [2]poly.Poly{c1, c2}
	return sk
}

func (ckks *CKKS) PKeyGen(sk [2]poly.Poly) [2]poly.Poly {

	//rg := ring.NewRing(ckks.Q, ckks.N)
	s := sk[1]
	a := random.RandomPol(ckks.N, ckks.Q)
	//fmt.Println("a :", a)
	e := poly.NewPoly(random.DG(ckks.N, 0.5))

	b := poly.MultMod(a, s, ckks.Cyclo)
	//fmt.Println("a*s :", b)
	b.Scale(-1)
	b = poly.Add(b, e)

	a.TakeCoefMod(ckks.Q)
	b.TakeCoefMod(ckks.Q)

	pk := [2]poly.Poly{b, a}

	return pk
}

func (ckks *CKKS) EvKeyGen(sk [2]poly.Poly) [2]poly.Poly {

	//rg := ring.NewRing(ckks.Q*ckks.P, ckks.N)
	s := sk[1]
	//fmt.Println("s :", s)
	a := random.RandomPol(ckks.N, ckks.Q*ckks.P)
	//fmt.Println("a :", a)
	e := poly.NewPoly(random.DG(ckks.N, 0.5))

	b := poly.MultMod(a, s, ckks.Cyclo)
	//fmt.Printf("a*s : %0.f \n", b)

	b.Scale(-1)
	//fmt.Printf("-a*s : %0.f \n", b)

	//fmt.Printf("e : %0.f \n", e)
	b = poly.Add(b, e)
	//fmt.Printf("-as + e : %0.f \n", b)

	s2 := poly.MultMod(s, s, ckks.Cyclo)
	//fmt.Printf("s2 : %0.f \n", s2)
	s2.Scale(float64(ckks.P))
	//fmt.Printf("P* s2 : %0.f \n", s2)

	b = poly.Add(b, s2)

	//fmt.Printf("-as + e + P * s2 : %0.f \n", b)

	a.TakeCoefMod(ckks.Q * ckks.P)
	b.TakeCoefMod(ckks.Q * ckks.P)
	evk := [2]poly.Poly{b, a}

	return evk
}

func (ckks *CKKS) Encrypt(pt PT, pk [2]poly.Poly) CT {

	//rg := ring.NewRing(ckks.Q, ckks.N)

	v := poly.NewPoly(random.ZO(ckks.N, 0.5))
	e0 := poly.NewPoly(random.DG(ckks.N, 0.5))
	e1 := poly.NewPoly(random.DG(ckks.N, 0.5))

	res0 := poly.MultMod(pk[0], v, ckks.Cyclo)
	res0 = poly.Add(res0, pt.Pol)
	res0 = poly.Add(res0, e0)

	res1 := poly.MultMod(pk[1], v, ckks.Cyclo)
	res1 = poly.Add(res1, e1)

	res0.TakeCoefMod(ckks.Q)
	res1.TakeCoefMod(ckks.Q)

	CT := NewCT(res1, res0, ckks.Q, pt.Scale)

	//return [2]poly.Poly{res0, res1}
	//fmt.Println("check :", res0, res1)
	return CT
}

func (ckks *CKKS) Decrypt(ct CT, sk [2]poly.Poly) PT {

	//rg := ring.NewRing(ckks.Q, ckks.N)

	pt := poly.MultMod(ct.A, sk[1], ckks.Cyclo)
	pt = poly.Add(pt, ct.B)
	pt.TakeCoefMod(ckks.Q)

	res := PT{Pol: pt, Scale: ct.Scale}
	return res
}

// returns de sum of cyphertexts ct1 and ct2
func (ckks *CKKS) CTAdd(ct1, ct2 CT) CT {
	//rg := ring.NewRing(ckks.Q, ckks.N)

	a := poly.Add(ct1.A, ct2.A)
	b := poly.Add(ct1.B, ct2.B)

	a.TakeCoefMod(ckks.Q)
	b.TakeCoefMod(ckks.Q)

	sum := NewCT(a, b, ckks.Q, ct1.Scale)

	return sum
}

// adds b to a
func (a *PT) PTAdd(b PT) *PT {

	a.Pol = poly.Add(a.Pol, b.Pol)

	return a
}

// scales a by a factor k used with a precision of nbd decimals
// returns an approximate encoding of k*dec(a)

func (a *PT) PTScale(k float64, nbd int) *PT {

	k = math.Round(k * math.Pow(10, float64(nbd)))
	a.Pol = poly.Scale(a.Pol, k)
	a.Scale = a.Scale * complex(math.Pow(10, float64(nbd)), 0)

	return a
}

// returns the product of pt1 and pt2 as PT
// BOF... LE FAIT D'APPELER SUR UN CKKS EST PAS TRES INTUITIF... MAIS SINON FAUT SE DONNER EXPLICITEMENT LE CYCLO
func (ckks *CKKS) PTMult(pt1, pt2 PT) PT {

	pol := poly.MultMod(pt1.Pol, pt2.Pol, ckks.Cyclo)
	scale := pt1.Scale * pt2.Scale
	res := NewPT(pol, scale)

	return res
}

// returns the scaling of ct by k
//func (ckks *CKKS) CTScale(ct CT, ptk PT, pk, evk [2]poly.Poly) CT {
//	ctk := ckks.Encrypt(ptk, pk)
//	return ckks.CTMult(ct, ctk, evk)
//}

// scales ct by a real k, with a precision of nbd decimals
func (ckks *CKKS) CTScale(ct *CT, k float64, nbd int) *CT {
	k = math.Round(k * math.Pow(10, float64(nbd)))

	ct.A = poly.Scale(ct.A, k)
	ct.B = poly.Scale(ct.B, k)
	ct.A.TakeCoefMod(ckks.Q) // bon... si on dépasse le modulo, on est de tt façon ko... pour la beauté du geste seulement
	ct.B.TakeCoefMod(ckks.Q)

	ct.Scale = ct.Scale * complex(math.Pow(10, float64(nbd)), 0)

	return ct
}

func (ckks *CKKS) CTMult(ct1, ct2 CT, evk [2]poly.Poly) CT {

	//rg := ring.NewRing(ckks.Q, ckks.N)

	d0 := poly.MultMod(ct1.B, ct2.B, ckks.Cyclo)
	d1 := poly.Add(poly.MultMod(ct1.A, ct2.B, ckks.Cyclo), poly.MultMod(ct1.B, ct2.A, ckks.Cyclo))
	d2 := poly.MultMod(ct1.A, ct2.A, ckks.Cyclo)

	d0.TakeCoefMod(ckks.Q)
	d1.TakeCoefMod(ckks.Q)
	d2.TakeCoefMod(ckks.Q)

	res0 := poly.MultMod(d2, evk[0], ckks.Cyclo)
	//fmt.Printf("res0 : %0.f \n", res0)
	res0.Scale(1 / float64(ckks.P))
	//fmt.Printf("res0 : %f \n", res0)
	res0.Round()
	//fmt.Printf("rounded res0 : %1.f \n", res0)
	res0 = poly.Add(res0, d0)
	//rg.ToRing(&res0)

	res1 := poly.MultMod(d2, evk[1], ckks.Cyclo)
	res1.Scale(1 / float64(ckks.P))
	res1.Round()
	res1 = poly.Add(res1, d1)

	res0.TakeCoefMod(ckks.Q)
	res1.TakeCoefMod(ckks.Q)

	prod := NewCT(res1, res0, ckks.Q, ct1.Scale*ct2.Scale)

	return prod
}

// returns the mean of a list of CTs as a CT
// ATTENTION : DEVRAIT CHECKER QUE LES DATA ONT LE MEME SCALING FACTOR
func (ckks *CKKS) Mean(data []CT) CT {
	N := ckks.N
	n := len(data)
	//scale := data[0].Scale

	res := NewCT(poly.ZeroPoly(ckks.N-1), poly.ZeroPoly(N-1), data[0].Mod, data[0].Scale)
	for i := 0; i < n; i++ {
		res = ckks.CTAdd(res, data[i])
	}

	k := 1. / float64(n)
	ckks.CTScale(&res, k, int(math.Ceil(math.Log10(float64(n)))))
	return res
}

// TO DO :
// - SUM FOR DIFFERENT SCALES
