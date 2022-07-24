package ckks

import (
	"math"
	"math/big"

	"kazat.ch/lbcrypto/encoder"
	"kazat.ch/lbcrypto/poly"
	"kazat.ch/lbcrypto/random"
)

// A MODIFIER : IL FAUDRAIT LUI DONNER UN Q_0 ET UN DELTA, LE CKKS SE CALCULE SON Q TOUT SEUL
// Contitent les paramètres du schéma nécessaires à l'exécution des différentes procédures
// Sert essentiellement à éviter l'utilisation de variables globales.
type CKKS struct {
	N     int       // for the X^N + 1
	Q     *big.Int  // initial coef modulus
	P     *big.Int  // for modulus in evaluation key
	Cyclo poly.Poly // the cyclotomic polynomial for the underlying ring (créé à partir de N)
	L     int       //number of levels
	H     int       // for the HWT distribution
	s2    float64   //variance of the DG distribution
}

// Représente un ciphertext (B, A)
// Mélange des notations de l'article original sur CKKS et des autres sources... à modifier
type CT struct {
	A     poly.Poly
	B     poly.Poly
	Mod   *big.Int
	Scale complex128
	L     int // current level
}

// retourne un objet de type ckks contenant tous les paramètres d'une instance du schéma
func NewCKKS(Q, P *big.Int, N, H, L int, s2 float64) CKKS {
	//creating the cyclotomic polynomial
	cyclo := poly.ZeroPoly(N)
	cyclo.Coefs[0].Set(big.NewInt(1))
	cyclo.Coefs[N].Set(big.NewInt(1))

	CKKS := CKKS{N: N, Q: Q, P: P, H: H, Cyclo: cyclo, L: L, s2: s2}
	return CKKS
}

// retourne un ciphertext (b, a) de module, échelle et niveau (mod, scale, L)
func NewCT(a, b poly.Poly, mod *big.Int, scale complex128, L int) CT {

	modulus := new(big.Int).Set(mod) // pour que différents ct crés avec le même mod ne le partagent pas !

	CT := CT{A: a, B: b, Mod: modulus, Scale: scale, L: L}
	return CT
}

// retourne une secret key (1, s)
func (ckks *CKKS) SKeyGen() [2]poly.Poly {

	c1 := poly.NewPoly([]*big.Int{big.NewInt(1)})
	c2 := poly.NewPoly(random.Hwt(ckks.N, ckks.H))

	sk := [2]poly.Poly{c1, c2}
	return sk
}

// retourne une public key (-as + e, a) liée à la sk
func (ckks *CKKS) PKeyGen(sk [2]poly.Poly) [2]poly.Poly {

	//r := ckks.r
	s := sk[1]
	a := random.RandomPol(ckks.N, ckks.Q)
	e := poly.NewPoly(random.DG(ckks.N, ckks.s2, ckks.Q))

	b := poly.MultMod(a, s, ckks.Cyclo)
	b.Scale(big.NewInt(-1))
	b = poly.Add(b, e)

	a.TakeCoefMod(ckks.Q)
	b.TakeCoefMod(ckks.Q)

	pk := [2]poly.Poly{b, a}

	return pk
}

// retourne une evaluation key liée à la sk
func (ckks *CKKS) EvKeyGen(sk [2]poly.Poly) [2]poly.Poly {

	s := sk[1]
	prod := new(big.Int).Mul(ckks.Q, ckks.P)
	a := random.RandomPol(ckks.N, prod)
	e := poly.NewPoly(random.DG(ckks.N, ckks.s2, ckks.Q))

	b := poly.MultMod(a, s, ckks.Cyclo)
	b.Scale(big.NewInt(-1))
	b = poly.Add(b, e)

	s2 := poly.MultMod(s, s, ckks.Cyclo)
	s2.Scale(ckks.P)

	b = poly.Add(b, s2)

	a.TakeCoefMod(prod)
	b.TakeCoefMod(prod)

	evk := [2]poly.Poly{b, a}

	return evk
}

// retourne un ciphertext chiffrant le plaintext pt à l'aide de la clé pk
func (ckks *CKKS) Encrypt(pt encoder.PT, pk [2]poly.Poly) CT {

	v := poly.NewPoly(random.ZO(ckks.N, 0.5))
	e0 := poly.NewPoly(random.DG(ckks.N, ckks.s2, ckks.Q))
	e1 := poly.NewPoly(random.DG(ckks.N, ckks.s2, ckks.Q))

	res0 := poly.MultMod(pk[0], v, ckks.Cyclo)
	res0 = poly.Add(res0, pt.Pol)
	res0 = poly.Add(res0, e0)

	res1 := poly.MultMod(pk[1], v, ckks.Cyclo)
	res1 = poly.Add(res1, e1)

	res0.TakeCoefMod(ckks.Q)
	res1.TakeCoefMod(ckks.Q)

	CT := NewCT(res1, res0, ckks.Q, pt.Scale, ckks.L)

	return CT
}

// retourne un plaintext correspondant au ciphertext ct à l'aide de la clé sk
func (ckks *CKKS) Decrypt(ct CT, sk [2]poly.Poly) encoder.PT {

	pt := poly.MultMod(ct.A, sk[1], ckks.Cyclo)
	pt = poly.Add(pt, ct.B)
	pt.TakeCoefMod(ct.Mod)

	res := encoder.PT{Pol: pt, Scale: ct.Scale}
	return res
}

// returns de sum of cyphertexts ct1 and ct2
// A AJOUTER : VERIFICATION QUE LE SCALING FACTOR EST LE MEME
func (ckks *CKKS) CTAdd(ct1, ct2 CT) CT {

	a := poly.Add(ct1.A, ct2.A)
	b := poly.Add(ct1.B, ct2.B)

	a.TakeCoefMod(ct1.Mod)
	b.TakeCoefMod(ct1.Mod)

	mod := new(big.Int) // nécessaire sinon les mod de ct1 et du prod pointent vers la meme adresse
	mod.Set(ct1.Mod)

	sum := NewCT(a, b, mod, ct1.Scale, ct1.L)

	return sum
}

// return the CT corresponding to the constant vector (k, ..., k) at the given scale
func (ckks *CKKS) ConstToCT(k float64, mod *big.Int, scale complex128, pk [2]poly.Poly) CT {
	enc := encoder.NewEncoder(ckks.N, scale)
	pt := enc.ConstToPT(k)
	ct := ckks.Encrypt(pt, pk)
	res := NewCT(ct.A, ct.B, mod, scale, ckks.L)
	return res
}

//Increases the scale of the reciever by a factor k
//This does not change the underlying message
//AJOUTER ERREUR SI K < 0
func (ct *CT) CTIncScale(k *big.Int) *CT {
	ct.CTScale(k)
	ct.Scale = ct.Scale * complex(float64(k.Int64()), 0)
	return ct
}

// scales the reciever by a integer k. This does not work for non integer k.
// this method changes the underlying message
func (ct *CT) CTScale(k *big.Int) *CT {

	ct.A = poly.Scale(ct.A, k)
	ct.B = poly.Scale(ct.B, k)

	return ct
}

// retourne le produit de ct1 et de ct1 en utilisant l'evaluation key evk
//ON POURRAIT/DEVRAIT CHANGER LE TYPE DE RECIEVER.
//ICI, LE RECIEVER CKKS EST JUSTE UTILISE POUR SON Q, MAIS C EST CELUI DES CT QUI COMPTE !
//DEVRAIT CHECKER QUE LES DEUX CT ONT LE MEME MODULUS
func (ckks *CKKS) CTMult(ct1, ct2 CT, evk [2]poly.Poly) CT {

	d0 := poly.MultMod(ct1.B, ct2.B, ckks.Cyclo)
	d1 := poly.Add(poly.MultMod(ct1.A, ct2.B, ckks.Cyclo), poly.MultMod(ct1.B, ct2.A, ckks.Cyclo))
	d2 := poly.MultMod(ct1.A, ct2.A, ckks.Cyclo)

	d0.TakeCoefMod(ct1.Mod)
	d1.TakeCoefMod(ct1.Mod)
	d2.TakeCoefMod(ct1.Mod)

	res0 := poly.MultMod(d2, evk[0], ckks.Cyclo)
	res0 = poly.ScaleDiv(res0, ckks.P)
	res0 = poly.Add(res0, d0)

	res1 := poly.MultMod(d2, evk[1], ckks.Cyclo)
	res1 = poly.ScaleDiv(res1, ckks.P)
	res1 = poly.Add(res1, d1)

	res0.TakeCoefMod(ct1.Mod)
	res1.TakeCoefMod(ct1.Mod)

	prod := NewCT(res1, res0, ct1.Mod, ct1.Scale*ct2.Scale, ct1.L)
	return prod
}

// returns a CT corresponding to the mean of the ciphertexts in the data list
// computations are done under the pk and evk keys
// ATTENTION : DEVRAIT CHECKER QUE LES DATA ONT LE MEME SCALING FACTOR
func (ckks *CKKS) Mean(data []CT, pk, evk [2]poly.Poly) CT {
	N := ckks.N
	n := len(data)

	mod := new(big.Int) // nécessaire sinon les mod de ct1 et du prod pointent vers la meme adresse
	mod.Set(data[0].Mod)

	res := NewCT(poly.ZeroPoly(ckks.N-1), poly.ZeroPoly(N-1), mod, data[0].Scale, data[0].L)
	for i := 0; i < n; i++ {
		res = ckks.CTAdd(res, data[i])
	}

	k := 1.0 / float64(n)
	scale := complex(math.Log(float64(n))*100000000, 0.) // marche mieux que le scale basé sur data[0]
	ctk := ckks.ConstToCT(k, mod, scale, pk)
	res = ckks.CTMult(ctk, res, evk)
	return res
}

// returns a CT corresponding to the var of the ciphertexts in the data list
// paramètre delta pour l'appel de RS : à modifier : CKKS devrait contenir un champ delta
func (ckks *CKKS) Var(data []CT, pk, evk [2]poly.Poly, delta *big.Int) CT {

	mean := ckks.Mean(data, pk, evk)
	// on ne RS pas ici sinon mean n'aura pas le même modulus que les data

	scalesRatio := big.NewInt(int64(real(mean.Scale / data[0].Scale)))

	//puts the (x - bar(x))^2 into the data variable
	for i := range data {
		data[i].CTIncScale(scalesRatio)
		data[i].CTScale(big.NewInt(-1))
		data[i] = ckks.CTAdd(data[i], mean)
		data[i] = ckks.CTMult(data[i], data[i], evk)
		ckks.RS(&data[i], delta)
	}

	res := ckks.Mean(data, pk, evk)
	ckks.RS(&res, delta)
	return res
}

// returns a rescaled version of ct by factor delta
// DEVRAIT VERIFIER QUE LE SCALE DU CT EST ASSEZ GRAND POUR SUBIR LE RS
func (ckks *CKKS) RS(ct *CT, delta *big.Int) *CT {

	ct.A = poly.ScaleDiv(ct.A, delta)
	ct.B = poly.ScaleDiv(ct.B, delta)

	ct.Mod.Div(ct.Mod, delta)

	ct.A.TakeCoefMod(ct.Mod)
	ct.B.TakeCoefMod(ct.Mod)

	temp := new(big.Float).SetInt(delta)
	ttemp, _ := temp.Float64()
	ct.Scale = ct.Scale / complex(ttemp, 0)

	ct.L = ct.L - 1
	return ct
}
