package main

import (
	"fmt"
	"math/rand"
	"testing"
	"time"

	"kazat.ch/lbcrypto/cMat"
	"kazat.ch/lbcrypto/ckks"
	"kazat.ch/lbcrypto/encoder"
	"kazat.ch/lbcrypto/poly"
	"kazat.ch/lbcrypto/ring"
)

func aTestMod(t *testing.T) {
	r := ring.NewRing(13, 4)
	coefs := []float64{1, 20, 33}
	pol := poly.NewPoly(coefs)
	fmt.Println("Pol before:", pol)
	r.TakeCoefMod(&pol)
	fmt.Println("Pol after:", pol)
}

func aTestToRing(t *testing.T) {
	r := ring.NewRing(13, 4)
	coefs := []float64{0, 0, 0, 0, 0, 1}
	pol := poly.NewPoly(coefs)
	r.ToRing(&pol)
	fmt.Println(pol)
}

func aTestSum(t *testing.T) {
	fmt.Println("Testing the sum")
	r := ring.NewRing(13, 4)
	coefsu := []float64{1, 6, 8, 3}
	coefsv := []float64{0, 9, 8, 11}

	u := poly.NewPoly(coefsu)
	v := poly.NewPoly(coefsv)

	s := r.Add(u, v)
	fmt.Println(s)
}

// tests the product in a ring
func aTestProd(t *testing.T) {
	fmt.Println("Testing the product")

	r := ring.NewRing(13, 4)
	coefsu := []float64{1, 6, 8, 3}
	coefsv := []float64{0, 9, 8, 11}

	u := poly.NewPoly(coefsu)
	v := poly.NewPoly(coefsv)

	s := r.Mult(u, v)
	fmt.Println(s)
}

// tests the scaling in a ring
func aTestScaling(t *testing.T) {
	fmt.Println("Testing the scaling")

	r := ring.NewRing(13, 4)
	coefsu := []float64{1, 6, 8, 3}

	u := poly.NewPoly(coefsu)
	fmt.Println(u)
	r.Scale(2, &u)
	fmt.Println(u)
}

// tests the encoder
func aTestEncoding(t *testing.T) {
	c1 := 0.5 + 0.5i
	c2 := 0.5 - 0.5i

	N := 4
	mod := 999983
	scale := 64 + 0i

	//********************************
	enc := encoder.NewEncoder(N, scale, mod)
	va := cMat.NewCMat(2, 1, []complex128{c1, c2})
	fmt.Println("Original vector :")
	va.PPrint()

	ya := enc.Encode(&va)
	fmt.Println("Polynomial encoding of the vector a:")
	fmt.Println(ya)

	xa := enc.Decode(ya)
	fmt.Println("Decoded encoded vector :")
	xa.PPrint()

}

// tests the encoder
func aTestDecoding(t *testing.T) {

	N := 4
	mod := 100
	scale := 64 + 0i

	//********************************
	enc := encoder.NewEncoder(N, scale, mod)

	pol := poly.NewPoly([]float64{1, 1, -10, 50})
	fmt.Println("Poly to decode :", pol)

	xa := enc.Decode(pol)
	fmt.Println("Decoded encoded vector :")
	xa.PPrint()

}

// tests the keygeneration
func aTestKeyGen(t *testing.T) {
	fmt.Println("Testing the key generation")
	rand.Seed(time.Now().UnixNano())
	N := 10
	Q := 5
	h := 3
	ckks := ckks.CKKS{N: N, Q: Q, H: h, Scale: 64 + 0i}

	//********************************
	fmt.Println("GENERATING KEYS")
	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)
	fmt.Println("sk :", sk[1])
	fmt.Println("pk", pk)
}

// tests the encryption
func aTestEncrypt(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	fmt.Println("Testing the encryption")

	mod := 999983
	scale := 64 + 0i
	N := 4
	h := 3
	ckks := ckks.CKKS{N: N, Q: mod, H: h, Scale: scale}
	c1 := 3 + 5i
	c2 := 2 + 1i
	va := cMat.NewCMat(2, 1, []complex128{c1, c2})

	//********************************
	fmt.Println("ENCODING")
	enc := encoder.NewEncoder(ckks.N, ckks.Scale, ckks.Q)
	pta := enc.Encode(&va)

	//********************************
	fmt.Println("GENERATING KEYS")

	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)
	fmt.Println("sk :", sk[1])

	//********************************
	fmt.Println("ENCRYPTING")
	ct := ckks.Encrypt(pta, pk)
	fmt.Println("ct :", ct)
}

// tests the encryption
func aTestEncryptDecrypt(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	fmt.Println("Testing the encryption - decryption")

	mod := 2000000
	scale := 1024 + 0i
	N := 4
	h := 3
	ckks := ckks.CKKS{N: N, Q: mod, H: h, Scale: scale}

	c1 := -350 + 152i
	c2 := 55.13 - 110i

	va := cMat.NewCMat(2, 1, []complex128{c1, c2})
	fmt.Println("message :")
	va.PPrint()

	//********************************
	fmt.Println("ENCODING :")
	enc := encoder.NewEncoder(ckks.N, ckks.Scale, ckks.Q)
	pta := enc.Encode(&va)

	//********************************
	fmt.Println("KEYS GENERATION :")
	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)

	//********************************
	fmt.Println("ENCRYPTING :")
	ct := ckks.Encrypt(pta, pk)
	fmt.Println("ct :", ct[0])

	//********************************
	fmt.Println("DECRYPTING :")
	pt := ckks.Decrypt(ct, sk)
	fmt.Println("pt :", pt)

	//********************************
	fmt.Println("DECODING")
	vect := enc.Decode(pt)
	vect.PPrint()

	//********************************
	fmt.Println("SQUARED NORM OF ERRORS")
	vect = *vect.Scale(-1)
	vect.Add(&vect, &va)
	fmt.Println(vect.SquaredNorm())
}

func TestHomomAdd(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	fmt.Println("Testing homomorphism on +")

	mod := 2000000
	scale := 1024 + 0i
	N := 4
	h := 3
	ckks := ckks.CKKS{N: N, Q: mod, H: h, Scale: scale}

	//********************************

	c1 := -100 + 5i
	c2 := 200 - 50i
	c3 := -35 - 7i
	c4 := 6 - 40i

	va := cMat.NewCMat(2, 1, []complex128{c1, c2})
	vb := cMat.NewCMat(2, 1, []complex128{c3, c4})
	vs := cMat.NewCMat(2, 1, []complex128{0, 0})

	vs.Add(&va, &vb)

	fmt.Println("sum of messages :")
	vs.PPrint()

	//********************************
	fmt.Println("ENCODING")
	enc := encoder.NewEncoder(ckks.N, ckks.Scale, ckks.Q)
	pta := enc.Encode(&va)
	ptb := enc.Encode(&vb)

	//********************************
	fmt.Println("GENERAING KEYS")
	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)

	//********************************
	fmt.Println("ENCRYPTING")
	cta := ckks.Encrypt(pta, pk)
	ctb := ckks.Encrypt(ptb, pk)
	cts := [2]poly.Poly{poly.PolyAdd(cta[0], ctb[0]), poly.PolyAdd(cta[1], ctb[1])}
	fmt.Println("cts :", cts[0])

	//********************************
	fmt.Println("DECRYPTING :")
	pts := ckks.Decrypt(cts, sk)
	fmt.Println("pts :", pts)

	//********************************
	fmt.Println("DECODING")
	vect := enc.Decode(pts)
	vect.PPrint()

	//********************************
	fmt.Println("SQUARED NORM OF ERRORS")
	vect = *vect.Scale(-1)
	vect.Add(&vect, &vs)
	fmt.Println(vect.SquaredNorm())
}
