package main

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"

	"kazat.ch/lbcrypto/cMat"
	"kazat.ch/lbcrypto/ckks"
	"kazat.ch/lbcrypto/encoder"
	"kazat.ch/lbcrypto/poly"
	"kazat.ch/lbcrypto/ring"
)

// return a random complex128 with real et imaginary parts in (-b, b)
func randComplex(b float64) complex128 {
	//rand.Seed(time.Now().UnixNano())

	re := math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64()
	im := math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64()
	cpx := complex(re, im)
	return cpx
}

// return a random complex128 with real et imaginary parts in (-b, b)
func randIntComplex(b float64) complex128 {
	//rand.Seed(time.Now().UnixNano())
	re := math.Round(math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64())
	im := math.Round(math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64())
	cpx := complex(re, im)
	return cpx
}

// return a vector of n complex128 with parts in (-b, b) as []complex128
func randComplexVect(n int, b float64) []complex128 {

	//rand.Seed(time.Now().UnixNano())
	vect := make([]complex128, n)
	for i, _ := range vect {
		vect[i] = randComplex(b)
	}
	return vect
}

// return a vector of n random course means (as complex numbers)
func randGradesVect(n int) []complex128 {

	rand.Seed(time.Now().UnixNano())
	vect := make([]complex128, n)
	for i, _ := range vect {
		vect[i] = complex(float64(rand.Intn(11)/2+1), 0)
	}
	return vect
}

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

	s := poly.Add(u, v)
	r.TakeCoefMod(&s)
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
	u.Scale(2)
	r.TakeCoefMod(&u)
	fmt.Println(u)
}

// tests the encoder
func aTestEncoding(t *testing.T) {

	N := 4
	mod := 999983
	scale := 64 + 0i

	//********************************
	enc := encoder.NewEncoder(N, scale, mod)
	va := cMat.NewCMat(2, 1, randComplexVect(2, 100))
	va.PPrint()

	ya := enc.Encode(&va)
	fmt.Println("Polynomial encoding of the vector a:")
	fmt.Println(ya)

	xa := enc.Decode(ya)
	fmt.Println("Decoded encoded vector :")
	xa.PPrint()

}

func aTestHomomEncoding(t *testing.T) {
	rand.Seed(time.Now().UnixNano())

	Q := 999983
	baseScale := 64 + 0i
	N := 4
	h := 3
	P := 10
	ckks := ckks.NewCKKS(N, Q, P, h)

	//********************************
	enc := encoder.NewEncoder(N, baseScale, Q)
	va := cMat.NewCMat(2, 1, randComplexVect(2, 10))
	vb := cMat.NewCMat(2, 1, randComplexVect(2, 10))
	vprod := cMat.NewCMat(2, 1, randComplexVect(2, 2))
	vprod.CoefWiseProd(&va, &vb)

	vprod.PPrint()

	pta := enc.Encode(&va)
	ptb := enc.Encode(&vb)

	ptprod := ckks.PTMult(pta, ptb)
	decoded := enc.Decode(ptprod)
	decoded.PPrint()

}

// tests the encoder
func aTestDecoding(t *testing.T) {

	N := 4
	mod := 100
	baseScale := 64 + 0i

	//********************************
	enc := encoder.NewEncoder(N, baseScale, mod)

	pol := poly.NewPoly([]float64{1, 1, -10, 50})
	pt := ckks.PT{Pol: pol, Scale: baseScale}
	fmt.Println("Poly to decode :", pol)

	xa := enc.Decode(pt)
	fmt.Println("Decoded encoded vector :")
	xa.PPrint()

}

// tests the encryption
func aTestEncrypt(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	fmt.Println("Testing the encryption")

	Q := 999983
	baseScale := 64 + 0i
	N := 4
	h := 3
	P := 10
	ckks := ckks.NewCKKS(N, Q, P, h)

	va := cMat.NewCMat(N/2, 1, randComplexVect(N/2, 100))

	//********************************
	fmt.Println("ENCODING")
	enc := encoder.NewEncoder(ckks.N, baseScale, ckks.Q)
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

	Q := 2000000
	baseScale := 1024 + 0i
	N := 4
	h := 3
	P := 10
	ckks := ckks.NewCKKS(N, Q, P, h)

	va := cMat.NewCMat(N/2, 1, randComplexVect(N/2, 100))

	fmt.Println("message :")
	va.PPrint()

	//********************************
	fmt.Println("ENCODING :")
	enc := encoder.NewEncoder(ckks.N, baseScale, ckks.Q)
	pta := enc.Encode(&va)

	//********************************
	fmt.Println("KEYS GENERATION :")
	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)

	//********************************
	fmt.Println("ENCRYPTING :")
	ct := ckks.Encrypt(pta, pk)
	fmt.Println("ct :", ct.B)

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
	fmt.Printf("%f \n", real(vect.SquaredNorm()))
}

func aTestHomomAdd(t *testing.T) {
	for i := 0; i < 1000; i = i + 1 {
		rand.Seed(time.Now().UnixNano())
		//fmt.Println("Testing homomorphism on +")

		Q := 200000

		baseScale := 64 + 0i
		N := 16
		h := 3
		P := 10
		ckks1 := ckks.NewCKKS(N, Q, P, h)

		//********************************
		v1 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, 10))
		v2 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, 10))
		vs := cMat.NewCMat(N/2, 1, make([]complex128, N/2))

		vs.Add(&v1, &v2)

		//fmt.Println("sum of messages :")
		//vs.PPrint()

		//********************************
		//fmt.Println("ENCODING")
		enc := encoder.NewEncoder(ckks1.N, baseScale, ckks1.Q)
		pt1 := enc.Encode(&v1)
		pt2 := enc.Encode(&v2)
		//fmt.Println("pts :", enc.Encode(&vs))

		//********************************
		//fmt.Println("GENERAING KEYS")
		sk := ckks1.SKeyGen()
		//fmt.Println("sk :", sk)
		pk := ckks1.PKeyGen(sk)
		//fmt.Println("pk :", pk)

		//********************************
		//fmt.Println("ENCRYPTING")
		ct1 := ckks1.Encrypt(pt1, pk)
		ct2 := ckks1.Encrypt(pt2, pk)

		//fmt.Println("ct1 :", ct1)
		//fmt.Println("ct2 :", ct2)

		cts := ckks1.CTAdd(ct1, ct2)
		//fmt.Println("cts :", cts)

		//********************************
		//fmt.Println("DECRYPTING :")
		pts := ckks1.Decrypt(cts, sk)
		//fmt.Println("pts :", pts)

		//********************************
		//fmt.Println("DECODING")
		vect := enc.Decode(pts)
		//vect.PPrint()

		//********************************
		//fmt.Println("SQUARED NORM OF ERRORS")
		vect = *vect.Scale(-1)
		vect.Add(&vect, &vs)
		err := real(vect.SquaredNorm())
		if err > 0.01 {
			fmt.Printf("%f \n", err)
		}
	}
}

// tests the keygeneration
func aTestKeyGen(t *testing.T) {
	fmt.Println("Testing the key generation")
	rand.Seed(time.Now().UnixNano())
	N := 4
	Q := 5
	P := 10
	h := 3
	ckks := ckks.NewCKKS(N, Q, P, h)

	//********************************
	fmt.Println("GENERATING KEYS")
	sk := ckks.SKeyGen()
	fmt.Println("sk :", sk)

	pk := ckks.PKeyGen(sk)
	fmt.Println("pk", pk)

	evk := ckks.EvKeyGen(sk)
	fmt.Println("evk", evk)

}

func aTestHomomMult(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	fmt.Println("Testing homomorphism on *")

	Q := int(2e9)         //200'000 // 2'000'000
	baseScale := 128 + 0i //64 / 128
	N := 8
	h := 3
	P := Q

	inputVectorBound := 5

	ckks1 := ckks.NewCKKS(N, Q, P, h)

	//********************************
	v1 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, float64(inputVectorBound)))
	v2 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, float64(inputVectorBound)))
	vp := cMat.NewCMat(N/2, 1, make([]complex128, N/2))

	v1.PPrint()
	v2.PPrint()

	vp.CoefWiseProd(&v1, &v2)

	fmt.Println("coefwise product of messages :")
	//v1.PPrint()
	//v2.PPrint()
	vp.PPrint()

	//********************************
	fmt.Println("ENCODING")
	enc := encoder.NewEncoder(ckks1.N, baseScale, ckks1.Q)
	pt1 := enc.Encode(&v1)
	pt2 := enc.Encode(&v2)
	fmt.Println("ptp :", enc.Encode(&vp))

	//********************************
	fmt.Println("GENERAING KEYS")
	sk := ckks1.SKeyGen()
	pk := ckks1.PKeyGen(sk)
	evk := ckks1.EvKeyGen(sk)

	//********************************
	fmt.Println("ENCRYPTING")
	ct1 := ckks1.Encrypt(pt1, pk)
	//fmt.Println("cts :", ct1)

	ct2 := ckks1.Encrypt(pt2, pk)
	//fmt.Println("cts :", ct2)

	fmt.Println("ct1 : ", ct1)
	fmt.Println("ct2 : ", ct2)

	//********************************
	fmt.Println("PRODUCT OF CYPHERTEXTS")
	ctp := ckks1.CTMult(ct1, ct2, evk)
	fmt.Println("ctp :", ctp)

	//********************************
	fmt.Println("DECRYPTING PRODUCT:")
	ptp := ckks1.Decrypt(ctp, sk)
	fmt.Println("pt :", ptp)

	//********************************
	fmt.Println("DECODING")
	enc2 := encoder.NewEncoder(ckks1.N, baseScale, ckks1.Q)

	vect := enc2.Decode(ptp)
	vect.PPrint()

	//********************************
	fmt.Println("SQUARED NORM OF ERRORS")
	vect = *vect.Scale(-1)
	vect.Add(&vect, &vp)
	fmt.Printf("%f \n", real(vect.SquaredNorm()))
}

func aTestConstant(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	Q := 200000

	baseScale := 1 + 0i
	N := 4
	enc := encoder.NewEncoder(N, baseScale, Q)

	k := 4.0
	pt := enc.ConstToPT(k, baseScale)
	fmt.Println(pt)
}

func aTestHomomScaling(t *testing.T) {
	fmt.Println("TEST HOMOM SCALING")
	rand.Seed(time.Now().UnixNano())

	Q := 200000000

	baseScale := 64 + 0i
	N := 4
	h := 3
	P := Q
	ckks1 := ckks.NewCKKS(N, Q, P, h)
	enc := encoder.NewEncoder(ckks1.N, baseScale, ckks1.Q)

	k := 0.52735
	//********************************
	v := cMat.NewCMat(N/2, 1, randComplexVect(N/2, 10))
	vscaled := v.Copy()
	vscaled.Scale(complex(k, 0))
	v.PPrint()
	vscaled.PPrint()

	//********************************
	fmt.Println("GENERAING KEYS")
	sk := ckks1.SKeyGen()
	pk := ckks1.PKeyGen(sk)
	//evk := ckks1.EvKeyGen(sk)
	//********************************
	fmt.Println("ENCODING")
	pt := enc.Encode(&v)
	scaledpt := enc.Encode(&vscaled)
	fmt.Println("pt :", pt)
	fmt.Println("pt of scaled vect", scaledpt)

	//********************************
	fmt.Println("ENCRYPTING")
	ct := ckks1.Encrypt(pt, pk)
	ctok := ckks1.Encrypt(scaledpt, pk)
	//********************************
	fmt.Println("PERFORMING COMPUTATIONS")
	fmt.Println("ct before :", ct)
	//ct.A.Scale(k)
	//ct.B.Scale(k)
	ckks1.CTScale(&ct, k, 5)
	fmt.Println("ct after :", ct)
	//********************************
	fmt.Println("DECRYPTING :")
	pt2 := ckks1.Decrypt(ct, sk)
	fmt.Println("pt :", pt2)
	ptok := ckks1.Decrypt(ctok, sk)
	fmt.Println("should be :", ptok)
	//********************************
	fmt.Println("DECODING")
	vect := enc.Decode(pt2)
	vect.PPrint()

	//********************************
	fmt.Println("SQUARED NORM OF ERRORS")
	vect = *vect.Scale(-1)
	vect.Add(&vect, &vscaled)
	err := real(vect.SquaredNorm())
	fmt.Printf("%f \n", err)

}

func aTestPTScale(t *testing.T) {
	rand.Seed(time.Now().UnixNano())

	Q := 999983
	baseScale := 128 + 0i
	N := 8

	//********************************
	enc := encoder.NewEncoder(N, baseScale, Q)
	v := cMat.NewCMat(N/2, 1, randComplexVect(N/2, 10))
	pt := enc.Encode(&v)

	k := 1. / 3
	v.Scale(complex(k, 0))
	fmt.Println("scaled v :")
	v.PPrint()

	fmt.Println("pt = enc(v) :", pt)

	pt.PTScale(k, 2)
	fmt.Println("scaled pt :", pt)

	dec := enc.Decode(pt)
	fmt.Println("dec(k*pt) :")
	dec.PPrint()

}
func TestMean(t *testing.T) {
	rand.Seed(time.Now().UnixNano())

	nbStudents := 200
	nbCourses := 15

	Q := 2000000

	baseScale := 512 + 0i
	N := 2 * nbStudents
	h := 3
	P := Q

	ckks1 := ckks.NewCKKS(N, Q, P, h)

	//********************************

	fmt.Println("Number of students :", nbStudents)
	fmt.Println("Number of courses :", nbCourses)

	vectors := make([]cMat.CMat, nbCourses)
	mean := cMat.NewCMat(N/2, 1, randComplexVect(N/2, 0))
	for i := 0; i < nbCourses; i++ {
		vectors[i] = cMat.NewCMat(nbStudents, 1, randGradesVect(nbStudents))
		mean.Add(&mean, &vectors[i])
		fmt.Println("Means for course n° ", i)
		vectors[i].PPrint()
		fmt.Println("")
	}

	mean.Scale(1 / complex(float64(nbCourses), 0))
	fmt.Println("Mean :")
	mean.PPrint()

	//********************************
	fmt.Println("ENCODING-ENCODING DATA")

	enc := encoder.NewEncoder(ckks1.N, baseScale, ckks1.Q)
	pt := make([]ckks.PT, nbCourses)
	for i := 0; i < nbCourses; i++ {
		pt[i] = enc.Encode(&vectors[i])
	}

	//********************************
	sk := ckks1.SKeyGen()
	pk := ckks1.PKeyGen(sk)
	//evk := ckks1.EvKeyGen(sk)
	//********************************
	ct := make([]ckks.CT, nbCourses)

	for i := 0; i < nbCourses; i++ {
		ct[i] = ckks1.Encrypt(pt[i], pk)
	}
	//********************************
	fmt.Println("PERFORMING COMPUTATIONS IN CT SPACE")
	ctMean := ckks1.Mean(ct)
	fmt.Println(ctMean)

	//********************************
	fmt.Println("DECRYPTING - DECODING")
	ptMean := ckks1.Decrypt(ctMean, sk)

	//********************************
	msgMean := enc.Decode(ptMean)
	fmt.Println("Mean computed on CTs :")
	msgMean.PPrint()

	//********************************
	fmt.Println("MAX NORM OF ERRORS")
	msgMean = *msgMean.Scale(-1)
	msgMean.Add(&mean, &msgMean)
	//err := real(msgMean.SquaredNorm())
	err := msgMean.MaxNorm()
	fmt.Printf("%f \n", err)

}
