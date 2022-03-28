package main

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"testing"
	"time"

	"kazat.ch/lbcrypto/cMat"
	"kazat.ch/lbcrypto/ckks"
	"kazat.ch/lbcrypto/encoder"
	"kazat.ch/lbcrypto/poly"
	"kazat.ch/lbcrypto/random"
)

// ATTENTION : peut-être UN PROBLEMENT DANS LA CONVERSION DES DELTA ENTRE BIG.INT -> FLOAT64
var tolerance = 0.1
var N = 128
var s2 = 4.0 // variance of the DG distribution

var q0_nb_bits = 1000 //1000
var delta_nb_bits = 20
var nb_levels = 10

var Q0 = q0(q0_nb_bits)
var QL = q0(q0_nb_bits + nb_levels*delta_nb_bits)

var deltaString = q0(delta_nb_bits)
var deltaBigInt, _ = (new(big.Int)).SetString(deltaString, 2)
var deltaBigFloat = (new(big.Float)).SetInt(deltaBigInt)
var deltaFloat64, _ = deltaBigFloat.Float64()

var baseScale = complex(deltaFloat64, 0)

var boundForVectorEntries = 50.0
var boundForGradeEntries = 50.0

var h = N / 2 // must be < N

func q0(n int) string {
	q0 := "1"
	for i := 0; i < n; i++ {
		q0 += "0"
	}
	return q0
}

func _TestParams(t *testing.T) {
	fmt.Println("bscale = ", baseScale)
	Q := new(big.Int)
	Q.SetString(QL, 2)
	fmt.Println("Q = ", Q)
	fmt.Println("QL = ", QL)
	fmt.Println("DeltaBigFloat and DeltaFloat64 : ", deltaBigFloat, deltaFloat64)

}

// return a random complex128 with real et imaginary parts in (-b, b)
func randComplex(b float64) complex128 {
	rand.Seed(time.Now().UnixNano())
	re := math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64()
	im := math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64()
	cpx := complex(re, im)
	return cpx
}

// return a random complex128 with real et imaginary parts in (-b, b)
func randIntComplex(b float64) complex128 {
	rand.Seed(time.Now().UnixNano())
	re := math.Round(math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64())
	im := math.Round(math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64())
	cpx := complex(re, im)
	return cpx
}

// return a vector of n complex128 with parts in (-b, b) as []complex128
func randComplexVect(n int, b float64) []complex128 {

	rand.Seed(time.Now().UnixNano())
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
		vect[i] = complex(float64(rand.Intn(int(boundForGradeEntries)*2-1))/2+1, 0)
	}
	return vect
}

func compare(a, b cMat.CMat) float64 {
	b = *b.Scale(-1)
	a.Add(&a, &b)
	return a.MaxNorm()
}

func _TestDG(t *testing.T) {
	N := 10000
	max := new(big.Int)
	max.SetString(QL, 2)
	fmt.Println(random.DG(N, s2, max))
}

func _TestPolyCopy(t *testing.T) {
	fmt.Println("Testing the polynomial division")

	coefsu := []*big.Int{big.NewInt(5), big.NewInt(10), big.NewInt(2), big.NewInt(3)}
	u := poly.NewPoly(coefsu)
	v := poly.Copy(u)

	fmt.Println("u and v: ", u, v)
	u.Coefs[0].Set(big.NewInt(20))
	fmt.Println("u and v: ", u, v)
	v.Coefs[0].Set(big.NewInt(50))
	fmt.Println("u and v: ", u, v)
}

func _TestDeg(t *testing.T) {
	fmt.Println("Testing degree")

	coefsu := []*big.Int{big.NewInt(5), big.NewInt(10), big.NewInt(-6), big.NewInt(0)}
	u := poly.NewPoly(coefsu)

	fmt.Println("deg(u):", u.Degree())
}

func _TestPolyMult(t *testing.T) {
	fmt.Println("Testing the polynomial multiplication")

	coefsu := []*big.Int{big.NewInt(0), big.NewInt(-1), big.NewInt(-1), big.NewInt(1)}
	coefsv := []*big.Int{big.NewInt(-1), big.NewInt(-6), big.NewInt(6), big.NewInt(-4)}

	u := poly.NewPoly(coefsu)
	v := poly.NewPoly(coefsv)

	fmt.Println("u and v befor mul: ", u, v)
	p := poly.Mult(u, v)

	fmt.Println("u and v after mult: ", u, v)

	fmt.Println("prod: ", p)
}

func _TestPolyDiv(t *testing.T) {
	fmt.Println("Testing the polynomial division")

	coefsu := []*big.Int{big.NewInt(5), big.NewInt(10), big.NewInt(2), big.NewInt(3), big.NewInt(0), big.NewInt(-4)}
	coefsv := []*big.Int{big.NewInt(1), big.NewInt(0), big.NewInt(-1)}

	u := poly.NewPoly(coefsu)
	v := poly.NewPoly(coefsv)

	fmt.Println("u and v befor division: ", u, v)
	q, r := poly.PolyDiv(u, v)
	fmt.Println("u and v after division: ", u, v)

	fmt.Println("q and r: ", q, r)
}

// tests the scaling of a pol in place
func _TestScaling(t *testing.T) {
	fmt.Println("Testing the scaling")
	coefsu := []*big.Int{big.NewInt(1), big.NewInt(6), big.NewInt(8), big.NewInt(3)}
	u := poly.NewPoly(coefsu)
	fmt.Println(u)
	u.Scale(big.NewInt(-20))
	fmt.Println(u)
}

//tests the scaling of a pol in place

func _TestScaleDiv(t *testing.T) {
	fmt.Println("Testing the division of a pol")
	coefsu := []*big.Int{big.NewInt(14914851435130), big.NewInt(643252345200), big.NewInt(8234534250), big.NewInt(12344354325643630)}
	u := poly.NewPoly(coefsu)
	fmt.Println(u)
	f := big.NewInt(10)
	u = poly.ScaleDiv(u, f)
	fmt.Println(u)
}

// tests the encoder
func _TestEncoding(t *testing.T) {

	enc := encoder.NewEncoder(N, baseScale)
	va := cMat.NewCMat(2, 1, randComplexVect(2, boundForVectorEntries))
	va.PPrint()

	ya := enc.Encode(&va)
	fmt.Println("Polynomial encoding of the vector a:")
	fmt.Println(ya)

	xa := enc.Decode(ya)
	fmt.Println("Decoded encoded vector :")
	xa.PPrint()

}

// tests the encryption
func _TestEncrypt(t *testing.T) {
	fmt.Println("Testing the encryption")

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)

	va := cMat.NewCMat(N/2, 1, randComplexVect(N/2, boundForVectorEntries))

	//********************************
	fmt.Println("ENCODING")
	enc := encoder.NewEncoder(ckks.N, baseScale)
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

// tests the encryption - decryption
func TestEncryptDecrypt(t *testing.T) {
	fmt.Println("TESTING THE ENCRYPTION - DECRYPTION")

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)

	va := cMat.NewCMat(N/2, 1, randComplexVect(N/2, boundForVectorEntries))

	//fmt.Println("message :")
	//va.PPrint()

	//********************************
	//fmt.Println("ENCODING :")
	enc := encoder.NewEncoder(ckks.N, baseScale)
	pta := enc.Encode(&va)

	//********************************
	//fmt.Println("KEYS GENERATION :")
	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)

	//********************************
	//fmt.Println("ENCRYPTING :")
	ct := ckks.Encrypt(pta, pk)
	//fmt.Println("ct :", ct.B)

	//********************************
	//fmt.Println("DECRYPTING :")
	pt := ckks.Decrypt(ct, sk)
	//fmt.Println("pt :", pt)

	//********************************
	//fmt.Println("DECODING")
	vect := enc.Decode(pt)
	//vect.PPrint()

	//********************************
	err := compare(vect, va)
	fmt.Printf("max norm of errors : %f \n", err)
	if err > tolerance {
		t.Fail()
	}
}

func TestHomomAdd(t *testing.T) {
	for i := 0; i < 1; i = i + 1 {
		fmt.Println("TESTING HOMOMORPHISM ON +")

		Q := new(big.Int)
		Q.SetString(QL, 2)
		P := Q
		ckks := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)

		//********************************
		v1 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, boundForVectorEntries))
		v2 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, boundForVectorEntries))
		vs := cMat.NewCMat(N/2, 1, make([]complex128, N/2))

		vs.Add(&v1, &v2)

		//fmt.Println("sum of messages :")
		//vs.PPrint()

		//********************************
		//fmt.Println("ENCODING")
		enc := encoder.NewEncoder(ckks.N, baseScale)
		pt1 := enc.Encode(&v1)
		pt2 := enc.Encode(&v2)
		//fmt.Println("pts :", enc.Encode(&vs))

		//********************************
		//fmt.Println("GENERAING KEYS")
		sk := ckks.SKeyGen()
		//fmt.Println("sk :", sk)
		pk := ckks.PKeyGen(sk)
		//fmt.Println("pk :", pk)

		//********************************
		//fmt.Println("ENCRYPTING")
		ct1 := ckks.Encrypt(pt1, pk)
		ct2 := ckks.Encrypt(pt2, pk)

		//fmt.Println("ct1 :", ct1)
		//fmt.Println("ct2 :", ct2)

		cts := ckks.CTAdd(ct1, ct2)
		//fmt.Println("cts :", cts)

		//********************************
		//fmt.Println("DECRYPTING :")
		pts := ckks.Decrypt(cts, sk)
		//fmt.Println("pts :", pts)

		//********************************
		//fmt.Println("DECODING")
		vect := enc.Decode(pts)
		//vect.PPrint()

		//********************************
		err := compare(vect, vs)
		fmt.Printf("max norm of errors : %f \n", err)
		if err > tolerance {
			t.Fail()
		}
	}
}

// tests the keygeneration
func _TestKeyGen(t *testing.T) {
	fmt.Println("Testing the key generation")

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)

	//********************************
	fmt.Println("GENERATING KEYS")
	sk := ckks.SKeyGen()
	fmt.Println("sk :", sk)

	pk := ckks.PKeyGen(sk)
	fmt.Println("pk", pk)

	evk := ckks.EvKeyGen(sk)
	fmt.Println("evk", evk)

}

func TestHomomMult(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	fmt.Println("TESTING HOMOMORPHISM ON *")

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)
	enc := encoder.NewEncoder(ckks.N, baseScale)

	//********************************
	v1 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, float64(boundForVectorEntries)))
	v2 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, float64(boundForVectorEntries)))
	vp := cMat.NewCMat(N/2, 1, make([]complex128, N/2))

	vp.CoefWiseProd(&v1, &v2)

	//fmt.Println("should get:")
	//vp.PPrint()

	//********************************
	//fmt.Println("ENCODING")
	pt1 := enc.Encode(&v1)
	pt2 := enc.Encode(&v2)
	//fmt.Println("pt1 :", pt1)
	//fmt.Println("pt2 :", pt2)

	//fmt.Println("pt of prod :", enc.Encode(&vp))

	//********************************
	//fmt.Println("GENERAING KEYS")
	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)
	evk := ckks.EvKeyGen(sk)

	//********************************
	//fmt.Println("ENCRYPTING")
	ct1 := ckks.Encrypt(pt1, pk)
	ct2 := ckks.Encrypt(pt2, pk)

	//********************************
	//fmt.Println("PRODUCT OF CYPHERTEXTS")
	ctp := ckks.CTMult(ct1, ct2, evk)
	//fmt.Println("ctp :", ctp)
	//fmt.Println("ctp scale :", ctp.Scale)
	//********************************
	//fmt.Println("DECRYPTING PRODUCT:")
	ptp := ckks.Decrypt(ctp, sk)

	vect := enc.Decode(ptp)
	//vect.PPrint()

	//********************************
	err := compare(vect, vp)
	fmt.Printf("max norm of errors : %f \n", err)
	if err > tolerance {
		t.Fail()
	}

}

func _TestConstant(t *testing.T) {
	k := 4.0

	enc := encoder.NewEncoder(N, baseScale)

	pt := enc.ConstToPT(k)
	fmt.Println("Encoding of constant vector (k=4) with baseScale = 64")
	fmt.Println(pt)
}

func TestHomScalingBis(t *testing.T) {
	fmt.Println("TEST HOMOM SCALING")

	k := -107

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)
	enc := encoder.NewEncoder(ckks.N, baseScale)

	//********************************
	v := cMat.NewCMat(N/2, 1, randComplexVect(N/2, boundForVectorEntries))

	//v.PPrint()

	//fmt.Println("---------")
	vscaled := v.Copy()
	vscaled.Scale(complex(float64(k), 0))
	//vscaled.PPrint()

	//********************************
	//fmt.Println("GENERAING KEYS")
	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)

	//********************************
	//fmt.Println("ENCODING")
	pt := enc.Encode(&v)

	//********************************
	//fmt.Println("ENCRYPTING")
	ct := ckks.Encrypt(pt, pk)

	//********************************
	//fmt.Println("PERFORMING COMPUTATIONS")

	//fmt.Println("ct before :", ct)
	ct.CTScale(big.NewInt(int64(k)))
	//fmt.Println("ct after :", ct)

	//********************************
	//fmt.Println("DECRYPTING :")
	pt2 := ckks.Decrypt(ct, sk)

	//********************************
	//fmt.Println("DECODING")
	vect := enc.Decode(pt2)
	//vect.PPrint()

	//********************************
	err := compare(vect, vscaled)
	fmt.Printf("max norm of errors : %f \n", err)
	if err > tolerance {
		t.Fail()
	}

}

func TestCTIncScale(t *testing.T) {
	fmt.Println("TEST CT INCREASE SCALE")

	k := new(big.Int)
	k.SetString("10", 2)

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q

	ckks := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)
	enc := encoder.NewEncoder(ckks.N, baseScale)

	//********************************
	v := cMat.NewCMat(N/2, 1, randComplexVect(N/2, boundForVectorEntries))

	//v.PPrint()

	//********************************
	//fmt.Println("GENERAING KEYS")
	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)

	//********************************
	//fmt.Println("ENCODING")
	pt := enc.Encode(&v)

	//********************************
	//fmt.Println("ENCRYPTING")
	ct := ckks.Encrypt(pt, pk)

	//********************************
	//fmt.Println("PERFORMING COMPUTATIONS")

	//fmt.Println("ct before :", ct)
	ct.CTIncScale(k)
	//fmt.Println("ct after :", ct)

	//********************************
	//fmt.Println("DECRYPTING :")
	pt2 := ckks.Decrypt(ct, sk)

	//********************************
	//fmt.Println("DECODING")
	vect := enc.Decode(pt2)
	//vect.PPrint()

	//********************************
	err := compare(vect, v)
	fmt.Printf("max norm of errors : %f \n", err)
	if err > tolerance {
		t.Fail()
	}

}
func TestMean(t *testing.T) {
	fmt.Println("TEST MEAN")

	nbStudents := 50
	nbCourses := 10

	N := 2 * nbStudents
	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q

	ckks1 := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)
	enc := encoder.NewEncoder(ckks1.N, baseScale)

	//********************************

	//fmt.Println("Number of students :", nbStudents)
	//fmt.Println("Number of courses :", nbCourses)

	vectors := make([]cMat.CMat, nbCourses)
	mean := cMat.NewCMat(N/2, 1, randComplexVect(N/2, 0))
	for i := 0; i < nbCourses; i++ {
		vectors[i] = cMat.NewCMat(nbStudents, 1, randGradesVect(nbStudents))
		mean.Add(&mean, &vectors[i])
		//fmt.Println("Means for course n° ", i)
		//vectors[i].PPrint()
		//fmt.Println("")
	}

	mean.Scale(1 / complex(float64(nbCourses), 0))
	//fmt.Println("Mean :")
	//mean.PPrint()

	//********************************
	//fmt.Println("ENCODING-ENCODING DATA")

	pt := make([]encoder.PT, nbCourses)
	for i := 0; i < nbCourses; i++ {
		pt[i] = enc.Encode(&vectors[i])
	}

	//********************************
	sk := ckks1.SKeyGen()
	pk := ckks1.PKeyGen(sk)
	evk := ckks1.EvKeyGen(sk)
	//********************************
	ct := make([]ckks.CT, nbCourses)

	for i := 0; i < nbCourses; i++ {
		ct[i] = ckks1.Encrypt(pt[i], pk)
	}
	//********************************
	//fmt.Println("PERFORMING COMPUTATIONS IN CT SPACE")
	ctMean := ckks1.Mean(ct, pk, evk, deltaBigInt)
	//fmt.Println("ctmean scale :", ctMean.Scale)
	//********************************
	//fmt.Println("DECRYPTING - DECODING")
	ptMean := ckks1.Decrypt(ctMean, sk)

	//********************************
	msgMean := enc.Decode(ptMean)
	//fmt.Println("Mean computed on CTs after decryption:")
	//msgMean.PPrint()

	//********************************
	err := compare(mean, msgMean)
	fmt.Printf("max norm of errors : %f \n", err)
	if err > tolerance {
		t.Fail()
	}
}

func TestRS(t *testing.T) {
	fmt.Println("TESTING RS")

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q

	//floatScale := math.Pow(2, 10)
	//delta := complex(floatScale, 0)
	//baseScale := delta * delta
	//bigBaseScale := big.NewInt(int64(floatScale))
	//Q.Mul(Q, bigBaseScale)
	//Q.Mul(Q, bigBaseScale)

	ckks := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)
	enc := encoder.NewEncoder(ckks.N, baseScale)

	va := cMat.NewCMat(N/2, 1, randComplexVect(N/2, boundForVectorEntries))

	//va.PPrint()

	//********************************
	//fmt.Println("ENCODING")
	pt := enc.Encode(&va)
	//fmt.Println("pt: ", pt)

	//********************************
	//fmt.Println("GENERATING KEYS")

	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)

	//********************************
	//fmt.Println("ENCRYPTING")
	ct := ckks.Encrypt(pt, pk)
	//fmt.Println("ct :", ct)

	ct.CTIncScale(deltaBigInt)
	ct.CTIncScale(deltaBigInt)
	ct.CTIncScale(deltaBigInt)

	ckks.RS(&ct, deltaBigInt)

	//fmt.Println("ct :", ct)

	//fmt.Println("DECRYPTING")

	pt2 := ckks.Decrypt(ct, sk)

	//fmt.Println("pt2: ", pt2)
	res := enc.Decode(pt2)
	//res.PPrint()

	//********************************
	err := compare(res, va)
	fmt.Printf("max norm of errors : %f \n", err)
	if err > tolerance {
		t.Fail()
	}
}

func TestVar(t *testing.T) {

	fmt.Println("TEST VAR")
	rand.Seed(time.Now().UnixNano())

	nbStudents := 50
	nbCourses := 10

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	N := 2 * nbStudents

	ckks1 := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)

	//********************************

	//fmt.Println("Number of students :", nbStudents)
	//fmt.Println("Number of courses :", nbCourses)

	//each element of vectors is a vector of grades (on per course)
	vectors := make([]cMat.CMat, nbCourses)
	mean := cMat.NewCMat(N/2, 1, randComplexVect(N/2, 0))

	for i := 0; i < nbCourses; i++ {
		vectors[i] = cMat.NewCMat(nbStudents, 1, randGradesVect(nbStudents))
		mean.Add(&mean, &vectors[i])
	}

	//********************************
	//fmt.Println("ENCODING-ENCODING DATA")

	enc := encoder.NewEncoder(ckks1.N, baseScale)
	pts := make([]encoder.PT, nbCourses)
	for i := 0; i < nbCourses; i++ {
		pts[i] = enc.Encode(&vectors[i])
	}

	//********************************
	sk := ckks1.SKeyGen()
	pk := ckks1.PKeyGen(sk)
	evk := ckks1.EvKeyGen(sk)
	//********************************
	cts := make([]ckks.CT, nbCourses)

	for i := 0; i < nbCourses; i++ {
		cts[i] = ckks1.Encrypt(pts[i], pk)
	}
	//********************************

	//fmt.Println("PERFORMING COMPUTATIONS IN CT SPACE")
	ctVar := ckks1.Var(cts, sk, pk, evk, deltaBigInt)

	//ckks1.RS(&ctVar, deltaBigInt)
	//********************************
	//fmt.Println("DECRYPTING - DECODING")
	ptVar := ckks1.Decrypt(ctVar, sk)
	//fmt.Println("computed pt: ", ptVar)
	//fmt.Println("pt as decryption of ct: ", ptVar)

	//********************************
	msgVar := enc.Decode(ptVar)
	//fmt.Println("Var computed on CTs after decryption:")
	//msgVar.PPrint()
	//fmt.Println("---")

	//********************************

	//calcul de la vraie variance
	mean.Scale(-1 / complex(float64(nbCourses), 0))

	sqrdiff := make([]cMat.CMat, nbCourses)
	for i := 0; i < nbCourses; i++ {
		sqrdiff[i] = cMat.NewCMat(nbStudents, 1, randGradesVect(nbStudents))
		sqrdiff[i].Add(&vectors[i], &mean)
		sqrdiff[i].CoefWiseProd(&sqrdiff[i], &sqrdiff[i])
	}

	sig2 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, 0))
	for i := 0; i < nbCourses; i++ {
		sig2.Add(&sig2, &sqrdiff[i])
	}
	sig2.Scale(1 / complex(float64(nbCourses), 0))

	//ptcheck := enc.Encode(&sig2)
	//fmt.Println("encoding of tru var: ", ptcheck)
	//fmt.Println("Var computed on original data :")
	//sig2.PPrint()
	//fmt.Println("---------------------------------")

	//---------------------------------
	err := compare(sig2, msgVar)
	fmt.Printf("max norm of errors : %f \n", err)
	if err > 0.1 {
		t.Fail()
	}
}

func TestRSBis(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	fmt.Println("TESTING RS AFTER *")

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, N, h, nb_levels, s2)
	enc := encoder.NewEncoder(ckks.N, baseScale)

	//********************************
	v1 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, float64(boundForVectorEntries)))
	v2 := cMat.NewCMat(N/2, 1, randComplexVect(N/2, float64(boundForVectorEntries)))
	vp := cMat.NewCMat(N/2, 1, make([]complex128, N/2))

	vp.CoefWiseProd(&v1, &v2)

	//fmt.Println("should get:")
	//vp.PPrint()

	//********************************
	//fmt.Println("ENCODING")
	pt1 := enc.Encode(&v1)
	pt2 := enc.Encode(&v2)
	//fmt.Println("pt1 :", pt1)
	//fmt.Println("pt2 :", pt2)

	//fmt.Println("pt of prod :", enc.Encode(&vp))

	//********************************
	//fmt.Println("GENERAING KEYS")
	sk := ckks.SKeyGen()
	pk := ckks.PKeyGen(sk)
	evk := ckks.EvKeyGen(sk)

	//********************************
	//fmt.Println("ENCRYPTING")
	ct1 := ckks.Encrypt(pt1, pk)
	ct2 := ckks.Encrypt(pt2, pk)

	//********************************
	//fmt.Println("PRODUCT OF CYPHERTEXTS")
	ctp := ckks.CTMult(ct1, ct2, evk)
	//fmt.Println("ctp befors RS: ", ctp.Mod)
	ckks.RS(&ctp, deltaBigInt)
	//fmt.Println("ctp after RS: ", ctp.Mod)

	//fmt.Println("ctp :", ctp)
	//fmt.Println("ctp scale :", ctp.Scale)
	//********************************
	//fmt.Println("DECRYPTING PRODUCT:")
	ptp := ckks.Decrypt(ctp, sk)
	vect := enc.Decode(ptp)
	//vect.PPrint()

	//********************************
	err := compare(vect, vp)
	fmt.Printf("max norm of errors : %f \n", err)
	if err > tolerance {
		t.Fail()
	}

}
