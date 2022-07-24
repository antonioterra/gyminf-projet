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

/* ----------------------------------------------------------*/
/* --------------------- USER PARAMETERS --------------------*/
/* ----------------------------------------------------------*/

/* ----------------- testing unit params --------------------*/

var tolerance = 0.1              // allowed diff betweed results and expected results
var boundForVectorEntries = 50.0 // bound on data values
var boundForGradeEntries = 6.0   // like boundForVectorEntries, but for vectors representing grades

/* ----------------- ckks related params, set by user --------------------*/

var NN = 64 // used below for the polynomial modulus x^NN + 1

var q0_nb_bits = 1000  // used to get q_L = q_0 * Delta^L
var delta_nb_bits = 20 // precision and modulus related param (greater means more precision and bigger q_L modulus)
var nb_levels = 10     // number of levels

var h = NN / 2 // sk distribution related param (must be < N)
var s2 = 2.0   // random errors related param

/* ----------------------------------------------------------*/
/* ------------------- COMPUTED PARAMETERS ------------------*/
/* ----------------------------------------------------------*/

var Q0 = q0(q0_nb_bits)
var QL = q0(q0_nb_bits + nb_levels*delta_nb_bits)

var deltaString = q0(delta_nb_bits)
var deltaBigInt, _ = (new(big.Int)).SetString(deltaString, 2)
var deltaBigFloat = (new(big.Float)).SetInt(deltaBigInt)
var deltaFloat64, _ = deltaBigFloat.Float64()

var baseScale = complex(deltaFloat64, 0)

/* ----------------------------------------------------------*/
/* ----------------------------------------------------------*/
/* ----------------------------------------------------------*/

// used for generation of moduli
func q0(n int) string {
	q0 := "1"
	for i := 0; i < n; i++ {
		q0 += "0"
	}
	return q0
}

// returns a random complex128 with real et imaginary parts in (-b, b)
func randComplex(b float64) complex128 {
	rand.Seed(time.Now().UnixNano())
	re := math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64()
	im := math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64()
	cpx := complex(re, im)
	return cpx
}

// returns a random complex128 z such that |z| < b
func randComplexBoundedNorm(b float64) complex128 {
	rand.Seed(time.Now().UnixNano())
	arg := 2 * math.Pi * rand.Float64()
	//fmt.Println("arg and COS :", arg, math.Cos(arg))
	re := b * math.Cos(arg)
	im := b * math.Sin(arg)
	cpx := complex(re, im)
	return cpx
}

// returns a random complex128 with real et imaginary parts in (-b, b)
func randIntComplex(b float64) complex128 {
	rand.Seed(time.Now().UnixNano())
	re := math.Round(math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64())
	im := math.Round(math.Pow(-1, float64(rand.Intn(2))) * b * rand.Float64())
	cpx := complex(re, im)
	return cpx
}

// returns a vector of n complex128 with parts in (-b, b) as []complex128
func randComplexVect(n int, b float64) []complex128 {

	rand.Seed(time.Now().UnixNano())
	vect := make([]complex128, n)
	for i, _ := range vect {
		vect[i] = randComplex(b)
	}
	return vect
}

// returns a vector of n complex128 z such that |z| < b as []complex128
func randComplexVectBoundedNorm(n int, b float64) []complex128 {

	rand.Seed(time.Now().UnixNano())
	vect := make([]complex128, n)
	for i, _ := range vect {
		vect[i] = randComplexBoundedNorm(b)
	}
	return vect
}

// returns a vector of n random course grades (as complex numbers)
func randGradesVect(n int) []complex128 {

	rand.Seed(time.Now().UnixNano())
	vect := make([]complex128, n)
	for i, _ := range vect {
		vect[i] = complex(float64(rand.Intn(int(boundForGradeEntries)*2-1))/2+1, 0)
	}
	return vect
}

// returns the max norm of <a> - <b>
func compare(a, b cMat.CMat) float64 {
	b = *b.Scale(-1)
	a.Add(&a, &b)
	return a.MaxNorm()
}

// tests the random.DG function
func _TestDG(t *testing.T) {
	N := 10000
	max := new(big.Int)
	max.SetString(QL, 2)
	fmt.Println(random.DG(N, s2, max))
}

// tests the copy of a poly.Poly returning a copy of a polynomial
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

// tests the poly.Deg function returning the degree of a polynomial
func _TestDeg(t *testing.T) {
	fmt.Println("Testing degree")

	coefsu := []*big.Int{big.NewInt(5), big.NewInt(10), big.NewInt(-6), big.NewInt(0)}
	u := poly.NewPoly(coefsu)

	fmt.Println("deg(u):", u.Degree())
}

// tests the poly.Mult function returning the product of two polynomials
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

// tests the poly.PolyDiv function returning the quotient and residue of a polynomial division
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

// tests the poly.Scale method scaling a pol in place
func _TestScaling(t *testing.T) {
	fmt.Println("Testing the scaling")
	coefsu := []*big.Int{big.NewInt(1), big.NewInt(6), big.NewInt(8), big.NewInt(3)}
	u := poly.NewPoly(coefsu)
	fmt.Println(u)
	u.Scale(big.NewInt(-20))
	fmt.Println(u)
}

//tests the poly.ScaleDiv function returning the quotient and residue of a polynomial division with roundings
func _TestScaleDiv(t *testing.T) {
	fmt.Println("Testing the division of a pol")
	coefsu := []*big.Int{big.NewInt(14914851435130), big.NewInt(643252345200), big.NewInt(8234534250), big.NewInt(12344354325643630)}
	u := poly.NewPoly(coefsu)
	fmt.Println(u)
	f := big.NewInt(10)
	u = poly.ScaleDiv(u, f)
	fmt.Println(u)
}

// tests the encoder.Encode function encoding a complex vector into a polynomial
func _TestEncoding(t *testing.T) {

	enc := encoder.NewEncoder(NN, baseScale)
	va := cMat.NewCMat(NN/2, 1, randComplexVect(NN/2, boundForVectorEntries))
	va.PPrint()

	ya := enc.Encode(&va)
	fmt.Println("Polynomial encoding of the vector a:")
	fmt.Println(ya)

	xa := enc.Decode(ya)
	fmt.Println("Decoded encoded vector :")
	xa.PPrint()

}

// used for message -- plaintext magnitude changes observations
func _TestEncodingBis(t *testing.T) {

	NN = 4
	bound := 100000.0
	baseScale = 1

	enc := encoder.NewEncoder(NN, baseScale)

	maxErr := 0.0

	for i := 0; i <= 100; i++ {
		data := randComplexVectBoundedNorm(NN/2, bound)
		va := cMat.NewCMat(NN/2, 1, data)
		ya := enc.Encode(&va)
		fmt.Println(ya)
		za := enc.Decode(ya)
		diff := compare(va, za)
		if diff > maxErr {
			maxErr = diff
		}

	}
	fmt.Println("----------------------------")
	fmt.Println(maxErr)
	fmt.Println("----------------------------")

}

// tests the ckks.Encrypt function encrypting a plaintext
func _TestEncrypt(t *testing.T) {
	fmt.Println("Testing the encryption")

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, NN, h, nb_levels, s2)

	va := cMat.NewCMat(NN/2, 1, randComplexVect(NN/2, boundForVectorEntries))

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

// tests the correctness through the ckks.Encrypt - ckks.Decrypt cycle
func TestEncryptDecrypt(t *testing.T) {
	fmt.Println("TESTING THE ENCRYPTION - DECRYPTION")

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, NN, h, nb_levels, s2)

	va := cMat.NewCMat(NN/2, 1, randComplexVect(NN/2, boundForVectorEntries))

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

// tests the ckks.CTAdd function returning the sum of two ciphertexts
func TestHomomAdd(t *testing.T) {
	for i := 0; i < 1; i = i + 1 {
		fmt.Println("TESTING HOMOMORPHISM ON +")

		Q := new(big.Int)
		Q.SetString(QL, 2)
		P := Q
		ckks := ckks.NewCKKS(Q, P, NN, h, nb_levels, s2)

		//********************************
		v1 := cMat.NewCMat(NN/2, 1, randComplexVect(NN/2, boundForVectorEntries))
		v2 := cMat.NewCMat(NN/2, 1, randComplexVect(NN/2, boundForVectorEntries))
		vs := cMat.NewCMat(NN/2, 1, make([]complex128, NN/2))

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

// tests the ckks.SKeyGen, ckks.PKeyGen and ckks.EvKeyGen functions returning sk, pk and evk
func _TestKeyGen(t *testing.T) {
	fmt.Println("Testing the key generation")

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, NN, h, nb_levels, s2)

	//********************************
	fmt.Println("GENERATING KEYS")
	sk := ckks.SKeyGen()
	fmt.Println("sk :", sk)

	pk := ckks.PKeyGen(sk)
	fmt.Println("pk", pk)

	evk := ckks.EvKeyGen(sk)
	fmt.Println("evk", evk)

}

// tests the ckks.CTMult function returning the product of two ciphertexts

func TestHomomMult(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	fmt.Println("TESTING HOMOMORPHISM ON *")

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, NN, h, nb_levels, s2)
	enc := encoder.NewEncoder(ckks.N, baseScale)

	//********************************
	v1 := cMat.NewCMat(NN/2, 1, randComplexVect(NN/2, float64(boundForVectorEntries)))
	v2 := cMat.NewCMat(NN/2, 1, randComplexVect(NN/2, float64(boundForVectorEntries)))
	vp := cMat.NewCMat(NN/2, 1, make([]complex128, NN/2))

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

// tests the encoder.ConstToPt function returning the plaintext corresponding to a vector having all entries set to the same constant
func _TestConstant(t *testing.T) {
	k := 4.0

	enc := encoder.NewEncoder(NN, baseScale)

	pt := enc.ConstToPT(k)
	fmt.Println("Encoding of constant vector (k=4) with baseScale = 64")
	fmt.Println(pt)
}

// tests the ckks.CTScale method changing the scale of a ciphertext without modifying its scale
func TestHomScaling(t *testing.T) {
	fmt.Println("TESTING HOMOM SCALING")

	k := -10

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q
	ckks := ckks.NewCKKS(Q, P, NN, h, nb_levels, s2)
	enc := encoder.NewEncoder(ckks.N, baseScale)

	//********************************
	v := cMat.NewCMat(NN/2, 1, randComplexVect(NN/2, boundForVectorEntries))

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

// tests the ckks.CTIncScale method encreasing the scale of a ciphertext and modifying its scale
func TestCTIncScale(t *testing.T) {
	fmt.Println("TESTING CT INCREASE SCALE")

	k := new(big.Int)
	k.SetString("10", 2)

	Q := new(big.Int)
	Q.SetString(QL, 2)
	P := Q

	ckks := ckks.NewCKKS(Q, P, NN, h, nb_levels, s2)
	enc := encoder.NewEncoder(ckks.N, baseScale)

	//********************************
	v := cMat.NewCMat(NN/2, 1, randComplexVect(NN/2, boundForVectorEntries))

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

// tests the ckks.Mean method computing the mean of a list of ciphertexts
func TestMean(t *testing.T) {
	fmt.Println("TESTING MEAN")

	nbStudents := NN
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
	ctMean := ckks1.Mean(ct, pk, evk) //, deltaBigInt)
	ckks1.RS(&ctMean, deltaBigInt)
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

// tests the ckks.RS method rescaling a ciphertext
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

	ckks := ckks.NewCKKS(Q, P, NN, h, nb_levels, s2)
	enc := encoder.NewEncoder(ckks.N, baseScale)

	va := cMat.NewCMat(NN/2, 1, randComplexVect(NN/2, boundForVectorEntries))

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

// tests the ckks.Var method computing the variance of a list of ciphertexts
func TestVar(t *testing.T) {

	fmt.Println("TESTING VAR")
	rand.Seed(time.Now().UnixNano())

	nbStudents := NN
	nbCourses := 50

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
	ctVar := ckks1.Var(cts, pk, evk, deltaBigInt)

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
