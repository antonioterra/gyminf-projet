package main

import (
	"fmt"
	"math/big"
	"math/rand"
	"time"

	"kazat.ch/lbcrypto/cMat"
	"kazat.ch/lbcrypto/ckks"
	"kazat.ch/lbcrypto/encoder"
)

// Returns a vector of n random course grades (as complex numbers)
func randGradesVect6(n int) []complex128 {

	rand.Seed(time.Now().UnixNano())
	vect := make([]complex128, n)
	for i, _ := range vect {
		vect[i] = complex(float64(rand.Intn(int(6)*2-1))/2+1, 0)
	}
	return vect
}

// Returns a binary string representation of 2^n
func binStrRep(n int) string {
	q0 := "1"
	for i := 0; i < n; i++ {
		q0 += "0"
	}
	return q0
}

// Returns the scaling factor as a *big.Int
func GetDelta(delta_nb_bits int) *big.Int {

	deltaString := binStrRep(delta_nb_bits)
	deltaBigInt, _ := (new(big.Int)).SetString(deltaString, 2)

	return deltaBigInt
}

// Returns an encoder for data vectors of length <vectLen>, with a base scaling factor given by <delta_nb_bits>
func GetEncoder(vectLen, delta_nb_bits int) encoder.Encoder {

	deltaBigInt := GetDelta(delta_nb_bits)
	deltaBigFloat := (new(big.Float)).SetInt(deltaBigInt)
	deltaFloat64, _ := deltaBigFloat.Float64()

	baseScale := complex(deltaFloat64, 0)

	return encoder.NewEncoder(vectLen, baseScale)
}

// Returns a CKKS scheme with h = N/2 and s2 = 3.2
// Should be defined otherwise to give security control to users
func GetCKKS(N, delta_nb_bits, q0_nb_bits, nb_levels int) ckks.CKKS {

	Q := binStrRep(q0_nb_bits + nb_levels*delta_nb_bits)
	QL := new(big.Int)
	QL.SetString(Q, 2)

	P := QL

	h := N / 2
	s2 := 3.2

	return ckks.NewCKKS(QL, P, N, h, nb_levels, s2)
}

func main() {

	// Situation defined parameters
	nbStudents := 16 // includes fake students when not a power of 2
	nbCourses := 10

	// Generating data based on the above parameters
	// vectors[i] will contains the grades for the course #i
	vectors := make([]cMat.CMat, nbCourses)
	for i := 0; i < nbCourses; i++ {
		vectors[i] = cMat.NewCMat(nbStudents, 1, randGradesVect6(nbStudents))
	}

	// Setting parameters for CKKS on client side
	N := 2 * nbStudents // fake students can be added for security and getting a power of 2
	delta_nb_bits := 20
	q0_nb_bits := 100
	nb_levels := 10

	// Generating encoder and CKKS instance based on the above parameters on client side
	delta := GetDelta(delta_nb_bits) // used for the RS procedure in Var
	enc := GetEncoder(N, delta_nb_bits)
	ckks1 := GetCKKS(N, delta_nb_bits, q0_nb_bits, nb_levels)

	// Generating keys on client side
	sk := ckks1.SKeyGen()
	pk := ckks1.PKeyGen(sk)
	evk := ckks1.EvKeyGen(sk)

	// Encoding
	// pts[i] will contain the plaintext corresponding to the grades of course #i
	pts := make([]encoder.PT, nbCourses)
	for i := 0; i < nbCourses; i++ {
		pts[i] = enc.Encode(&vectors[i])
	}

	// Encryption
	// cts[i] will contain the plaintext corresponding to the grades of course #i
	cts := make([]ckks.CT, nbCourses)
	for i := 0; i < nbCourses; i++ {
		cts[i] = ckks1.Encrypt(pts[i], pk)
	}

	// Computing mean and var on ciphertexts on server side
	meanCT := ckks1.Mean(cts, pk, evk)
	varCT := ckks1.Var(cts, pk, evk, delta)

	// Decryption and decoding on client side
	meanPT := ckks1.Decrypt(meanCT, sk)
	varPT := ckks1.Decrypt(varCT, sk)
	meanClear := enc.Decode(meanPT)
	varClear := enc.Decode(varPT)

	// Printing grades, mean and var for each studend
	for i := 0; i < nbStudents; i++ {
		for j := 0; j < nbCourses; j++ {
			fmt.Print(real(vectors[j].GetData()[i]), " ")
		}
		fmt.Print(real(meanClear.GetData()[i]), " ", real(varClear.GetData()[i]), "\n")
	}
}
