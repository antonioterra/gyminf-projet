package main

import (
	"fmt"

	"kazat.ch/lbcrypto/cMat"
	"kazat.ch/lbcrypto/encoder"
)

func main() {

	// MODIFIER POUR AVOIR DU RANDOM ROUNDING

	encoder := encoder.NewEncoer(8.0, 64, 727)
	va := cMat.NewCMat(2, 1, []complex128{3 + 4i, 2 - 1i})
	fmt.Println("Original vector :")
	va.PPrint()
	x := encoder.Encode(&va)
	fmt.Println("Polynomial encoding of the vector :", x)
	y := encoder.Decode(x)
	fmt.Println("Decoded encoded vector :")
	y.PPrint()
}
