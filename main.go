package main

import (
	"fmt"
	"math/rand"

	"kazat.ch/lbcrypto/polyMod"
)

func Equal(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i, v := range a {
		if v != b[i] {
			return false
		}
	}
	return true
}

func main() {

	mod := 101
	nbtests := 0
	for i := 0; i < 10000; i++ {
		borne := 100
		u := polyMod.NewPolyInt(rand.Perm(mod)[:rand.Intn(borne)], mod)
		v := polyMod.NewPolyInt(rand.Perm(mod)[:rand.Intn(borne)], mod)
		p := polyMod.PolyIntMultMod(&u, &v)
		r := polyMod.NewPolyInt(rand.Perm(mod)[:v.Degree()], mod)

		if p.Degree() > r.Degree() && v.Degree() > r.Degree() {

			nbtests++
			p.PolyAddMod(&r)
			q, re, erCode := polyMod.PolyIntDivMod(&p, &v)

			//fmt.Println(i, erCode, p.Degree(), v.Degree(), r.Degree())
			if (!Equal(u.Coefs, q.Coefs) || !Equal(r.Coefs, re.Coefs)) && erCode == nil {
				fmt.Println("-----------------------")
				fmt.Println("u :", u, "deg", u.Degree())
				fmt.Println("v :", v, "deg", v.Degree())
				fmt.Println("r :", r, "deg", r.Degree())
				fmt.Println("u*v + r =", p, "deg", p.Degree())
				fmt.Println("q :", q)
				fmt.Println("rest :", re)
				fmt.Println("-----------------------")
				//polyMod.PolyIntDivModTestVerb(&p, &v)
				break
			}
		}

	}
	fmt.Println(nbtests, "tests done")

	/* var z = 1 + 1i
	cmplx.Conj(z)
	fmt.Println(cmplx.Conj(z))

	var pt = []complex128{1 + 1i, -1 + 1i, -1 - 1i, 1 - 1i} */

}
