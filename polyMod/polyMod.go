package polyMod

import (
	"errors"
	"fmt"
)

// représente des polynômes à coefficients entiers
// les coefficients seront toujours pris dans Z_p avec p premier
// {1, 2, 3} correspond à 1 + 2x + 3x^2
// TOUJOURS CREER LES POLYNOMES AVEC LE CONSTRUCTEUR CI-DESSOUS POUR GARANTIR QUE LES COEFS SOIENT OK
type PolyInt struct {
	Coefs []int
	mod   int
}

// retourne le degré d'un polynôme
func (poly *PolyInt) Degree() int {
	var mod = poly.mod
	var deg = 0
	for i := range poly.Coefs {
		if poly.Coefs[len(poly.Coefs)-1-i]%mod != 0 {
			deg = len(poly.Coefs) - i - 1
			break
		}
	}
	return deg
}

// Enlève les 0*x^n dominants
// pas exportée : quand les poly sont créés avec le constructeur, il sont déjà deflatés
func (pol *PolyInt) deflate() {
	pol.Coefs = pol.Coefs[:pol.Degree()+1]
}

// retourne un nouveau PolyInt avec coefficients dans Z_mod
func NewPolyInt(coefs []int, mod int) PolyInt {
	for i := range coefs {
		coefs[i] = (coefs[i]%mod + mod) % mod
	}
	var poly = PolyInt{Coefs: coefs, mod: mod}
	poly.deflate()
	return poly
}

// modifie les coefficients de pol en en prenant le modulo
func (pol *PolyInt) Mod(mod int) *PolyInt {
	for i, c := range pol.Coefs {
		x := (c%mod + mod) % mod
		pol.Coefs[i] = x
	}
	return pol
}

// renvoie une copie modularisées de pol
func Mod(pol *PolyInt, mod int) PolyInt {
	var res = Copy(pol)
	res.Mod(pol.mod)
	return res
}

// returns x, y and gcd in ax + by = gcd(a, b)
// piqué sur le web
func extanded_gcd(a, b int) (int, int, int) {
	var oldr, r = a, b
	var olds, s = 1, 0
	var oldt, t = 0, 1

	for r != 0 {
		var quotient = oldr / r
		oldr, r = r, oldr-quotient*r
		olds, s = s, olds-quotient*s
		oldt, t = t, oldt-quotient*t
	}
	return olds, oldt, oldr
}

// A utiliser avec mod premier
func ModularInv(a, mod int) int {
	var x, _, _ = extanded_gcd(a, mod)
	return (x%mod + mod) % mod
}

//retourne une copie de pol
func Copy(pol *PolyInt) PolyInt {
	var coefs = make([]int, len(pol.Coefs), cap(pol.Coefs))
	copy(coefs, pol.Coefs)
	var res = NewPolyInt(coefs, pol.mod)
	return res
}

// modifie pol en le multipliant par l'entier f
func (pol *PolyInt) PolyScale(f int) {
	for i, c := range pol.Coefs {
		pol.Coefs[i] = (c * f)
	}
}

// renvoie f*pol
func PolyScale(pol *PolyInt, f int) PolyInt {
	res := Copy(pol)
	res.PolyScale(f)
	return res
}

// modifie pol en le multipliant par l'entier f, modulo pol.moy
func (pol *PolyInt) PolyScaleMod(f int) {
	pol.PolyScale(f)
	pol.Mod(pol.mod)
}

// renvoie f*pol, modulo pol.moy
func PolyScaleMod(pol *PolyInt, f int) PolyInt {
	var res = Copy(pol)
	res.PolyScaleMod(f)
	return res
}

// ajoute v à u
func (u *PolyInt) PolyAdd(v *PolyInt) *PolyInt {
	if u.Degree() >= v.Degree() {
		for i := range v.Coefs {
			u.Coefs[i] += v.Coefs[i]
		}
	} else {
		fmt.Println(u.Degree(), v.Degree())
		for i := range u.Coefs {
			u.Coefs[i] += u.Coefs[i]
		}
	}
	u.deflate()
	return u
}

// renvoie u + v
func PolyAdd(u, v *PolyInt) PolyInt {
	var res = Copy(u)
	res.PolyAdd(v)
	return res
}

// ajoute v à u, modulo u.mod
func (u *PolyInt) PolyAddMod(v *PolyInt) (*PolyInt, error) {
	mod := u.mod
	if v.mod != mod {
		return u, errors.New("adding polynomials with different modulus, returning first pol")
	} else {
		u.PolyAdd(v)
		u.Mod(mod)
		return u, nil
	}
}

// retourne u + v, modulo u.mod
func PolyAddMod(u, v *PolyInt) (PolyInt, error) {
	mod := u.mod
	if v.mod != mod {
		return *u, errors.New("adding polynomials with different modulus, returning first pol")
	} else {
		var res = Copy(u)
		res.PolyAddMod(v)
		return res, nil
	}
}

// renvoie une slice contenant n occurences de l'entier c
func ConstantIntSlice(n uint, c int) []int {
	var res []int
	for i := 1; i <= n; i++ {
		res = append(res, c)
	}
	return res
}

// multiplie pol par x^n
func (pol *PolyInt) PolyIntShift(n uint) *PolyInt {
	shift := ConstantIntSlice(n, 0)
	pol.Coefs = append(shift, pol.Coefs...)
	return pol
}

// revoie pol * x^n
func PolyIntShift(pol *PolyInt, n uint) PolyInt {
	//var res = Copy(pol)
	//res.PolyIntShift(n)
	return Copy(pol).PolyIntShift(n)
}

// A REECRIRE
// FAIRE UNE VARIANTE POUR DES POLYNOMES
// ICI, ON A QUE LA COMPARAISON DE LEURS COEFFICIENTS
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

// retourne le produit de u et v
func PolyIntMultMod(u, v *PolyInt) PolyInt {

	var mod = u.mod
	var du, dv = u.Degree(), v.Degree()
	var dprod = du + dv
	var product = []int{}
	var newCoeff int

	for i := 0; i <= dprod; i++ {
		newCoeff = 0
		for iu := 0; iu <= i; iu++ {
			if iu <= du && i-iu <= dv {
				newCoeff += (u.Coefs[iu] * v.Coefs[i-iu]) % mod
			}
		}
		product = append(product, newCoeff)
	}
	return NewPolyInt(product, mod)
}

// ajoute v à u sans deflate, modulo u.mod : pour faciliter la division
func (u *PolyInt) specialAdd(v *PolyInt) *PolyInt {
	if u.Degree() >= v.Degree() {
		for i := range v.Coefs {
			u.Coefs[i] += v.Coefs[i]
		}
	} else {
		for i := range v.Coefs {
			u.Coefs[i] += u.Coefs[i]
		}
	}
	u.Mod(u.mod)
	return u
}

// retourne le quotient et le reste de la division de u par v, modulo u.mod
func PolyIntDivMod(u, v *PolyInt) (PolyInt, PolyInt, error) {
	if Equal(v.Coefs, []int{0}) {
		return *u, *v, errors.New("division by 0, returnin u and v")
	} else {

		if u.Degree() < v.Degree() {
			return NewPolyInt([]int{0}, u.mod), *u, nil
		} else {
			var mod = u.mod
			var newu = Copy(u)
			var degNewu = newu.Degree()
			var dv = v.Degree()
			var quotient = []int{}
			var newCoef int
			var scaledv PolyInt

			for degNewu >= dv {

				newCoef = (newu.Coefs[degNewu] * ModularInv(v.Coefs[dv], mod)) % mod
				quotient = append([]int{newCoef}, quotient...)

				scaledv = PolyScaleMod(v, -newCoef)
				scaledv.PolyIntShift(uint(degNewu - dv))

				newu.specialAdd(&scaledv)
				degNewu = degNewu - 1
			}

			newu.deflate()
			polQuotient := NewPolyInt(quotient, mod)
			return polQuotient, newu, nil
			//
		}
	}
}
