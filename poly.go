package poly

import (
	"math/big"
)

// {1, 2, 3} correspond à 1 + 2x + 3x^2
type Poly struct {
	Coefs []*big.Int
}

// retourne le polynôme dont les coefficients sont donnés par <coefs>
// {1, 2, 3} correspond à 1 + 2x + 3x^2
func NewPoly(coefs []*big.Int) Poly {
	var poly = Poly{Coefs: coefs}
	return poly
}

// Returns a slice containing <n> pointers to the bigFloat <c>
func ConstantSlice(n int, c big.Int) []*big.Int {
	var res = []*big.Int{}
	for i := 1; i <= n; i++ {
		res = append(res, &c)
	}
	return res
}

// shifts the coefs of the reciever <pol> by <n> positions towards biggest degree
func (pol *Poly) PolyShift(n int) *Poly {
	shift := ConstantSlice(n, *big.NewInt(0.))
	pol.Coefs = append(shift, pol.Coefs...)
	return pol
}

// returns the zero polynomial of fake degrée <n>, i.e. : 0 + 0x + ... + 0x^n
func ZeroPoly(n int) Poly {
	coefs := make([]*big.Int, n+1)
	for i := range coefs {
		coefs[i] = big.NewInt(0)
	}

	return NewPoly(coefs)
}

// retourne le degré du polynome reciever <pol>
func (poly *Poly) Degree() int {
	deg := 0
	for i := range poly.Coefs {
		if poly.Coefs[len(poly.Coefs)-1-i].Cmp(big.NewInt(0)) != 0 {
			deg = len(poly.Coefs) - i - 1
			break
		}
	}
	return deg
}

// retourne une copie du polynôme <pol>
func Copy(pol Poly) Poly {
	coefs := make([]*big.Int, len(pol.Coefs), cap(pol.Coefs))
	for i := 0; i < len(pol.Coefs); i++ {
		coefs[i] = big.NewInt(0)
		coefs[i].Set(pol.Coefs[i])
	}
	return NewPoly(coefs)
}

// enlève du reciever les coefficients nuls de plus haut degré
// ex : 1 + 2x + 0x^2 + 0x^3 -> 1 + 2x
func (poly *Poly) Deflate() {
	poly.Coefs = poly.Coefs[:poly.Degree()+1]
}

// retourne (MODIFIE ?) le polynôme correspondant au produit de <poly> par le big.Int <f>
// IN PLACE ? PLUS CLAIR. A VERIRIER.
func (poly *Poly) Scale(f *big.Int) *Poly {
	for i, c := range poly.Coefs {
		poly.Coefs[i].Mul(c, f)
	}
	return poly
}

// revoie le scaling du Poly <pol> par le big.Int <f>
func Scale(pol Poly, f *big.Int) Poly {
	coefs := make([]*big.Int, len(pol.Coefs))
	for i, c := range pol.Coefs {
		coefs[i] = new(big.Int).Mul(c, f)
	}
	return NewPoly(coefs)
}

// revoie la division du Poly <poly> par le bit.Int <f> avec arrondi des coefficients
func ScaleDiv(pol Poly, f *big.Int) Poly {
	coefs := make([]*big.Int, len(pol.Coefs))
	for i := range pol.Coefs {
		coefFloat := new(big.Float).SetInt(pol.Coefs[i])
		fFloat := new(big.Float).SetInt(f)
		coefFloat.Quo(coefFloat, fFloat)
		coefs[i] = new(big.Int)
		coefs[i], _ = coefFloat.Int(nil)
	}
	return NewPoly(coefs)
}

// retourne le somme des polynômes <u> et <v>
func Add(u, v Poly) Poly {

	var res Poly
	if u.Degree() >= v.Degree() {

		coefs := make([]*big.Int, len(u.Coefs))

		for i := 0; i <= v.Degree(); i++ { //, c := range v.Coefs { // ne marche pas si v n'est pas deflaté et il n'est pas souhaitable qu'il le soit
			coefs[i] = new(big.Int).Add(u.Coefs[i], v.Coefs[i])

		}
		for i := v.Degree() + 1; i < len(u.Coefs); i++ {
			coefs[i] = new(big.Int).Set(u.Coefs[i])
		}
		res = NewPoly(coefs)
	} else {
		//fmt.Println("here !")
		coefs := make([]*big.Int, len(v.Coefs))
		for i := 0; i <= u.Degree(); i++ { //, c := range u.Coefs {// ne marche pas si u n'est pas deflaté et il n'est pas souhaitable qu'il le soit
			coefs[i] = new(big.Int).Add(v.Coefs[i], u.Coefs[i])
		}
		for i := u.Degree() + 1; i <= v.Degree(); i++ {
			coefs[i] = new(big.Int).Set(v.Coefs[i])
		}
		res = NewPoly(coefs)
	}
	return res
}

// retourne le quotient et le reste de la division du poly.Poly <u> par le poly.Poly <v>
//AJTOUER GESTION DES ERREUR SUR DIVISION PAR POLYNOME NUL
func PolyDiv(u, v Poly) (Poly, Poly) {

	if u.Degree() < v.Degree() {
		return NewPoly([]*big.Int{big.NewInt(0)}), u
	} else {
		newu := Copy(u)
		degNewu := newu.Degree()
		dv := v.Degree()
		quoDeg := u.Degree() - v.Degree() + 1
		quotient := make([]*big.Int, quoDeg)
		newCoef := big.NewInt(0)
		var scaledv Poly

		i := 0
		for degNewu >= dv {
			newCoef.Quo(newu.Coefs[degNewu], v.Coefs[dv])
			quotient[quoDeg-i-1] = new(big.Int).Set(newCoef)

			scaledv = Scale(v, newCoef.Neg(newCoef))
			scaledv.PolyShift(degNewu - dv)
			newu = Add(newu, scaledv)

			degNewu = degNewu - 1
			i++

		}
		return NewPoly(quotient), newu
	}

}

// retourne le produit des poly.Poly <u> et <v>
func Mult(u, v Poly) Poly {

	var du, dv = u.Degree(), v.Degree()
	var dprod = du + dv
	product := make([]*big.Int, dprod+1)

	for i := 0; i <= dprod; i++ {
		newCoeff := big.NewInt(0)
		for iu := 0; iu <= i; iu++ {
			if iu <= du && i-iu <= dv {
				prod := new(big.Int).Mul(u.Coefs[iu], v.Coefs[i-iu])

				newCoeff.Add(newCoeff, prod)
			}

		}
		product[i] = new(big.Int).Set(newCoeff)
	}
	return NewPoly(product)
}

// Réduit in place les coefficients du reciever <pol> modulo le big.Int <mod>
func (pol *Poly) TakeCoefMod(mod *big.Int) *Poly {
	q, m := new(big.Int), new(big.Int)
	halfMod := new(big.Int).Div(mod, big.NewInt(2))
	for i, c := range pol.Coefs {
		q.DivMod(c, mod, m)
		absC := new(big.Int).Abs(c)
		if absC.Cmp(halfMod) != -1 { //il y a peut-être une petite imprécision ici : faut avoir un intervalle semi-ouvert
			if m.Cmp(halfMod) == -1 {
				pol.Coefs[i].Set(m)
			} else {
				pol.Coefs[i].Sub(m, mod)
			}
		}
	}
	return pol
}

//Retourne le reste de la division du produit polynomial (<u>*<v>) par le polynôme <mod>
func MultMod(u, v, mod Poly) Poly {
	pro := Mult(u, v)
	_, pro = PolyDiv(pro, mod)
	pro.Deflate()
	return pro
}
