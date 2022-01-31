package poly

// {1, 2, 3} correspond à 1 + 2x + 3x^2
type Poly struct {
	Coefs []float64
}

// retourne un polynôme deflaté
// retourne un polynôme deflaté
func NewPoly(coefs []float64) Poly {
	var poly = Poly{Coefs: coefs}
	//poly.Deflate()
	return poly
}

// vérifie si un polynôme est nul
func IsZero(poly *Poly) bool {
	var res = false
	for _, c := range poly.Coefs {
		if c == 0 {
			res = true
			break
		}

	}
	return res
}

func ConstantSlice(n int, c float64) []float64 {
	var res = []float64{}
	for i := 1; i <= n; i++ {
		res = append(res, c)
	}
	return res
}

func SliceShift(n int, slice *[]float64) []float64 {
	shift := ConstantSlice(n, 0)
	return append(shift, *slice...)
}

func PolyShift(n int, poly *Poly) Poly {
	var res = poly
	shift := ConstantSlice(n, 0)
	res.Coefs = append(shift, res.Coefs...)
	return *res
}

// retourne l'élément neutre de R[x]_n
func ZeroPoly(n int) Poly {
	var res = Poly{Coefs: []float64{}}
	for i := 0; i <= n; i++ {
		res.Coefs = append(res.Coefs, 0.0)
	}
	return res
}

// retourne le degré d'un polynôme
func (poly *Poly) Degree() int {
	var deg = 0
	for i, _ := range poly.Coefs {
		if poly.Coefs[len(poly.Coefs)-1-i] != 0 {
			deg = len(poly.Coefs) - i - 1
			break
		}
	}
	return deg
}

//retourne une copie de pol
func Copy(pol *Poly) Poly {
	coefs := make([]float64, len(pol.Coefs), cap(pol.Coefs))
	copy(coefs, pol.Coefs)
	var res = NewPoly(coefs)
	return res
}

func (poly *Poly) Deflate() {
	poly.Coefs = poly.Coefs[:poly.Degree()+1]
}

// multiplie un polynôme par un scalaire
func (poly *Poly) Scale(f float64) {
	for i, c := range poly.Coefs {
		poly.Coefs[i] = c * f
	}
}

// revoie le scaling du Poly pol par le float f
func Scale(pol Poly, f float64) Poly {
	res := Copy(&pol)
	for i, c := range pol.Coefs {
		res.Coefs[i] = c * f
	}
	return res
}

// retourne le somme de deux polynômes
func PolyAdd(u, v Poly) Poly {
	var res Poly
	if u.Degree() >= v.Degree() {
		coefs := make([]float64, len(u.Coefs))
		copy(coefs, u.Coefs)
		for i, c := range v.Coefs {
			coefs[i] += c
		}
		res = NewPoly(coefs)
	} else {
		coefs := make([]float64, len(v.Coefs))
		copy(coefs, v.Coefs)

		for i, c := range u.Coefs {
			coefs[i] += c
		}
		res = NewPoly(coefs)
	}
	return res
}

// multiplie pol par x^n
func (pol *Poly) PolyShift(n int) *Poly {
	shift := ConstantSlice(n, 0)
	pol.Coefs = append(shift, pol.Coefs...)
	return pol
}

// retourne le quotient et le reste d'une division polynomiale
func PolyDiv(u, v Poly) (Poly, Poly) {
	//AJTOUER GESTION DES ERREUR SUR DIVISION PAR POLYNOME NUL
	//u.Deflate()
	//v.Deflate()

	if u.Degree() < v.Degree() {
		return NewPoly([]float64{0}), u
	} else {
		newu := Copy(&u)
		degNewu := newu.Degree()
		dv := v.Degree()
		quotient := []float64{}
		var newCoef float64
		var scaledv Poly

		for degNewu >= dv {

			newCoef = newu.Coefs[degNewu] / v.Coefs[dv]
			quotient = append([]float64{newCoef}, quotient...)

			scaledv = Scale(v, -newCoef)
			scaledv.PolyShift(degNewu - dv)

			newu = PolyAdd(newu, scaledv)
			degNewu = degNewu - 1

		}
		return NewPoly(quotient), newu
	}

}

func Mult(u, v Poly) Poly {

	var du, dv = u.Degree(), v.Degree()
	var dprod = du + dv
	var product = []float64{}
	var newCoeff float64

	for i := 0; i <= dprod; i++ {
		newCoeff = 0
		for iu := 0; iu <= i; iu++ {
			if iu <= du && i-iu <= dv {
				newCoeff += (u.Coefs[iu] * v.Coefs[i-iu])
			}
		}
		product = append(product, newCoeff)
	}
	return NewPoly(product)
}