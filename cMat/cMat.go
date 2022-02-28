package cMat

import (
	"fmt"
	"math"
	"math/cmplx"

	"gonum.org/v1/gonum/mat"
)

// structure sur laquelle seront définies de nouvelles opérations sur les matrices à coefficients complexes
type CMat struct {
	mat mat.CDense
}

func Copy(a *CMat) CMat {
	r, c := a.mat.Dims()
	data := make([]complex128, r*c)
	copy(data, a.GetData())
	return NewCMat(r, c, data)
}

func (a *CMat) Copy() CMat {
	r, c := a.mat.Dims()
	data := make([]complex128, r*c)
	copy(data, a.GetData())
	return NewCMat(r, c, data)
}

func (m *CMat) GetData() []complex128 {
	return m.mat.RawCMatrix().Data
}

func NewCMat(n, m int, data []complex128) CMat {
	ma := mat.NewCDense(n, m, data)
	return CMat{mat: *ma}
}

// retourne le # de ligne dans la matrice m sur lequel se trouve le ième coefficient
// compte des coefs commence à 0
// première ligne/col est ligne/col 0
func (ma *CMat) getRawOf(i int) int {
	m := ma.mat
	_, c := m.Dims()
	return int(math.Floor(float64(i) / float64(c)))
}

// retourne le # de colonne dans la matrice m sur lequel se trouve le ième coefficient
// compte des coefs commence à 0
// première ligne/col est ligne/col 0
func (ma *CMat) getColOf(i int) int {
	m := ma.mat
	_, c := m.Dims()
	return i % c
}

// AJOUTER GESTION DERREUR
// Retourne une copie de la ième colonne du reciever
// la première colonne est la colonne 0
func (ma *CMat) GetCol(i int) *CMat {

	r, c := ma.mat.Dims()
	resdata := make([]complex128, r)
	j := 0
	for k, coef := range ma.GetData() {
		if k%c == i {
			resdata[j] = coef
			j++
		}
	}
	res := NewCMat(r, 1, resdata)
	return &res
}

// AJOUTER GESTION DERREUR
// Retourne une copie de la ième ligne du reciever
// la première ligne est la ligne 0
func (ma *CMat) GetRaw(i int) *CMat {

	_, c := ma.mat.Dims()
	resdata := make([]complex128, c)
	copy(resdata, ma.GetData()[i*c:(i+1)*c])
	res := NewCMat(1, c, resdata)
	return &res
}

// retourne des matrices E et F telles que reciever = E + iF
func (ma *CMat) RI() (mat.Dense, mat.Dense) {
	m := ma.mat
	r, c := m.Dims()
	len := r * c
	rePart := make([]float64, len)
	imPart := make([]float64, len)
	originaldata := m.RawCMatrix().Data
	for k, coef := range originaldata {
		rePart[k] = real(coef)
		imPart[k] = imag(coef)
	}
	return *mat.NewDense(r, c, rePart), *mat.NewDense(r, c, imPart)
}

// retourne la matrice CMAT = E + iF
func Merge(E, F *mat.Dense) CMat {
	r, c := E.Dims()
	data := make([]complex128, r*c)
	for k := range E.RawMatrix().Data {
		data[k] = complex(E.RawMatrix().Data[k], F.RawMatrix().Data[k])
	}
	return NewCMat(r, c, data)
}

//retourne des matrices E et F telles que E + iF = inverse du reciever
// cf article sur inversion de matrices complexes
func (ma *CMat) Invert() (mat.Dense, mat.Dense) {
	m := ma.mat
	n, _ := m.Dims()
	n2 := int(math.Pow(float64(n), 2))
	data := make([]float64, 4*n2)
	for k, coef := range m.RawCMatrix().Data {
		baseIndex := ma.getRawOf(k)*n + k
		data[baseIndex] = real(coef)
		data[baseIndex+n] = imag(coef)
		data[baseIndex+2*n2] = -imag(coef)
		data[baseIndex+2*n2+n] = real(coef)
	}

	M := mat.NewDense(2*n, 2*n, data)
	var MInv mat.Dense
	MInv.Inverse(M)

	dataA, dataB := make([]float64, n2), make([]float64, n2)

	for k := 0; k < n; k++ {
		for l := 0; l < n; l++ {
			dataA[k*n+l] = MInv.At(k, l)
			dataB[k*n+l] = MInv.At(k, n+l)
		}
	}

	return *mat.NewDense(n, n, dataA), *mat.NewDense(n, n, dataB)
}

// Pretty print pour les CMat
func (ma *CMat) PPrint() {
	m := ma.mat
	r, c := m.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			fmt.Printf("%4.4f ", m.At(i, j))
			//fmt.Printf("%.2g", m.At(i, j))

		}
		fmt.Printf("\n")
	}
}

// Pretty print pour les mat.Dense
func PP(m *mat.Dense) {
	fc := mat.Formatted(m, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%.2g\n\n", fc)
}

// place la somme a + b dans le reciever
func (m *CMat) Add(a, b *CMat) *CMat {
	mR, mI := m.RI()
	aR, aI := a.RI()
	bR, bI := b.RI()

	mR.Add(&aR, &bR)
	mI.Add(&aI, &bI)

	*m = Merge(&mR, &mI)
	return m
}

// place le produit coef par coef a . b dans le reciever
// AJOUTER GESTION ERREURS SUIVANT LES TYPES
func (m *CMat) CoefWiseProd(a, b *CMat) *CMat {
	aR, aI := a.RI()
	bR, bI := b.RI()

	var s, t, u, v mat.Dense
	s.MulElem(&aI, &bI)
	s.Scale(-1, &s)
	t.MulElem(&aR, &bR)
	t.Add(&s, &t)

	u.MulElem(&aR, &bI)
	v.MulElem(&aI, &bR)
	u.Add(&u, &v)

	*m = Merge(&t, &u)
	return m
}

// mutiplie in place le reciever par c
func (m *CMat) Scale(c complex128) *CMat {
	ma := m.mat
	for i, k := range ma.RawCMatrix().Data {
		ma.RawCMatrix().Data[i] = k * c
	}
	return m
}

// place la différence a - b dans le reciever
func (m *CMat) Sub(a, b *CMat) *CMat {
	mR, mI := m.RI()
	aR, aI := a.RI()
	bR, bI := a.RI()
	bR.Scale(-1, &bR)
	bI.Scale(-1, &bI)

	mR.Add(&aR, &bR)
	mI.Add(&aI, &bI)

	*m = Merge(&mR, &mI)
	return m
}

// place le produit matriciel a*b dans le reciever
func (m *CMat) Mult(a, b *CMat) *CMat {
	mR, mI := m.RI()
	aR, aI := a.RI()
	bR, bI := b.RI()

	var aa, ab, ba, bb mat.Dense

	aa.Mul(&aR, &bR)
	bb.Mul(&aI, &bI)
	bb.Scale(-1, &bb)

	ab.Mul(&aR, &bI)
	ba.Mul(&aI, &bR)

	mR.Add(&aa, &bb)
	mI.Add(&ab, &ba)

	*m = Merge(&mR, &mI)
	return m
}

func (m *CMat) Transpose() {
	r, c := m.mat.Dims()
	data := make([]complex128, r*c)
	for i := 0; i < c; i++ {
		//fmt.Println("i = ", i)
		for j := 0; j < r; j++ {
			//fmt.Println("j = ", j)
			//fmt.Println("coef value = ", i*r+j)

			data[i*r+j] = m.mat.At(j, i)
			//fmt.Println("hop")
		}
	}
	m.mat = *mat.NewCDense(c, r, data)
}

func (m *CMat) Conjugate() *CMat {
	for k, coef := range m.GetData() {
		m.mat.RawCMatrix().Data[k] = cmplx.Conj(coef)
	}
	return m
}

//AJOUTER CONTROLE DES TAILLE DES VECTEURS
func DotProduct(a, b *CMat) complex128 {
	res := 0 + 0i
	dataa, datab := a.mat.RawCMatrix().Data, b.mat.RawCMatrix().Data
	for k, _ := range a.GetData() {
		res += dataa[k] * cmplx.Conj(datab[k])
	}
	return res
}

//Renvoie le carré de la norme
//Ajouter le controle de la dimension
func (m *CMat) SquaredNorm() complex128 {
	return DotProduct(m, m)
}

//Returns ||a-b||_infty
func (a *CMat) MaxNorm() float64 {
	ra, ca := a.mat.Dims()
	max := 0.
	for i := 0; i < ra; i++ {
		for j := 0; j < ca; j++ {
			challenger := cmplx.Abs(a.mat.At(i, j))
			if challenger > max {
				max = challenger
			}
		}
	}

	return max
}

//Renvoie la projection de a sur b en arrondissant les composantes
//Arrondi non random
//ne modifie ni a ni b
func ProjectOn(a, b *CMat) *CMat {
	coef := real(DotProduct(a, b))
	norm := real(b.SquaredNorm())
	//coef := DotProduct(a, b)
	//norm := b.SquaredNorm()
	res := Copy(b)
	//res.Scale(coef / norm)
	res.Scale(complex(math.Round(coef/norm), 0))
	return &res
}

//retourne la somme des projections de a sur les colonnes de b
func ProjectOnCols(a, b *CMat) *CMat {
	r, c := b.mat.Dims()

	res := NewCMat(r, 1, make([]complex128, r))

	var col, proj *CMat

	for k := 0; k < c; k++ {
		col = b.GetCol(k)
		proj = ProjectOn(a, col)
		res.Add(&res, proj)
	}

	return &res

}

//retourne la somme des projections de a sur les lignes de b
// sortie en tant que vecteur colonne
func ProjectOnRaws(a, b *CMat) *CMat {
	r, c := b.mat.Dims()
	res := NewCMat(c, 1, make([]complex128, c))

	var raw, proj *CMat

	for k := 0; k < r; k++ {
		raw = b.GetRaw(k)
		proj = ProjectOn(a, raw)
		proj.Transpose()
		res.Add(&res, proj)
	}

	return &res

}
