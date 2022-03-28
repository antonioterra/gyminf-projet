package ring

// RING ELEMENT
/* type Ring struct {
	CMod *big.Int  // coefs modulus
	PMod poly.Poly // pol modulus
} */

// returns a ring with coef modulus mod and pol modulus X^N + 1
/* func NewRing(mod *big.Int, N int) Ring {
	pol := poly.ZeroPoly(N)
	pol.Coefs[0].Set(big.NewInt(1))
	pol.Coefs[N].Set(big.NewInt(1))
	return Ring{CMod: mod, PMod: pol}
} */

/* // takes the mod of the coefs of pol (in place)
func (r *Ring) TakeCoefMod(pol *poly.Poly) *poly.Poly {
	pol.TakeCoefMod(r.CMod)
	return pol
} */

// changes pol to [pol % PMod] (in place) and returns pol
// DOES NOT AND SHOULD NOT TAKE THE MOD OF THE COEFS !!!
/* func (r *Ring) ToRing(pol *poly.Poly) *poly.Poly {
	_, reste := poly.PolyDiv(*pol, r.PMod)
	pol.Coefs = reste.Coefs
	return pol
} */

/* //returns the sum in ring r of pols u and v WITHOUT REDUCING THE COEF TO MODULUS
func (r *Ring) RAdd(u, v poly.Poly) poly.Poly {
	res := poly.Add(u, v)
	//res.TakeCoefMod(r.CMod)
	return res
} */

// returns the ring-product of u and v WITHOUT REDUCING THE COEF TO MODULUS
/* func (r *Ring) RMult(u, v poly.Poly) poly.Poly {
	pro := poly.Mult(u, v)
	_, pro = poly.PolyDiv(pro, r.PMod)
	pro.Deflate()
	//pro.TakeCoefMod(r.CMod)s
	return pro
} */
