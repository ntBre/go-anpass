package anpass

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/mat"
)

// threshold for considering values zero
const THR = 1e-10

// Global settings
var (
	Quiet bool
	Debug bool
	MAXIT = 100
)

type FC struct {
	Coord [4]int
	Val   float64
}

func toFloat(strs ...string) []float64 {
	ret := make([]float64, len(strs))
	var err error
	for i := range strs {
		ret[i], err = strconv.ParseFloat(strs[i], 64)
		if err != nil {
			panic(err)
		}
	}
	return ret
}

// ReadInput reads an anpass input file and returns the displacements,
// energies, and exponents as []float64s. exps could be integers, but
// you want them as floats for use in math.Pow
func ReadInput(filename string) (disps *mat.Dense, energies []float64,
	exps [][]int, biases []float64, stationary bool) {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	scanner := bufio.NewScanner(f)
	var (
		line      string
		fields    []string
		dispSlice []float64
		expsSlice []int
		ndisps    int
	)
	var handler func(string)
	dispHandler := func(line string) {
		fields = strings.Fields(line)
		ndisps = len(fields) - 1
		dispSlice = append(dispSlice, toFloat(fields[:ndisps]...)...)
		energies = append(energies, toFloat(fields[ndisps])...)
	}
	unkHandler := func(line string) {
		fields = strings.Fields(line)
		for _, d := range fields {
			v, _ := strconv.Atoi(d)
			expsSlice = append(expsSlice, v)
		}
	}
	statHandler := func(line string) {
		fields = strings.Fields(line)
		biases = append(biases, toFloat(fields...)...)
	}
	for scanner.Scan() {
		line = scanner.Text()
		switch {
		case len(line) > 0 && line[0] == '!':
			continue
		case strings.Contains(line, "("):
			handler = dispHandler
		case strings.Contains(line, "UNKNOWNS"):
			handler = nil
		case strings.Contains(line, "FUNCTION"):
			handler = unkHandler
		case strings.Contains(line, "END OF DATA"):
			handler = nil
		case strings.Contains(line, "STATIONARY POINT"):
			stationary = true
			handler = statHandler
		case handler != nil:
			handler(line)
		}
	}
	if biases == nil {
		biases = make([]float64, ndisps+1)
	}
	disps = mat.NewDense(len(dispSlice)/ndisps, ndisps, dispSlice)
	exps = Reshape(ndisps, len(expsSlice)/ndisps, expsSlice)
	return
}

// Dims returns the number of rows and cols in m, assuming each row in m has the
// same length as the first
func Dims(m [][]int) (rows, cols int) {
	return len(m), len(m[0])
}

func Reshape(rows, cols int, v []int) [][]int {
	ret := make([][]int, rows)
	for i := 0; i < rows; i++ {
		ret[i] = v[cols*i : cols*i+cols]
	}
	return ret
}

func Bias(disps *mat.Dense, energies, biases []float64) (*mat.Dense, []float64) {
	r, c := disps.Dims()
	ret := mat.NewDense(r, c, nil)
	enew := make([]float64, len(energies))
	last := len(biases)
	rbias := biases[:last-1]
	ebias := biases[last-1]
	for i, e := range energies {
		enew[i] = e - ebias
		for j := range rbias {
			ret.Set(i, j, disps.At(i, j)-rbias[j])
		}
	}
	return ret, enew
}

func Eval(Xi, C []float64, exps [][]int) float64 {
	var (
		prod, sum float64
		Ejk       int
	)
	for k := range C {
		prod = C[k]
		if math.Abs(prod) < THR {
			continue
		}
		for j := range Xi {
			Ejk = exps[j][k]
			if Ejk != 0 {
				prod *= math.Pow(Xi[j], float64(Ejk))
			}
		}
		sum += prod
	}
	return sum
}

func PrintMat(out io.Writer, matr mat.Matrix) {
	r, c := matr.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			fmt.Fprintf(out, "%12.5f", matr.At(i, j))
		}
		fmt.Fprint(out, "\n")
	}
}

func PrintVec(out io.Writer, v []float64) {
	for i := range v {
		if i > 0 && (i%6) == 0 {
			fmt.Fprint(out, "\n")
		}
		fmt.Fprintf(out, "%20.12E", v[i])
	}
	fmt.Fprint(out, "\n")
}

// Fit determines the coefficient vector using ordinary least squares
// and returns the solution vector along with the matrix describing
// the function
func Fit(disps *mat.Dense, energies []float64, exps [][]int) (
	soln, fn *mat.Dense) {
	_, coeffs := Dims(exps)
	pts := len(energies)
	X := mat.NewDense(pts, coeffs, nil)
	var (
		prod float64
		xijs []float64
		Ejk  int
	)
	tmp := make([]float64, coeffs)
	for i := 0; i < pts; i++ {
		xijs = disps.RawRowView(i)
		for k := 0; k < coeffs; k++ {
			prod = 1.0
			for j, xij := range xijs {
				Ejk = exps[j][k]
				switch Ejk {
				case 0.0:
				case 1.0:
					prod *= xij
				case 2.0:
					prod *= xij * xij
				case 3.0:
					prod *= xij * xij * xij
				case 4.0:
					prod *= xij * xij * xij * xij
				default:
					panic("didn't match exponent")
				}
			}
			tmp[k] = prod
		}
		X.SetRow(i, tmp)
	}
	y := mat.NewDense(pts, 1, energies)
	var XTX mat.Dense
	XTX.Mul(X.T(), X)
	var inv mat.Dense
	err := inv.Inverse(&XTX)
	if err != nil {
		if strings.Contains(err.Error(), "Inf") {
			panic(err)
		}
		if !Quiet {
			fmt.Fprintf(os.Stderr, "WARNING: %v\n", err)
		}
	}
	var mul mat.Dense
	mul.Mul(&inv, X.T())
	var sol mat.Dense
	sol.Mul(&mul, y)
	return &sol, X
}

// PrintResiduals computes and prints the residual for each point and
// returns the sum of squared residuals
func PrintResiduals(w io.Writer, x, A *mat.Dense, energies []float64) (
	sum float64) {
	fmt.Fprintf(w, "%5s%20s%20s%20s\n",
		"POINT", "COMPUTED", "OBSERVED", "RESIDUAL")
	var prod mat.Dense
	prod.Mul(A, x)
	var comp, resi float64
	for i, obsv := range energies {
		comp = prod.At(i, 0)
		resi = comp - obsv
		fmt.Fprintf(w, "%5d%20.12f%20.12f%20.8E\n",
			i+1, comp, obsv, resi)
		sum += resi * resi
	}
	fmt.Fprintf(w, "WEIGHTED SUM OF SQUARED RESIDUALS IS %17.8E\n", sum)
	return
}

// MakeFCs is a a copy of Make9903 without the Writer. It just returns the slice
// of force constants instead of formatting them for output
func MakeFCs(coeffs, exps *mat.Dense) (ret []FC) {
	c, r := exps.Dims()
	for i := 0; i < r; i++ {
		ifact := 1
		var ictmp [4]int
		iccount := 0
		for j := c - 1; j >= 0; j-- {
			iexpo := int(exps.At(j, i))
			switch iexpo {
			case 2:
				ifact *= 2
			case 3:
				ifact *= 6
			case 4:
				ifact *= 24
			}
			if iexpo > 0 {
				for k := 0; k < iexpo; k++ {
					ictmp[iccount+k] = j + 1
				}
				iccount += iexpo
			}
		}
		ffcc := coeffs.At(i, 0) * float64(ifact) * 4.359813653e0
		ret = append(ret, FC{ictmp, ffcc})
	}
	return
}

// Make9903 is a helper function for writing fort.9903 files, but it
// also returns the force constants in a more usable format for
// testing
func Make9903(w io.Writer, coeffs *mat.Dense, exps [][]int) (ret []FC) {
	c, r := Dims(exps)
	for i := 0; i < r; i++ {
		ifact := 1
		var ictmp [4]int
		iccount := 0
		for j := c - 1; j >= 0; j-- {
			iexpo := exps[j][i]
			ifact *= []int{1, 1, 2, 6, 24}[iexpo]
			if iexpo > 0 {
				for k := 0; k < iexpo; k++ {
					ictmp[iccount+k] = j + 1
				}
				iccount += iexpo
			}
		}
		ffcc := coeffs.At(i, 0) * float64(ifact) * 4.359813653e0
		for _, f := range ictmp {
			fmt.Fprintf(w, "%5d", f)
		}
		fmt.Fprintf(w, "%20.12f\n", ffcc)
		ret = append(ret, FC{ictmp, ffcc})
	}
	return
}

func Write9903(filename string, coeffs *mat.Dense, exps [][]int) []FC {
	f, err := os.Create(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	return Make9903(f, coeffs, exps)
}

func Grad(x []float64, coeffs *mat.Dense, exps [][]int) (grd []float64) {
	var sum float64
	nvbl, nunk := Dims(exps)
	grd = make([]float64, nvbl)
	for i := 0; i < nvbl; i++ {
		sum = 0.0
		for j := 0; j < nunk; j++ {
			fij := float64(exps[i][j])
			coj := coeffs.At(j, 0) * fij
			if math.Abs(coj) < THR {
				continue
			}
			if exps[i][j] != 1 {
				coj *= math.Pow(x[i], fij-1)
			}
			for k := 0; k < nvbl; k++ {
				ekj := exps[k][j]
				if k != i && ekj != 0 {
					coj *= math.Pow(x[k], float64(ekj))
				}
			}
			sum += coj
		}
		grd[i] = sum
	}
	return
}

func Hess(x []float64, coeffs *mat.Dense, exps [][]int) *mat.SymDense {
	nvbl, nunk := Dims(exps)
	var (
		il  int
		sum float64
	)
	coeffSlice := coeffs.RawMatrix().Data
	hess := make([]float64, nvbl*(nvbl-1))
	var (
		coj           float64
		eij, elj, ekj int
		fij, flj      float64
	)
	for i := 0; i < nvbl; i++ {
		for l := 0; l <= i; l++ {
			sum = 0.0
			if i != l { // => off-diagonal
				for j := 0; j < nunk; j++ {
					coj = coeffSlice[j]
					eij = exps[i][j]
					elj = exps[l][j]
					fij = float64(eij)
					flj = float64(elj)
					coj *= fij * flj
					if math.Abs(coj) < THR {
						continue
					}
					if eij != 1 {
						coj *= math.Pow(x[i], fij-1)
					}
					if elj != 1 {
						coj *= math.Pow(x[l], flj-1)
					}
					for k := 0; k < nvbl; k++ {
						if k != i && k != l {
							if ekj = exps[k][j]; ekj != 0 {
								coj *= math.Pow(x[k], float64(ekj))
							}
						}
					}
					sum += coj
				}
			} else { // => diagonal
				for j := 0; j < nunk; j++ {
					coj = coeffSlice[j]
					eij = exps[i][j]
					fij = float64(eij)
					coj *= fij * (fij - 1)
					if math.Abs(coj) < THR {
						continue
					}
					if exps[i][j] != 2 {
						coj *= math.Pow(x[i], fij-2)
					}
					for k := 0; k < nvbl; k++ {
						if k != i {
							if ekj = exps[k][j]; ekj != 0 {
								coj *= math.Pow(x[k], float64(ekj))
							}
						}
					}
					sum += coj
				}
			}
			hess[il] = sum
			il++
		}
	}
	ret := mat.NewSymDense(nvbl, nil)
	il = 0
	for i := 0; i < nvbl; i++ {
		for j := 0; j <= i; j++ {
			ret.SetSym(i, j, hess[il])
			il++
		}
	}
	return ret
}

// Newton uses the Newton-Raphson method to find the roots of the
// equation given by coeffs and exps
func Newton(coeffs *mat.Dense, exps [][]int) []float64 {
	nvbl, _ := Dims(exps)
	x := make([]float64, nvbl)
	// if exceed 100, give up, too many iterations
	for iter := 0; iter < MAXIT; iter++ {
		grad := Grad(x, coeffs, exps)
		hess := Hess(x, coeffs, exps)
		var invHess mat.Dense
		err := invHess.Inverse(hess)
		if err != nil {
			fmt.Fprintf(os.Stderr, "WARNING: %v\n", err)
		}
		// Newton-Raphson update with step of 0.5
		var (
			iconv bool = true
			sum   float64
		)
		_, c := hess.Dims()
		gradMat := mat.NewDense(c, 1, grad)
		var del mat.Dense
		del.Mul(&invHess, gradMat)
		del.Scale(0.5, &del)
		for _, i := range del.RawMatrix().Data {
			if i > 1.1e-8 {
				sum = i
				iconv = false
				break
			}
		}
		if iconv {
			return x
		} else if Debug {
			fmt.Printf("ITERATION %5d UPDATE VECTOR, RES = %+10.5e\n",
				iter+1, sum)
			for _, i := range del.RawMatrix().Data {
				fmt.Printf("%12.8f", i)
			}
			fmt.Print("\n")
		}
		for i := range x {
			x[i] -= del.At(i, 0)
		}
	}
	panic("TOO MANY NEWTON-RAPHSON ITERATIONS")
}

// Characterize stationary point x by computing the Hessian and
// determining its eigenvalues
func Characterize(x []float64, coeffs *mat.Dense, exps [][]int) (
	evals []float64, evecs *mat.Dense, kind Stat) {
	h := Hess(x, coeffs, exps)
	var eig mat.EigenSym
	eig.Factorize(h, true)
	evals = eig.Values(nil)
	var prod int
	for _, v := range evals {
		if v < 0 {
			prod--
		} else if v > 0 {
			prod++
		}
	}
	switch l := len(evals); prod {
	case -l:
		kind = MAX
	case l:
		kind = MIN
	default:
		kind = STAT
	}
	var vecs mat.Dense
	eig.VectorsTo(&vecs)
	evecs = &vecs
	return
}

// CopyAnpass copies the anpass1 input file in infile to an anpass2
// input file in outfile, adding the longLine describing the found
// stationary point
func CopyAnpass(infile, outfile string, longLine []float64) {
	in, err := os.Open(infile)
	defer in.Close()
	if err != nil {
		panic(err)
	}
	out, err := os.Create(outfile)
	defer out.Close()
	if err != nil {
		panic(err)
	}
	scanner := bufio.NewScanner(in)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.Contains(line, "END OF DATA") {
			fmt.Fprintln(out, "STATIONARY POINT")
			for _, v := range longLine {
				fmt.Fprintf(out, "%20.12f", v)
			}
			fmt.Fprint(out, "\n")
		}
		fmt.Fprintln(out, line)
	}
}

// Run runs anpass: it computes the coefficients that fit disps, energies, and
// exps; it then calls Newton to locate the stationary point and evaluates the
// function at the stationary point.
func Run(w io.Writer, dir string, disps *mat.Dense, energies []float64,
	exps [][]int) (longLine []float64, fcs []FC, stationary bool) {
	coeffs, fn := Fit(disps, energies, exps)
	PrintResiduals(w, coeffs, fn, energies)
	fcs = Write9903(filepath.Join(dir, "fort.9903"), coeffs, exps)
	x := Newton(coeffs, exps)
	// characterize stationary point found by Newton
	evals, evecs, kind := Characterize(x, coeffs, exps)
	fmt.Fprintf(w, "\n%s\n", kind)
	// print long line and intder steps
	e := Eval(x, coeffs.RawMatrix().Data, exps)
	fmt.Fprintf(w, "WHERE ENERGY IS %20.12f\n", e)
	at := "AT"
	for i := range x {
		fmt.Fprintf(w, "%12s", at)
		fmt.Fprintf(w, "%18.10f\n", x[i])
		at = ""
	}
	longLine = append(x, e)
	for _, v := range longLine {
		fmt.Fprintf(w, "%20.12f", v)
	}
	fmt.Fprint(w, "\n")
	fmt.Fprintln(w, "EIGENVALUE(S) OF HESSIAN, STARTING WITH LOWEST")
	fmt.Fprint(w, "\n")
	for i, v := range evals {
		fmt.Fprintf(w, "EIGENVALUE %5d%20.10E\n", i+1, v)
		for _, v := range evecs.RawRowView(i) {
			fmt.Fprintf(w, "%16.8f", v)
		}
		fmt.Fprint(w, "\n")
	}
	return longLine, fcs, stationary
}

func PrintBias(w io.Writer, biases []float64) {
	last := len(biases)
	rbias := biases[:last-1]
	ebias := biases[last-1]
	fmt.Fprintf(w, "INITIAL GUESS AT STATIONARY POINT IS %20.12f\n",
		ebias)
	for _, r := range rbias {
		fmt.Fprintf(w, "%12.8f", r)
	}
	fmt.Fprint(w, "\n")
}
