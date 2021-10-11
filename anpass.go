package anpass

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
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
	exps *mat.Dense, biases []float64, stationary bool) {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	scanner := bufio.NewScanner(f)
	var (
		line      string
		fields    []string
		dispSlice []float64
		expsSlice []float64
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
		expsSlice = append(expsSlice, toFloat(fields...)...)
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
	exps = mat.NewDense(ndisps, len(expsSlice)/ndisps, expsSlice)
	return
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

func Eval(Xi, C []float64, exps *mat.Dense) float64 {
	var prod, sum, Ejk float64
	for k := range C {
		prod = C[k]
		if math.Abs(prod) < THR {
			continue
		}
		for j := range Xi {
			Ejk = exps.At(j, k)
			if Ejk != 0.0 {
				prod *= math.Pow(Xi[j], Ejk)
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
func Fit(disps *mat.Dense, energies []float64, exps *mat.Dense) (
	soln, fn *mat.Dense) {
	_, coeffs := exps.Dims()
	pts := len(energies)
	X := mat.NewDense(pts, coeffs, nil)
	var prod, Ejk float64
	for i := 0; i < pts; i++ {
		for k := 0; k < coeffs; k++ {
			prod = 1.0
			for j, xij := range disps.RawRowView(i) {
				Ejk = exps.At(j, k)
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
			X.Set(i, k, prod)
		}
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
		fmt.Fprintf(os.Stderr, "WARNING: %v\n", err)
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

// Make9903 is a helper function for writing fort.9903 files, but it
// also returns the force constants in a more usable format for
// testing
func Make9903(w io.Writer, coeffs, exps *mat.Dense) (ret []FC) {
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
		for _, f := range ictmp {
			fmt.Fprintf(w, "%5d", f)
		}
		fmt.Fprintf(w, "%20.12f\n", ffcc)
		ret = append(ret, FC{ictmp, ffcc})
	}
	return
}

func Write9903(filename string, coeffs, exps *mat.Dense) []FC {
	f, err := os.Create(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	return Make9903(f, coeffs, exps)
}

func Grad(x []float64, coeffs, exps *mat.Dense) (grd []float64) {
	var sum float64
	nvbl, nunk := exps.Dims()
	grd = make([]float64, nvbl)
	for i := 0; i < nvbl; i++ {
		sum = 0.0
		for j := 0; j < nunk; j++ {
			coj := coeffs.At(j, 0) * exps.At(i, j)
			if math.Abs(coj) < THR {
				continue
			}
			if exps.At(i, j) != 1 {
				coj *= math.Pow(x[i], exps.At(i, j)-1)
			}
			for k := 0; k < nvbl; k++ {
				if k != i && exps.At(k, j) != 0 {
					coj *= math.Pow(x[k], exps.At(k, j))
				}
			}
			sum += coj
		}
		grd[i] = sum
	}
	return
}

func Hess(x []float64, coeffs, exps *mat.Dense) *mat.SymDense {
	nvbl, nunk := exps.Dims()
	var (
		ij  int
		sum float64
	)
	hess := make([]float64, nvbl*(nvbl-1))
	for i := 0; i < nvbl; i++ {
		for l := 0; l <= i; l++ {
			sum = 0.0
			if i != l { // => off-diagonal
				for j := 0; j < nunk; j++ {
					coj := coeffs.At(j, 0)
					coj *= exps.At(i, j) * exps.At(l, j)
					if math.Abs(coj) < THR {
						continue
					}
					if exps.At(i, j) != 1 {
						coj *= math.Pow(x[i],
							exps.At(i, j)-1)
					}
					if exps.At(l, j) != 1 {
						coj *= math.Pow(x[l],
							exps.At(l, j)-1)
					}
					for k := 0; k < nvbl; k++ {
						if k != i && k != l &&
							exps.At(k, j) != 0 {
							coj *= math.Pow(x[k],
								exps.At(k, j))
						}
					}
					sum += coj
				}
			} else { // => diagonal
				for j := 0; j < nunk; j++ {
					coj := coeffs.At(j, 0)
					coj *= exps.At(i, j)
					coj *= exps.At(i, j) - 1
					if math.Abs(coj) < THR {
						continue
					}
					if exps.At(i, j) != 2 {
						coj *= math.Pow(x[i],
							exps.At(i, j)-2)
					}
					for k := 0; k < nvbl; k++ {
						if k != i && exps.At(k, j) != 0 {
							coj *= math.Pow(x[k],
								exps.At(k, j))
						}
					}
					sum += coj
				}
			}
			hess[ij] = sum
			ij++
		}
	}
	ret := mat.NewSymDense(nvbl, nil)
	ij = 0
	for i := 0; i < nvbl; i++ {
		for j := 0; j <= i; j++ {
			ret.SetSym(i, j, hess[ij])
			ij++
		}
	}
	return ret
}

// Newton uses the Newton-Raphson method to find the roots of the
// equation given by coeffs and exps
func Newton(coeffs, exps *mat.Dense) []float64 {
	nvbl, _ := exps.Dims()
	x := make([]float64, nvbl)
	// if exceed 100, give up, too many iterations
	for iter := 0; iter < 100; iter++ {
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
		del := make([]float64, nvbl)
		for i := 0; i < nvbl; i++ {
			sum = 0.0
			for j := 0; j < nvbl; j++ {
				sum += 0.5 * invHess.At(i, j) * grad[j]
			}
			del[i] = sum
			if math.Abs(sum) > 1.1e-8 {
				iconv = false
			}
		}
		if iconv {
			return x
		} else if Debug {
			fmt.Printf("ITERATION %5d UPDATE VECTOR, RES = %12.9f\n",
				iter+1, sum)
			for i := range del {
				fmt.Printf("%12.8f", del[i])
			}
			fmt.Print("\n")
		}
		for i := range x {
			x[i] -= del[i]
		}
	}
	panic("TOO MANY NEWTON-RAPHSON ITERATIONS")
}

// Characterize stationary point x by computing the Hessian and
// determining its eigenvalues
func Characterize(x []float64, coeffs, exps *mat.Dense) (
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

func Run(w io.Writer, disps *mat.Dense, energies []float64, exps *mat.Dense) (
	longLine []float64, fcs []FC, stationary bool) {
	coeffs, fn := Fit(disps, energies, exps)
	PrintResiduals(w, coeffs, fn, energies)
	fcs = Write9903("fort.9903", coeffs, exps)
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
