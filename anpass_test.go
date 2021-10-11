package anpass

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"reflect"
	"strings"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func nearby(a, b, eps float64) bool {
	if math.Abs(a-b) > eps {
		return false
	}
	return true
}

func eql(a, b []float64, eps float64) bool {
	for i := range a {
		if !nearby(a[i], b[i], eps) {
			fmt.Printf("problem with value %d, diff = %g\n",
				i, a[i]-b[i])
			return false
		}
	}
	return true
}

func TestEval(t *testing.T) {
	want := []float64{
		1.00000, -0.00500, -0.00500, 0.00003, 0.00003, 0.00003,
		0.00010, -0.00000, -0.00000, -0.00000, -0.00000, -0.00000,
		-0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
		0.00000, 0.00000, 0.00000, 0.00000,
	}
	disps, _, exps, _, _ := ReadInput("testfiles/anpass.in")
	// nvbl, nunk := exps.Dims()
	_, nunk := exps.Dims()
	coeffs := make([]float64, nunk)
	got := make([]float64, len(want))
	i := 0
	for j := range coeffs {
		coeffs[j] = 1
		got[j] = Eval(disps.RawRowView(i), coeffs, exps)
		coeffs[j] = 0
	}
	// it only prints 5 sig figs so that's as close as I can check
	if !eql(got, want, 1e-5) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func loadSlice(filename string) (ret []float64) {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		if fields := strings.Fields(scanner.Text()); len(fields) > 0 {
			ret = append(ret, toFloat(fields...)...)
		}
	}
	return
}

func TestFit(t *testing.T) {
	disps, energies, exps, _, _ := ReadInput("testfiles/anpass.in")
	got, _ := Fit(disps, energies, exps)
	want := loadSlice("testfiles/coeffs.matrix")
	if !eql(got.RawMatrix().Data, want, 9e-11) {
		t.Error()
		fmt.Printf("%20s%20s%20s\n", "GOT", "WANT", "DELTA")
		for i, v := range got.RawMatrix().Data {
			fmt.Printf("%20.12f%20.12f%20.12E\n",
				v, want[i], want[i]-v)
		}
	}
}

func TestMake9903(t *testing.T) {
	disps, energies, exps, _, _ := ReadInput("testfiles/anpass.in")
	coeffs, _ := Fit(disps, energies, exps)
	got := Make9903(io.Discard, coeffs, exps)
	want := []FC{
		{[4]int{0, 0, 0, 0}, 0.000000000009},
		{[4]int{1, 0, 0, 0}, 0.000388752719},
		{[4]int{2, 0, 0, 0}, 0.000035616298},
		{[4]int{1, 1, 0, 0}, 8.358958979225},
		{[4]int{2, 1, 0, 0}, 0.364209957363},
		{[4]int{2, 2, 0, 0}, 0.705550725579},
		{[4]int{3, 3, 0, 0}, 8.560855485923},
		{[4]int{1, 1, 1, 0}, -41.630626601721},
		{[4]int{2, 1, 1, 0}, -0.611038496965},
		{[4]int{2, 2, 1, 0}, -0.447311911909},
		{[4]int{2, 2, 2, 0}, -0.701538295617},
		{[4]int{3, 3, 1, 0}, -41.476685904878},
		{[4]int{3, 3, 2, 0}, 0.392862290084},
		{[4]int{1, 1, 1, 1}, 181.917491654261},
		{[4]int{2, 1, 1, 1}, -0.290752353552},
		{[4]int{2, 2, 1, 1}, 0.372034224849},
		{[4]int{2, 2, 2, 1}, 1.034528413172},
		{[4]int{2, 2, 2, 2}, -0.656556198243},
		{[4]int{3, 3, 1, 1}, 182.205477416855},
		{[4]int{3, 3, 2, 1}, -1.231659580792},
		{[4]int{3, 3, 2, 2}, -0.821069269704},
		{[4]int{3, 3, 3, 3}, 183.625347325953},
	}
	for i, fc := range got {
		if !reflect.DeepEqual(fc.Coord, want[i].Coord) {
			t.Errorf("got %v, wanted %v at index %d\n",
				fc.Coord, want[i].Coord, i)
		}
		if !nearby(fc.Val, want[i].Val, 5e-7) {
			t.Errorf("got %v, wanted %v at index %d\n",
				fc.Val, want[i].Val, i)
		}
	}
}

func compFC(got, want []FC, eps float64) bool {
	for i, fc := range got {
		if !reflect.DeepEqual(fc.Coord, want[i].Coord) {
			fmt.Printf("problem with fc %d =>\n%v\n%v", i, fc, want[i])
			return false
		}
		if !nearby(fc.Val, want[i].Val, eps) {
			fmt.Printf("problem with value %d, diff = %g\n", i,
				fc.Val-want[i].Val)
			fmt.Printf("got %v, wanted %v\n", fc.Val, want[i].Val)
			return false
		}
	}
	return true
}

func deepError(t *testing.T, a, b interface{}) {
	if !reflect.DeepEqual(a, b) {
		t.Errorf("got %v, wanted %v\n", a, b)
	}
}

func TestReadInput(t *testing.T) {
	tests := []struct {
		infile   string
		disps    []float64
		energies []float64
		exps     []float64
		biases   []float64
		stat     bool
	}{
		{
			infile:   "testfiles/anpass.in",
			disps:    loadSlice("testfiles/anpass.disps"),
			energies: loadSlice("testfiles/anpass.nrg"),
			exps:     loadSlice("testfiles/anpass.exp"),
			biases:   []float64{0, 0, 0, 0},
			stat:     false,
		},
		{
			infile:   "testfiles/anpass2.in",
			disps:    loadSlice("testfiles/anpass2.disps"),
			energies: loadSlice("testfiles/anpass2.nrg"),
			exps:     loadSlice("testfiles/anpass.exp"),
			biases: []float64{
				-0.000045311426,
				-0.000027076533,
				0.000000000000,
				-0.000000002131,
			},
			stat: true,
		},
	}
	for _, test := range tests {
		disps, energies, exps, biases, stat := ReadInput(test.infile)
		deepError(t, disps.RawMatrix().Data, test.disps)
		deepError(t, energies, test.energies)
		deepError(t, exps.RawMatrix().Data, test.exps)
		deepError(t, biases, test.biases)
		deepError(t, stat, test.stat)
	}
}

func TestBias(t *testing.T) {
	disps := mat.NewDense(3, 4, []float64{
		0.001, 0.002, 0.003, 0.004,
		0.005, 0.006, 0.007, 0.008,
		0.009, 0.010, 0.011, 0.012,
	})
	energies := []float64{10, 20, 30}
	biases := []float64{
		0.001, 0.002, 0.003, 0.004, 5,
	}
	gdisp, gnrg := Bias(disps, energies, biases)
	wdisp := []float64{
		0.000, 0.000, 0.000, 0.000,
		0.004, 0.004, 0.004, 0.004,
		0.008, 0.008, 0.008, 0.008,
	}
	wnrg := []float64{5, 15, 25}
	deepError(t, gdisp.RawMatrix().Data, wdisp)
	deepError(t, gnrg, wnrg)
}

func TestPrintResiduals(t *testing.T) {
	disps, energies, exps, _, _ := ReadInput("testfiles/anpass.in")
	coeffs, fn := Fit(disps, energies, exps)
	got := PrintResiduals(io.Discard, coeffs, fn, energies)
	want := 2.73700953e-18
	if !nearby(got, want, 1e-21) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func BenchmarkFit(b *testing.B) {
	disps, energies, exps, _, _ := ReadInput("testfiles/anpass.in")
	for n := 0; n < b.N; n++ {
		Fit(disps, energies, exps)
	}
}
