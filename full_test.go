package anpass

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
	"testing"
	"time"
)

func load9903(filename string) (ret []FC) {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		hold := new(FC)
		if len(fields) == 5 {
			for i := range fields[:4] {
				hold.Coord[i], _ = strconv.Atoi(fields[i])
			}
			hold.Val, _ = strconv.ParseFloat(fields[4], 64)
			ret = append(ret, *hold)
		}
	}
	return
}

// want is a fort.9903 file from running the old anpass twice on
// infile
func TestFull(t *testing.T) {
	Quiet = true
	tests := []struct {
		infile string
		want   string
		lline  []float64
		lineps float64
		eps    float64
		out    bool
	}{
		{
			// from woods:/ddn/home1/r2518/Programs/pbqff/runtests/sic/freqs
			infile: "full_tests/h2o.in",
			want:   "full_tests/h2o.9903",
			lline: []float64{
				-0.000045311426, -0.000027076533,
				0.000000000000, -0.000000002131,
			},
			eps: 6e-9,
		},
		{
			// from woods:/ddn/home1/r2518/chem/c3h2/sic/v5/freqs
			infile: "full_tests/c3h2.in",
			want:   "full_tests/c3h2.9903",
			lline: []float64{
				-0.000124209618, 0.000083980449, -0.000036821098,
				-0.000117696241, 0.000000000000, 0.000000000000,
				0.000000000000, 0.000000000000, 0.000000000000,
				-0.000000022736,
			},
			eps: 2e-7,
		},
		{
			// from woods:/ddn/home1/r2518/chem/c3h2/sic/010/v5/freqs
			infile: "full_tests/c3h2_010.in",
			want:   "full_tests/c3h2_010.9903",
			lline: []float64{
				-0.000124206320, 0.000083982488, -0.000036826840,
				-0.000117696763, 0.000000000000, 0.000000000000,
				0.000000000000, 0.000000000000, 0.000000000000,
				-0.000000022724,
			},
			eps: 5e-8,
		},
		{
			// from woods:/ddn/home1/r2516/chem/grads/hoof/qff/sic/freqs
			infile: "full_tests/hoof.in",
			want:   "full_tests/hoof.9903",
			lline: []float64{
				-0.000033920092, -0.000015460374, -0.000022566774,
				-0.000031468053, -0.000106254368, 0.000025443044,
				-0.000000003701,
			},
			eps: 2e-7,
		},
		{
			// from woods:/ddn/home1/r2518/chem/C4H3/0050/freqs/new_anpass
			infile: "full_tests/c4h3.in",
			want:   "full_tests/c4h3.9903",
			lline: []float64{
				-0.000065487334, -0.000918813588, 0.006806221930,
				-0.000082390448, 0.003530873987, -0.007497973052,
				0.000000000000, 0.000000000000, 0.000000000000,
				0.000000000000, 0.000000000000, 0.000000000000,
				0.000000000000, 0.000000000000, 0.000000000000,
				-0.000007996132,
			},
			lineps: 3e-8,
			eps:    3e-3,
		},
		{
			// from woods:/ddn/home7/r2831/chem/c5h2/cs/dzccr/005/freqs
			infile: "full_tests/c5h2.in",
			want:   "full_tests/c5h2.9903",
			lline: []float64{
				-0.000301517579, 0.000175887534, -0.000562965241,
				-0.000421688382, -0.000952657555, -0.001308750799,
				0.005821510949, -0.000282705915, -0.005921867331,
				-0.008693097178, -0.010256478892, 0.000000000000,
				0.000000000000, 0.000000000000, 0.000000000000,
				-0.000003173285,
			},
			lineps: 8e-9,
			eps:    2.5e-3,
		},
	}
	if testing.Short() {
		tests = tests[:len(tests)-2]
	}
	for _, test := range tests {
		w := io.Discard
		if test.out {
			w = os.Stdout
		}
		now := time.Now()
		fmt.Printf("starting %s\n", test.infile)
		disps, energies, exps, biases, _ := ReadInput(test.infile)
		disps, energies = Bias(disps, energies, biases)
		longLine, _, _ := Run(w, disps, energies, exps)
		if test.lineps == 0 {
			test.lineps = 1e-12
		}
		if !eql(longLine, test.lline, test.lineps) {
			t.Fatalf("got %v, wanted %v\n", longLine, test.lline)
		}
		disps, energies = Bias(disps, energies, longLine)
		_, got, _ := Run(w, disps, energies, exps)
		want := load9903(test.want)
		if !compFC(got, want, test.eps) {
			t.Errorf("FAIL %s\n", test.infile)
		}
		fmt.Printf("finished %s in %v\n", test.infile, time.Since(now))
	}
}
