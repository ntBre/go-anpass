package main

import (
	"flag"
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/ntBre/anpass"
)

var (
	debug = flag.Bool("debug", false,
		"whether to print debugging information")
	quiet = flag.Bool("q", false,
		"print nothing, don't even make output file")
	once = flag.Bool("once", false,
		"only run one pass, don't refit to stationary point")
)

func main() {
	flag.Parse()
	args := flag.Args()
	var infile, outfile, infile2, outfile2 string
	switch len(args) {
	case 1:
		infile = args[0]
		outfile = strings.Replace(infile, "in", "out", -1)
	case 2:
		infile = args[0]
		outfile = args[1]
	default:
		panic("not enough args")
	}
	var out io.Writer
	if !*quiet {
		f, err := os.Create(outfile)
		defer f.Close()
		if err != nil {
			panic(err)
		}
		out = f
	} else {
		out = io.Discard
	}
	disps, energies, exps, biases, stationary := anpass.ReadInput(infile)
	anpass.PrintBias(out, biases)
	disps, energies = anpass.Bias(disps, energies, biases)
	dir := filepath.Dir(infile)
	longLine, _, stationary := anpass.Run(out, dir, disps, energies, exps)
	// pass the longline and do anpass2 if the first run wasn't on
	// a stationary point
	if !*once && !stationary {
		infile2 = "anpass2.in"
		anpass.CopyAnpass(infile, infile2, longLine)
		outfile2 = strings.Replace(infile2, "in", "out", -1)
		if !*quiet {
			f, err := os.Create(outfile2)
			defer f.Close()
			if err != nil {
				panic(err)
			}
			out = f
		} else {
			out = io.Discard
		}
		anpass.PrintBias(out, longLine)
		disps, energies = anpass.Bias(disps, energies, longLine)
		anpass.Run(out, dir, disps, energies, exps)
	}
}
