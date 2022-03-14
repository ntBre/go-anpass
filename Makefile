TFLAGS = -v -failfast -parallel 1

SHORT = 0
ifeq ($(SHORT),0)
TFLAGS += -short
endif

test:
	go test . $(TFLAGS)

cpu.prof: *.go
	go test -cpuprofile cpu.prof . -run Full -v

cover:
	go test . -v -short -coverprofile=/tmp/anpass.out
	go tool cover -html /tmp/anpass.out

debug:
	go build -gcflags "-N -l" .

anpass: *.go
	go build .

bench:
	go test . -bench=. -run BenchmarkFit

newton.prof: *.go
	go test . -bench Newton -run BenchmarkNewton -cpuprofile $@

fit.prof: *.go
	go test . -bench Fit -run BenchmarkFit -cpuprofile $@
