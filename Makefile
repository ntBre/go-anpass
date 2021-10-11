TFLAGS = -v -failfast -parallel 1

ifeq ($(SHORT),1)
TFLAGS += -short
endif

test:
	go test . $(TFLAGS)

cpu.prof: *.go
	go test -cpuprofile cpu.prof . -run Full -v

debug:
	go build -gcflags "-N -l" .

anpass: *.go
	go build .

bench:
	go test . -bench=. -run BenchmarkFit
