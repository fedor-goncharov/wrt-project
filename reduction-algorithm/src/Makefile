CC=gcc
LIBS=-lm -fopenmp -lgsl -lgslcblas
OPENMP_LIBS=fopenmp
OBJ=reduction.o init.o slicing.o
ONAME=../binary/radon_reduction

%.o: %.c
	$(CC) -c -o $@ $< $(LIBS)

$(ONAME): $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
