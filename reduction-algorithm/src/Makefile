CC=gcc
LIBS=-lm -fopenmp -lgsl -lgslcblas
OPENMP_LIBS=fopenmp
OBJ=reduction.o init.o slicing.o
ONAME=../binary/radon_reduction

%.o: %.c
	$(CC) -c -o $@ $< -O3 $(LIBS)

$(ONAME): $(OBJ)
	$(CC) -o $@ $^ -O3 $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
