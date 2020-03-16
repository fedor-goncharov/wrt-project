CC=gcc
CFLAGS=-I.
LIBS=-lm -fopenmp -lgsl -lgslcblas
DEPS=init.h integral.h test_function.h
TFUNC=test_function.c
OBJ=radon3d_analytic.o init.o integral.o ${TFUNC%.*}.o
ODIR=../binary
ONAME=(output name of binary here)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

$(ODIR)/$(ONAME): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
