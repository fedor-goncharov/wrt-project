CC=gcc
CFLAGS=-I.
LIBS=-lm -fopenmp -lgsl -lgslcblas
DEPS=init.h integral.h interp.h
OBJ=radon_grid.o init.o integral.o interp.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS) 

../binary/radon_grid: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
