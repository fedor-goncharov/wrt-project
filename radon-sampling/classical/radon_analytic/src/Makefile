CC=gcc
LIBS=-lm -fopenmp -lgsl -lgslcblas
OBJ=radon_analytic.o init.o integral.o test_function.o

%.o: %.c
	$(CC) -c -o $@ $< $(LIBS)

radon_analytic: $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
