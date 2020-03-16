CC=gcc
CFLAGS=-I.
LIBS=-lm -fopenmp
OBJ=wray_sampling_analytic.o init.o integral.o test_function.o weight.o
ONAME=../binary/wray_sampling_analytic

%.o: %.c
	$(CC) -c -o $@ $<

$(ONAME): $(OBJ)
	$(CC) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o
