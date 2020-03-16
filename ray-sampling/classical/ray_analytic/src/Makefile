CC=gcc
LIBS=-lm -fopenmp
OBJ=ray_analytic.o init.o integral.o test_function.o
ONAME=ray_analytic

%.o: %.c
	$(CC) -c -o $@ $< $(LIBS)

$(ONAME): $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
