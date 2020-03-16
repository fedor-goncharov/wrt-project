CC=gcc
CFLAGS=-I.
LIBS=-lm -fopenmp
DEPS=init.h integral.h interp.h
OBJ=ray_grid.o init.o integral.o interp.o

%.o: %.c $(DEPS)
	$(CC) -c $< -o $@ $(LIBS)

ray_grid: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
