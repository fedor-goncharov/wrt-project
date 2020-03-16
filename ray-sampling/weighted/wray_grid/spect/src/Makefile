CC=gcc
CFLAGS=-I.
LIBS=-lm -fopenmp
DEPS=init.h integral.h interp.h
OBJ=ray3d_grid.o init.o integral.o interp.o
ODIR=../binary

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

$(ODIR)/ray3d_grid: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
