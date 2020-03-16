CC=gcc 
LIBS=-lm -fopenmp -w
OBJ=spect_weight_expansions.o init.o integral.o interp.o aggregate.o
ODIR=../binary

%.o: %.c
	$(CC) -c $< -o $@ $(LIBS)

$(ODIR)/spect_weight_expansions: $(OBJ)
	$(CC) $^ -o $@ $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
