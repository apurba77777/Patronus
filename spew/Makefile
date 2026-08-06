CC=gcc
CLIBS= -lgsl -lgslcblas -lm
CFLAGS=-I. -fopenmp
c = $(CC) $(CFLAGS) 
s = $(CC) -O0 -shared -fPIC $(CFLAGS) 

CUBEOBJ=cubestatfns.so cubestatfns.c

cubestat:
	$s -o $(CUBEOBJ) $(CLIBS)

clean:
	rm *.so
