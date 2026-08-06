CC=gcc
CLIBS= -lgsl -lgslcblas -lm
CFLAGS=-I. -fopenmp
c = $(CC) $(CFLAGS) 
s = $(CC) -O0 -shared -fPIC $(CFLAGS) 

STATOBJ=cubestatfns.so cubestatfns.c

SRCHOBJ=cubesrchfns.so cubesrchfns.c

cubestat:
	$s -o $(STATOBJ) $(CLIBS)

cubesrch:
	$s -o $(SRCHOBJ) $(CLIBS)

clean:
	rm *.so
