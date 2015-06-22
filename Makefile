CC=gcc
CFLAGS=-Wall -O3 -DHAVE_INLINE
LDFLAGS=-lgsl -lgslcblas -lm
OBJFILES=bstrlib/bstrlib.o bstrlib/bstraux.o HTightBinding.o BandEnergy.o ctetra/dos.o ctetra/numstates.o ctetra/weights.o ctetra/sum.o

all: bstrlib.o bstraux.o HTightBinding.o BandEnergy.o HTightBinding_test.out BandEnergy_test.out

clean:
	rm *.o *.out

bstrlib.o:
	$(CC) $(CFLAGS) -c bstrlib/bstrlib.c -o bstrlib/bstrlib.o

bstraux.o:
	$(CC) $(CFLAGS) -c bstrlib/bstraux.c -o bstrlib/bstraux.o

HTightBinding.o: bstrlib/bstrlib.o bstraux.o HTightBinding.c HTightBinding.h
	$(CC) $(CFLAGS) -c HTightBinding.c

BandEnergy.o: HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o BandEnergy.c BandEnergy.h
	$(CC) $(CFLAGS) -c BandEnergy.c

HTightBinding_test.out: HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o HTightBinding_test.c
	$(CC) $(CFLAGS) -c HTightBinding_test.c -o HTightBinding_test.o
	$(CC) HTightBinding_test.o -o HTightBinding_test.out $(OBJFILES) $(LDFLAGS)

BandEnergy_test.out: BandEnergy.o HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o BandEnergy_test.c
	$(CC) $(CFLAGS) -c BandEnergy_test.c -o BandEnergy_test.o
	$(CC) BandEnergy_test.o -o BandEnergy_test.out $(OBJFILES) $(LDFLAGS)
