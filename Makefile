CC=gcc
CFLAGS=-Wall -O3 -DHAVE_INLINE
LDFLAGS=-lgsl -lgslcblas -lm
OBJFILES=bstrlib/bstrlib.o bstrlib/bstraux.o ParseSCF.o HTightBinding.o SpinOrbit.o BandEnergy.o ctetra/submesh.o ctetra/dos.o ctetra/numstates.o ctetra/fermi.o ctetra/weights.o ctetra/sum.o ctetra/ecache.o

all: bstrlib.o bstraux.o ParseSCF.o HTightBinding.o SpinOrbit.o BandEnergy.o ParseSCF_test.out HTightBinding_test.out BandEnergy_test.out

clean:
	rm *.o *.out

bstrlib.o:
	$(CC) $(CFLAGS) -c bstrlib/bstrlib.c -o bstrlib/bstrlib.o

bstraux.o:
	$(CC) $(CFLAGS) -c bstrlib/bstraux.c -o bstrlib/bstraux.o

ParseSCF.o: bstrlib/bstrlib.o bstraux.o
	$(CC) $(CFLAGS) -c ParseSCF.c

HTightBinding.o: bstrlib/bstrlib.o bstraux.o HTightBinding.c HTightBinding.h
	$(CC) $(CFLAGS) -c HTightBinding.c

SpinOrbit.o: SpinOrbit.c SpinOrbit.h HTightBinding.o
	$(CC) $(CFLAGS) -c SpinOrbit.c

BandEnergy.o: HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o BandEnergy.c BandEnergy.h
	$(CC) $(CFLAGS) -c BandEnergy.c

ParseSCF_test.out: ParseSCF.o bstrlib/bstrlib.o bstrlib/bstraux.o ParseSCF_test.c
	$(CC) $(CFLAGS) -c ParseSCF_test.c -o ParseSCF_test.o
	$(CC) ParseSCF_test.o -o ParseSCF_test.out $(OBJFILES) $(LDFLAGS)

HTightBinding_test.out: HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o HTightBinding_test.c
	$(CC) $(CFLAGS) -c HTightBinding_test.c -o HTightBinding_test.o
	$(CC) HTightBinding_test.o -o HTightBinding_test.out $(OBJFILES) $(LDFLAGS)

BandEnergy_test.out: BandEnergy.o HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o BandEnergy_test.c
	$(CC) $(CFLAGS) -c BandEnergy_test.c -o BandEnergy_test.o
	$(CC) BandEnergy_test.o -o BandEnergy_test.out $(OBJFILES) $(LDFLAGS)
