CC=gcc
CFLAGS=-Wall -O3 -DHAVE_INLINE
LDFLAGS=-lgsl -lgslcblas -lm
OBJFILES=bstrlib/bstrlib.o bstrlib/bstraux.o paths.o ParseSCF.o HTightBinding.o dos_util.o DosValues.o PartialDosValues.o PartialNumValues.o SpinOrbit.o BandEnergy.o ctetra/submesh.o ctetra/dos.o ctetra/partial.o ctetra/numstates.o ctetra/fermi.o ctetra/weights.o ctetra/sum.o ctetra/ecache.o ctetra/evcache.o ctetra/tetra.o

all: bstrlib.o bstraux.o paths.o ParseSCF.o HTightBinding.o DosValues.o PartialDosValues.o PartialNumValues.o SpinOrbit.o BandEnergy.o dos_util.o ParseSCF_test.out HTightBinding_test.out BandEnergy_test.out Anisotropy.out RunDosValues.out RunPartialDos.out RunPartialNum.out

clean:
	rm *.o *.out

bstrlib.o:
	$(CC) $(CFLAGS) -c bstrlib/bstrlib.c -o bstrlib/bstrlib.o

bstraux.o:
	$(CC) $(CFLAGS) -c bstrlib/bstraux.c -o bstrlib/bstraux.o

paths.o: bstrlib/bstrlib.o bstraux.o paths.c paths.h
	$(CC) $(CFLAGS) -c paths.c

ParseSCF.o: bstrlib/bstrlib.o bstraux.o ParseSCF.c ParseSCF.h
	$(CC) $(CFLAGS) -c ParseSCF.c

HTightBinding.o: bstrlib/bstrlib.o bstraux.o HTightBinding.c HTightBinding.h
	$(CC) $(CFLAGS) -c HTightBinding.c

SpinOrbit.o: SpinOrbit.c SpinOrbit.h HTightBinding.o
	$(CC) $(CFLAGS) -c SpinOrbit.c

BandEnergy.o: HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o BandEnergy.c BandEnergy.h
	$(CC) $(CFLAGS) -c BandEnergy.c

dos_util.o: dos_util.c dos_util.h
	$(CC) $(CFLAGS) -c dos_util.c

DosValues.o: HTightBinding.o DosValues.c DosValues.h
	$(CC) $(CFLAGS) -c DosValues.c

PartialDosValues.o: HTightBinding.o PartialDosValues.c PartialDosValues.h
	$(CC) $(CFLAGS) -c PartialDosValues.c

PartialNumValues.o: HTightBinding.o PartialDosValues.o PartialNumValues.c PartialNumValues.h
	$(CC) $(CFLAGS) -c PartialNumValues.c

ParseSCF_test.out: ParseSCF.o bstrlib/bstrlib.o bstrlib/bstraux.o ParseSCF_test.c
	$(CC) $(CFLAGS) -c ParseSCF_test.c -o ParseSCF_test.o
	$(CC) ParseSCF_test.o -o ParseSCF_test.out $(OBJFILES) $(LDFLAGS)

HTightBinding_test.out: HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o HTightBinding_test.c
	$(CC) $(CFLAGS) -c HTightBinding_test.c -o HTightBinding_test.o
	$(CC) HTightBinding_test.o -o HTightBinding_test.out $(OBJFILES) $(LDFLAGS)

BandEnergy_test.out: BandEnergy.o HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o BandEnergy_test.c
	$(CC) $(CFLAGS) -c BandEnergy_test.c -o BandEnergy_test.o
	$(CC) BandEnergy_test.o -o BandEnergy_test.out $(OBJFILES) $(LDFLAGS)

Anisotropy.out: paths.o ParseSCF.o SpinOrbit.o BandEnergy.o HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o Anisotropy.c
	$(CC) $(CFLAGS) -c Anisotropy.c -o Anisotropy.o
	$(CC) Anisotropy.o -o Anisotropy.out $(OBJFILES) $(LDFLAGS)

RunDosValues.out: HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o DosValues.o dos_util.o RunDosValues.c
	$(CC) $(CFLAGS) -c RunDosValues.c -o RunDosValues.o
	$(CC) RunDosValues.o -o RunDosValues.out $(OBJFILES) $(LDFLAGS)

RunPartialDos.out: HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o PartialDosValues.o dos_util.o RunPartialDos.c
	$(CC) $(CFLAGS) -c RunPartialDos.c -o RunPartialDos.o
	$(CC) RunPartialDos.o -o RunPartialDos.out $(OBJFILES) $(LDFLAGS)

RunPartialNum.out: HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o PartialNumValues.o dos_util.o RunPartialNum.c
	$(CC) $(CFLAGS) -c RunPartialNum.c -o RunPartialNum.o
	$(CC) RunPartialNum.o -o RunPartialNum.out $(OBJFILES) $(LDFLAGS)
