CC=gcc
CFLAGS=-Wall -O3 -DHAVE_INLINE

all: HTightBinding.o BandEnergy.o

clean:
	rm *.o

HTightBinding.o:
	$(CC) $(CFLAGS) -c HTightBinding.c

BandEnergy.o:
	$(CC) $(CFLAGS) -c BandEnergy.c
