ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
LDOBJS   =  Tabulated.o Range.o Random.o Coverage.o Pole.o
CFLAGS     = -O
#CFLAGS     = -pg -O

polelim:	polelim.o $(LDOBJS)
		g++ $(CFLAGS) -Wall $< $(LDOBJS) -ltclap -o $@
polelim.o:	polelim.cxx
		g++ $(CFLAGS) -Wall -c $<

polecov:	polecov.o $(LDOBJS)
		g++ $(CFLAGS) -Wall $< $(LDOBJS) -ltclap -o $@
polecov.o:	polecov.cxx
		g++ $(CFLAGS) -Wall -c $<


Pole.o:		Pole.cxx Pole.h Tabulated.h
		g++ $(CFLAGS) -Wall -c $<
Range.o:	Range.cxx Range.h
		g++ $(CFLAGS) -Wall -c $<
Coverage.o:	Coverage.cxx Coverage.h Range.h Tabulated.h Random.h Pole.h
		g++ $(CFLAGS) -Wall -c $<
Random.o:	Random.cxx Random.h
		g++ $(CFLAGS) -Wall -c $<
Tabulated.o:	Tabulated.cxx Tabulated.h
		g++ $(CFLAGS) -Wall -c $<
clean:
		rm -f *.o
