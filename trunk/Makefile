ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
LDOBJS   =  Pdf.o Range.o Random.o Coverage.o Pole.o
###CFLAGS     = -O
# Use -g to include debug info
# Use -pg for profiling info
CFLAGS     = -g -O
SRCTOOLS     = polelim.cxx polecov.cxx exptest.cxx plotexp.cxx
SRCLIB       = Pole.cxx Pole.h Coverage.cxx Coverage.h Random.cxx Random.h Range.cxx Range.h Pdf.cxx Pdf.h

all:		polelim polecov exptest plotexp

polelim:	polelim.o $(LDOBJS)
		g++ $(CFLAGS) -Wall $< $(LDOBJS) -ltclap -o $@
polelim.o:	polelim.cxx
		g++ $(CFLAGS) -Wall -c $<

polecov:	polecov.o $(LDOBJS)
		g++ $(CFLAGS) -Wall $< $(LDOBJS) -ltclap -o $@
polecov.o:	polecov.cxx
		g++ $(CFLAGS) -Wall -c $<

exptest:	exptest.o $(LDOBJS)
		g++ $(CFLAGS) -Wall $< $(LDOBJS) -ltclap -o $@
exptest.o:	exptest.cxx
		g++ $(CFLAGS) -Wall -c $<

plotexp:	plotexp.o
		g++ $(CFLAGS) -Wall $< $(ROOTLIBS) -ltclap -o $@
plotexp.o:	plotexp.cxx
		g++ $(CFLAGS) $(ROOTCFLAGS) -Wall -c $<

obstest:	obstest.o $(LDOBJS) Observable.o
		g++ $(CFLAGS) -Wall $< $(LDOBJS) Observable.o -ltclap -o $@
obstest.o:	obstest.cxx
		g++ $(CFLAGS) -Wall -c $<

Pole.o:		Pole.cxx Pole.h Pdf.h Range.h
		g++ $(CFLAGS) -Wall -c $<
Range.o:	Range.cxx Range.h
		g++ $(CFLAGS) -Wall -c $<
Coverage.o:	Coverage.cxx Coverage.h Range.h Random.h Pole.h Pdf.h
		g++ $(CFLAGS) -Wall -c $<
Random.o:	Random.cxx Random.h
		g++ $(CFLAGS) -Wall -c $<
Pdf.o:		Pdf.cxx Pdf.h
		g++ $(CFLAGS) -Wall -c $<
Observable.o:	Observable.cxx Observable.h Random.h
		g++ $(CFLAGS) -Wall -c $<

package:	$(SRCLIB) $(SRCTOOLS) Makefile
		tar -czf polelib.tgz $(SRCLIB) $(SRCTOOLS) Makefile

clean:
		rm -f *.o
