ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
LDOBJS   =  Pdf.o Range.o Random.o Coverage.o Pole.o
###CFLAGS     = -O
# Use -g to include debug info
# Use -pg for profiling info
CFLAGS     = -g -O
SRCTOOLS     = polelim.cxx polecov.cxx exptest.cxx plotexp.cxx estbelt.cxx polebelt.cxx poleconst.cxx
INCFILES     = Pole.h Coverage.h Random.h Range.h Pdf.h BeltEstimator.h
SRCFILES     = Pole.cxx Coverage.cxx Random.cxx Range.cxx Pdf.cxx

all:		polelim polebelt polecov exptest plotexp estbelt poleconst

estbelt:	estbelt.o
		g++ $(CFLAGS) -Wall $< -o $@
estbelt.o:	estbelt.cxx
		g++ $(CFLAGS) -Wall -c $<

poleconst:	poleconst.o $(LDOBJS)
		g++ $(CFLAGS) -Wall $< $(LDOBJS)  -o $@
poleconst.o:	poleconst.cxx
		g++ $(CFLAGS) -Wall -c $<

polelim:	polelim.o $(LDOBJS)
		g++ $(CFLAGS) -Wall $< $(LDOBJS)  -o $@
polelim.o:	polelim.cxx
		g++ $(CFLAGS) -Wall -c $<

polebelt:	polebelt.o $(LDOBJS)
		g++ $(CFLAGS) -Wall $< $(LDOBJS)  -o $@
polebelt.o:	polebelt.cxx
		g++ $(CFLAGS) -Wall -c $<

polecov:	polecov.o $(LDOBJS)
		g++ $(CFLAGS) -Wall $< $(LDOBJS)  -o $@
polecov.o:	polecov.cxx
		g++ $(CFLAGS) -Wall -c $<

exptest:	exptest.o $(LDOBJS)
		g++ $(CFLAGS) -Wall $< $(LDOBJS)  -o $@
exptest.o:	exptest.cxx
		g++ $(CFLAGS) -Wall -c $<

plotexp:	plotexp.o
		g++ $(CFLAGS) -Wall $< $(ROOTLIBS)  -o $@
plotexp.o:	plotexp.cxx
		g++ $(CFLAGS) $(ROOTCFLAGS) -Wall -c $<

obstest:	obstest.o $(LDOBJS) Observable.o
		g++ $(CFLAGS) -Wall $< $(LDOBJS) Observable.o  -o $@
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
		tar -czf polelib.tgz $(SRCFILES) $(INCFILES) $(SRCTOOLS) Makefile

clean:
		rm -f *.o
