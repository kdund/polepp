ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
LDOBJS     =  Pdf.o Random.o Pole.o Combination.o Tools.o Coverage.o
LDVERSION  = .1
###LDOBJS     =  Pdf.o Random.o Coverage.o Pole.o Combination.o Combine.o
#### Power.o
LDSMOBJS     =  Pdf.o Random.o
###LDFLAGS      = -fPIC -shared
LDFLAGS      = -shared
###CFLAGS     = -O
# Use -g to include debug info
# Use -pg for profiling info
# Use -O3 for optimizing
CFLAGS       = -g -O0
CFLAGSARG    = -g -O0
SRCTOOLS     = polelim.cxx polecov.cxx argsPole.cxx argsCoverage.cxx
SRCEXTRAS    = exptest.cxx plotexp.cxx estbelt.cxx polebelt.cxx poleconst.cxx polepow.cxx
INCFILES     = Pole.h Coverage.h Random.h Range.h Pdf.h BeltEstimator.h Combination.h Measurement.h \
	       Power.h Observable.h Combine.h Tools.h
SRCFILES     = Pole.cxx Coverage.cxx Random.cxx Pdf.cxx Combination.cxx Power.cxx Combine.cxx Tools.cxx
LDFILE	     = libPole++.so
LDSMFILE     = libPoleSM++.so
LIBS         = -L./ -lPole++
LIBSSM       = -L./ -lPoleSM++

all:		$(LDFILE) polelim polecov
extras:         polebelt exptest plotexp estbelt poleconst polecomb polepow
small:          $(LDSMFILE)

$(LDFILE):	$(LDOBJS)
		g++ $(LDFLAGS) $(LDOBJS) -o $(LDFILE)

$(LDSMFILE):	$(LDSMOBJS)
		g++ $(LDFLAGS) $(LDSMOBJS) -o $(LDSMFILE)

argsPole.o:	argsPole.cxx
		g++ $(CFLAGSARG) -Wall -c $<
argsCoverage.o:	argsCoverage.cxx
		g++ $(CFLAGSARG) -Wall -c $<

polepow:	polepow.o
		g++ $(CFLAGS) -Wall $< $(LIBS)  -o $@
polepow.o:	polepow.cxx
		g++ $(CFLAGS) -Wall -c $<

polecomb:       polecomb.o
		g++ $(CFLAGS) -Wall $< $(LIBS)  -o $@
polecomb.o:     polecomb.cxx
		g++ $(CFLAGS) -Wall -c $<

estbelt:	estbelt.o
		g++ $(CFLAGS) -Wall $< $(LIBS)  -o $@
estbelt.o:	estbelt.cxx
		g++ $(CFLAGS) -Wall -c $<

poleconst:	poleconst.o
		g++ $(CFLAGS) -Wall $< $(LIBS)  -o $@
poleconst.o:	poleconst.cxx
		g++ $(CFLAGS) -Wall -c $<

polelim:	polelim.o argsPole.o
		g++ $(CFLAGS) -Wall $^ $(LIBS)  -o $@
polelim.o:	polelim.cxx
		g++ $(CFLAGS) -Wall -c $<

polebelt:	polebelt.o argsPole.o
		g++ $(CFLAGS) -Wall $^ $(LIBS)  -o $@
polebelt.o:	polebelt.cxx
		g++ $(CFLAGS) -Wall -c $<

polecov:	polecov.o argsCoverage.o
		g++ $(CFLAGS) -Wall $^ $(LIBS)  -o $@
polecov.o:	polecov.cxx
		g++ $(CFLAGS) -Wall -c $<

exptest:	exptest.o
		g++ $(CFLAGS) -Wall $< $(LIBS)  -o $@
exptest.o:	exptest.cxx
		g++ $(CFLAGS) -Wall -c $<

plotexp:	plotexp.o
		g++ $(CFLAGS) -Wall $< $(ROOTLIBS)  -o $@
plotexp.o:	plotexp.cxx
		g++ $(CFLAGS) $(ROOTCFLAGS) -Wall -c $<

obstest:	obstest.o
		g++ $(CFLAGS) -Wall $< $(LIBS)  -o $@
obstest.o:	obstest.cxx
		g++ $(CFLAGS) -Wall -c $<

pdftst:		pdftst.o
		g++ $(CFLAGS) -Wall $^ $(LIBS)  -o $@
pdftst.o:	pdftst.cxx
		g++ $(CFLAGS) -Wall -c $<

pdftstold:	pdftstold.o
		g++ $(CFLAGS) -Wall $< $(LIBS)  -o $@
pdftstold.o:	pdftstold.cxx
		g++ $(CFLAGS) -Wall -c $<

meastst:	meastst.o
		g++ $(CFLAGS) -Wall $< $(LIBSSM)  -o $@
meastst.o:	meastst.cxx
		g++ $(CFLAGS) -Wall -c $<

Combine.o:      Combine.cxx Combine.h
		g++ $(CFLAGS) -Wall -c $<

Pole.o:		Pole.cxx Pole.h Pdf.h Range.h
		g++ $(CFLAGS) -Wall -c $<
Coverage.o:	Coverage.cxx Coverage.h Range.h Random.h Pole.h Pdf.h
		g++ $(CFLAGS) -Wall -c $<
Random.o:	Random.cxx Random.h
		g++ $(CFLAGS) -Wall -c $<
Pdf.o:		Pdf.cxx Pdf.h
		g++ $(CFLAGS) -Wall -c $<
Combination.o:	Combination.cxx Combination.h
		g++ $(CFLAGS) -Wall -c $<
Power.o:	Power.cxx Power.h
		g++ $(CFLAGS) -Wall -c $<
Tools.o:	Tools.cxx Tools.h
		g++ $(CFLAGS) -Wall -c $<

package:	$(SRCLIB) $(SRCTOOLS) Makefile
		tar -czf polelib.tgz $(SRCFILES) $(INCFILES) $(SRCTOOLS) $(SRCEXTRAS) release.notes Makefile Doxyfile README

clean:
		rm -f *.o
