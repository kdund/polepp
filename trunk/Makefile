 ######################################################################
 # Project: Pole++                                                    #
 ###################################################################### 
.SUFFIXES: .cxx .h .o
MAKEFLAGS = --no-print-directory -r -s
#MAKEFLAGS = --warn-undefined-variables --debug

include Makefile.arch

# Internal configuration
PACKAGE=Pole++
LIBVERS=.1
LD_LIBRARY_PATH:=.:$(ROOTSYS)/lib:$(LD_LIBRARY_PATH)
OBJDIR=/opt/scratch/obj
DEPDIR=$(OBJDIR)/dep
VPATH= $(OBJDIR)
INCLUDES += -I./  
ROOTSYS  ?= ERROR_RootSysIsNotDefined

TOOLLIST = polelim.cxx polecov.cxx polebelt.cxx exptest.cxx plotexp.cxx poleconst.cxx polepow.cxx
LIBLIST  = Pole.cxx  Coverage.cxx Random.cxx Pdf.cxx Combination.cxx Tools.cxx
ARGLIST  = argsPole.cxx argsCoverage.cxx
SKIPLIST = Combine.cxx Power.cxx Tabulated.cxx pdftst.cxx polecomb.cxx plotexp.cxx poleconst.cxx polepow.cxx
SKIPLIBLIST  = $(SKIPLIST) $(TOOLLIST) $(ARGLIST)
SKIPTOOLLIST = $(SKIPLIST) $(LIBLIST)
SKIPHLIST = Tabulated.h
LIBFILE      = lib$(PACKAGE).a${LIBVERS}
SHLIBFILE    = lib$(PACKAGE).so${LIBVERS}

UNAME = $(shell uname)


default: shlib


# List of all include files
HLIST   = $(filter-out $(SKIPHLIST),$(wildcard *.h))

# source for library
LIBCPPLIST   =  $(filter-out $(SKIPLIBLIST),$(LIBLIST))

# source for tools
TOOLCPPLIST  =  $(filter-out $(SKIPTOOLLIST),$(TOOLLIST))

# List of all object files to build
LIBOLIST =$(patsubst %.cxx,%.o,$(LIBCPPLIST))
TOOLOLIST=$(patsubst %.cxx,%.o,$(TOOLCPPLIST))
TOOLELIST=$(patsubst %.cxx,%,$(TOOLCPPLIST))

# list of scripts
SCRIPTLIST = log2confbelt.csh log2construct.csh plotconfbelt.C plotconst.C

# list of files for package
PACKAGELIST = $(LIBCPPLIST) $(TOOLCPPLIST) $(HLIST) $(ARGLIST) $(SCRIPTLIST) \
		Makefile Makefile.arch release.notes Doxyfile README

# rules to compile tools

tools: $(TOOLELIST)

polelim:	polelim.o argsPole.o
		$(CXX) $(CXXFLAGS) -Wall $(addprefix $(OBJDIR)/,$(notdir $^)) $(SHLIBFILE)  -o $@
polecov:	polecov.o argsCoverage.o
		$(CXX) $(CXXFLAGS) -Wall $(addprefix $(OBJDIR)/,$(notdir $^)) $(SHLIBFILE)  -o $@
polebelt:	polebelt.o argsPole.o
		$(CXX) $(CXXFLAGS) -Wall $(addprefix $(OBJDIR)/,$(notdir $^)) $(SHLIBFILE)  -o $@
poleconst:	poleconst.o argsPole.o
		$(CXX) $(CXXFLAGS) -Wall $(addprefix $(OBJDIR)/,$(notdir $^)) $(SHLIBFILE)  -o $@
exptest:	exptest.o
		$(CXX) $(CXXFLAGS) -Wall $(addprefix $(OBJDIR)/,$(notdir $^)) $(SHLIBFILE)  -o $@

# Explicit rule for args*.o
argsPole.o:     argsPole.cxx
		@echo -n "Compiling $< ... "
		@mkdir -p $(OBJDIR)
		$(CXX) $(INCLUDES) -fPIC -Wall -g -c $< -o $(OBJDIR)/$(notdir $@)
		@echo "Done"

argsCoverage.o: argsCoverage.cxx
		@echo -n "Compiling $< ... "
		@mkdir -p $(OBJDIR)
		$(CXX) $(INCLUDES) -fPIC -Wall -g -c $< -o $(OBJDIR)/$(notdir $@)
		@echo "Done"

# Implicit rule to compile all classes
%.o : %.cxx 
	@echo -n "Compiling $< ... "
	@mkdir -p $(OBJDIR)
	$(CXX) $(INCLUDES) $(CXXFLAGS) -g -c $< -o $(OBJDIR)/$(notdir $@)
	@echo "Done"

##############################
# The dependencies section   
# - the purpose of the .d files is to keep track of the
#   header file dependence
# - this can be achieved using the makedepend command 
##############################
# .d tries to pre-process .cc
-include $(foreach var,$(LIBCPPLIST:.$(SrcSuf)=.d),$(DEPDIR)/$(var))
-include $(foreach var,$(TOOLCPPLIST:.$(SrcSuf)=.d),$(DEPDIR)/$(var))
#-include $(foreach var,$(CPPLIST:.$(SrcSuf)=.d),$(DEPDIR)/$(var)) /dev/null


$(DEPDIR)/%.d: %.$(SrcSuf)
	@mkdir -p $(DEPDIR)
	if test -f $< ; then \
		echo -n "Building $(@F) ... "; \
		$(SHELL) -ec '$(CPP) -MM $(INCLUDES) $(CXXFLAGS) $< | sed '\''/Cstd\/rw/d'\'' > $@'; \
		echo "Done"; \
	fi

# Rule to combine objects into a library
$(LIBFILE): $(LIBOLIST)
	@echo -n "Making static library $(LIBFILE) ... "
	@rm -f $(LIBFILE)
	@ar q $(LIBFILE) $(addprefix $(OBJDIR)/,$(LIBOLIST))
	@ranlib $(LIBFILE)
	@echo "Done"

# Rule to combine objects into a unix shared library
$(SHLIBFILE): $(LIBOLIST)
	@echo -n "Building shared library $(SHLIBFILE) ... "
	@rm -f $(SHLIBFILE)
	$(LD) $(SOFLAGS) $(addprefix $(OBJDIR)/,$(LIBOLIST)) -o $(SHLIBFILE)
	@echo "Done"

# Rule to make package
package:
	@echo -n "Creating package... "
	@tar -czf polelib.tgz $(PACKAGELIST)
	@echo "Done"

# Useful build targets
lib: $(LIBFILE) 
shlib: $(SHLIBFILE)
all: $(SHLIBFILE) $(TOOLELIST)
clean:
	rm -f $(SHLIBFILE)
	rm -f $(OBJDIR)/*.o
	rm -f $(DEPDIR)/*.d
	rm -f $(LIBFILE)
	rm -f $(SHLIBFILE)

distclean:
	rm -rf obj 
	rm -f *~
	rm -f $(SHLIBFILE)
	rm -f $(LIBFILE)
	rm -f $(SHLIBFILE)

.PHONY : shlib lib default clean polelim


