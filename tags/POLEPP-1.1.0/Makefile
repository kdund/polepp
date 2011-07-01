DIRS	= ./src ./tools/polelim ./tools/polecov

.PHONY: default all clean_all doc clean_doc clean help $(DIRS)

default: help

$(DIRS):
	@echo "**** $@ ****"
	@$(MAKE) -C $@ $(MAKE-TARGET)

help:
	@echo ""
	@echo "make all      make all binaries"
	@echo "make clean    cleanup"
	@echo "make debug    make with added debug info"
	@echo "make doc      create documentation"
	@echo ""
	@echo "User settable variables:"
	@echo "  GSL_INCL   : GSL root directory   (def: /usr/include)"
	@echo "  GSL_LIB    : GSL library path     (def: /usr/lib)"
	@echo "  TCLAP_ROOT : TCLAP root directory (def: /usr/include)"
	@echo ""
	@echo "Examples:"
	@echo "  make all                             - make using all default values"
	@echo "  make GSL_LIB=/usr/lib64 all          - change GSL path"
	@echo ""

debug:
	USE_DEBUG = 1
	@$(MAKE) MAKE-TARGET=all $(DIRS)
all:
	@$(MAKE) MAKE-TARGET=all $(DIRS)

clean_all:
	@$(MAKE) MAKE-TARGET=clean $(DIRS)
	@$(RM) -r bin lib data

doc:
	@echo "  Doxygenating"
	@(export DIRS="$(DIRS)" ; doxygen)

clean_doc:
	@$(RM) -r doc

clean: clean_all clean_doc
