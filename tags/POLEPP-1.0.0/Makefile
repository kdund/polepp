DIRS	= ./src ./tools/polelim

.PHONY: default all clean_all doc clean_doc clean $(DIRS)

default: all

$(DIRS):
	@echo "**** $@ ****"
	@$(MAKE) -C $@ $(MAKE-TARGET)

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
