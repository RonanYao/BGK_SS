TOPSRCDIR =.
include $(TOPSRCDIR)/make.inc
include $(TOPSRCDIR)/make.sufrule


all: main

main:ECHO
	(cd src && $(MAKE))
	(cd examples && $(MAKE))

ECHO:
	@echo begin compile


.PHONY: clean
clean: 
	(cd src && $(MAKE) clean)
	(cd examples && $(MAKE) clean)
	
