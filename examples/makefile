TOPSRCDIR = ..
include $(TOPSRCDIR)/make.inc
include $(TOPSRCDIR)/make.sufrule

#MKL
mkl_include = -I$(MKL_DIR)/include
mkl_lib = -L$(MKL_DIR)/lib/intel64

#matlab -version 2021b
#If the version of matlab was smaller than 2021b, the path of libstdc++.so.6 wou-ld change 
Matlab_stdc = $(MATLAB_DIR)/sys/os/glnxa64/orig/libstdc++.so.6
matlab_include = -I$(MATLAB_DIR)/extern/include -I$(MATLAB_DIR)/simulink/include
matlab_lib = -Wl,-rpath-link,$(MATLAB_DIR)/bin/glnxa64  -L$(MATLAB_DIR)/sys/os/glnxa64/orig -L$(MATLAB_DIR)/bin/glnxa64 -lmx -lmex -lmat -lm

#Local
local_lib = -L$(TOPSRCDIR)/lib -lbgk_ss

FCCFLAG= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lmkl_blas95_ilp64 -liomp5

DEF_MPI = -DMPI
INCLUDE = -I$(TOPSRCDIR)/src\
        $(mkl_include) \
        $(matlab_include) 
LINK =  $(local_lib) \
        $(mkl_lib) $(FCCFLAG) \
        $(Matlab_stdc) \
        $(matlab_lib) 



all: ex1

ex1: ex1.o
	@echo  "Building" $@ "... "
	$(FC) $(FFLAG) $^ $(LINK) -o $@

.PHONY: clean
clean:
	rm -f *.o  
