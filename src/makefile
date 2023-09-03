include ../make.sufrule
include ../make.inc

TOPSRCDIR = ..

ifdef USE_MPI
ifneq ($(USE_MPI), 0)
DEF_MPI = -DMPI
endif
endif

#MKL
mkl_include = -I$(MKL_DIR)/include
mkl_lib = -L$(MKL_DIR)/lib/intel64

#matlab -version 2021b
#If the version of matlab was smaller than 2021b, the path of libstdc++.so.6 wou-ld change
Matlab_stdc = $(MATLAB_DIR)/sys/os/glnxa64/orig/libstdc++.so.6
matlab_include = -I$(MATLAB_DIR)/extern/include -I$(MATLAB_DIR)/simulink/include
matlab_lib = -Wl,-rpath-link,$(MATLAB_DIR)/bin/glnxa64  -L$(MATLAB_DIR)/sys/os/glnxa64/orig -L$(MATLAB_DIR)/bin/glnxa64 -lmx -lmex -lmat -lm

#Local
local_lib = -L$(TOPSRCDIR)/lib -lzpares -lmatrix

FCCFLAG= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lmkl_blas95_ilp64

# INCLUDE = -I$(TOPSRCDIR)/include \
#                $(mkl_include) \
#                $(matlab_include) 
INCLUDE = -I$(TOPSRCDIR)/src\
        $(mkl_include) \
        $(matlab_include) 
LINK =  $(local_lib) \
        $(mkl_lib) $(FCCFLAG) \
        $(Matlab_stdc) \
        $(matlab_lib) 

BGK_OBJS = BGK_SS.o 

all : BGK_SS
BGK_SS : $(BGK_OBJS)
	ar r ../lib/libbgk_ss.a $(BGK_OBJS); \
	cp bgk_ss.mod ../include;\

BGK_SS.o      : BGK_SS.f90 \
                BGK_SS_Comm.f90 \
                BGK_SS_DM.f90 \
                BGK_SS_DBLA.f90 \
                BGK_SS_interface.f90

# BGK_SS_Comm.o : BGK_SS_Comm.f90 BGK_SS.o    
# BGK_SS_DM.o   : BGK_SS_DM.f90 BGK_SS.o
# BGK_SS_DBLA.o    : BGK_SS_DBLA.f90 BGK_SS.o
# BGK_SS_Interface : BGK_SS_Interface.o




clean :
	rm -f *.o *.mod ../include/bgk_ss.mod ../lib/libbgk_ss.a  
