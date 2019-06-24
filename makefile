
SRCDIR=src
ODIR=build
ODIRCUDA=build_cuda
IDIR=include
LDIR_LOCAL=lib

PYDIR=gw_analysis_tools_py
PYSRC=mcmc_routines_ext.pyx waveform_generator_ext.pyx
PROJ_PYSRC=$(addprefix $(PYDIR)/src/,$(PYSRC))

PYLIB=mcmc_routines_ext.cpp waveform_generator_ext.cpp
PROJ_PYLIB=$(addprefix $(PYDIR)/,$(PYLIB))

LIBS=-ladolc -lgsl -lgslcblas -lfftw3 -llal -lcuda -lcudart 
LOCAL_LIB=libgwanalysistools.a
PROJ_LIB=$(addprefix $(LDIR_LOCAL)/,$(LOCAL_LIB))

TESTSRC=testing/test.cpp
TESTOBJ=testing/test.o
TESTDIR=bin
TEST=$(addprefix $(TESTDIR)/,exe.a)

TESTFISHERSRC=testing/fisher_comparison.cpp
TESTFISHEROBJ=testing/fisher_comparison.o
TESTFISHER=$(addprefix $(TESTDIR)/,exefisher.a)


#CFLAGS=-I$(IDIR) -I/opt/lalsuite/lalsimulation/src -I/opt/lalsuite/include -Wall -fPIC -g -O3 -std=c++11
CFLAGS=-I$(IDIR) -fopenmp -fPIC -g -O2 -std=c++11
#LFLAGS= -L/usr/local/cuda/lib64 -fopenmp 
LFLAGS= -fopenmp 
CFLAGSCUDA=-I$(IDIR) -shared -Xcompiler -fpic -O2 -std=c++11 
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(ODIR)/%,$(SOURCES:.$(SRCEXT)=.o))

SRCEXTCUDA := cu
SOURCESCUDA := $(shell find $(SRCDIR) -type f -name *.$(SRCEXTCUDA))
OBJECTSCUDA := $(patsubst $(SRCDIR)/%,$(ODIRCUDA)/%,$(SOURCESCUDA:.$(SRCEXTCUDA)=.o))

IEXT := h
DEPS:= $(shell find $(IDIR) -type f -name *.$(IEXT))

#CUDA specific header files -- not meant for external use
IEXTCUDA := hu
DEPSCUDA:= $(shell find $(IDIR) -type f -name *.$(IEXTCUDA))

#CC=g++-7
CC=g++
CCCUDA=nvcc
#CC=nvcc

.PHONY: all
all:  Doxyfile $(PROJ_LIB) $(PROJ_PYLIB)

$(ODIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(DEPS) 
	$(CC) $(CFLAGS) -c -o $@ $<

$(ODIRCUDA)/%.o : $(SRCDIR)/%.$(SRCEXTCUDA) $(DEPS) $(DEPSCUDA)
	$(CCCUDA) $(CFLAGSCUDA) -c -o $@ $<

$(OBJECTS): | $(ODIR)

$(OBJECTSCUDA): | $(ODIRCUDA)

$(ODIR):
	mkdir -p $(ODIR)

$(ODIRCUDA):
	mkdir -p $(ODIRCUDA)

$(TESTOBJ): $(TESTSRC)
	$(CC) $(CFLAGS) -c -o $@ $<

$(TESTFISHEROBJ): $(TESTFISHERSRC)
	$(CC) $(CFLAGS) -c -o $@ $<

$(PROJ_LIB) : $(OBJECTS) $(OBJECTSCUDA) | $(LDIR_LOCAL)
	ar rcs $(LDIR_LOCAL)/$(LOCAL_LIB) $(OBJECTS) $(OBJECTSCUDA)

$(LDIR_LOCAL):
	mkdir -p $(LDIR_LOCAL)

$(PROJ_PYLIB): $(PROJ_LIB) $(PROJ_PYSRC)
	make -C $(PYDIR) 
	
$(TESTFISHER) : $(OBJECTS) $(TESTFISHEROBJ) | $(TESTDIR)
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

$(TEST) : $(OBJECTS) $(OBJECTSCUDA) $(TESTOBJ) | $(TESTDIR)
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

$(TESTDIR):
	mkdir -p $(TESTDIR)

Doxyfile: $(OBJECTS) $(OBJECTSCUDA)
	doxygen Doxyfile
 
.PHONY: c
c: $(PROJ_LIB)

.PHONY: test
test: $(TEST) $(PROJ_LIB) $(PROJ_PYLIB) 

.PHONY: testfisher
testfisher: $(TESTFISHER) $(PROJ_LIB) 

.PHONY: testc
testc: $(TEST) $(PROJ_LIB)  

.PHONY: clean
clean:
	rm build/*.o 
	rm build_cuda/*.o 
.PHONY: remove
remove:
	rm build/*.o build_cuda/* $(TEST) lib/*.a
	make -C $(PYDIR) remove
