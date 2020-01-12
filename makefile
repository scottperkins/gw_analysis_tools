SRCDIR=src
ODIR=build
ODIRCUDA=build_cuda
IDIR=include/gwat/
LDIR_LOCAL=lib

PYDIR=gwatpy
PYSRC=mcmc_routines_ext.pyx waveform_generator_ext.pyx
PROJ_PYSRC=$(addprefix $(PYDIR)/src/,$(PYSRC))

PYLIB=mcmc_routines_ext.cpp waveform_generator_ext.cpp
PROJ_PYLIB=$(addprefix $(PYDIR)/,$(PYLIB))

LOCAL_LIB=libgwat.a
LOCAL_SHARED_LIB=libgwat.so
PROJ_LIB=$(addprefix $(LDIR_LOCAL)/,$(LOCAL_LIB))
PROJ_SHARED_LIB=$(addprefix $(LDIR_LOCAL)/,$(LOCAL_SHARED_LIB))
EXEDIR:=bin
MCMC_TOOL:=$(EXEDIR)/mcmc_gw_tool

TESTSRC=testing/test.cpp
TESTOBJ=testing/test.o
TESTDIR=bin
TEST=$(addprefix $(TESTDIR)/,exe.a)

TESTFISHERSRC=testing/fisher_comparison.cpp
TESTFISHEROBJ=testing/fisher_comparison.o
TESTFISHER=$(addprefix $(TESTDIR)/,exefisher.a)

CONFIGFILE:=include/gwat/GWATConfig.h

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

############################################################################
#CUDA OPTIONS
#LIBS=-ladolc -lgsl -lgslcblas -lfftw3 -llal  -lcudart 
#OBJECTSCUDA := $(patsubst $(SRCDIR)/%,$(ODIRCUDA)/%,$(SOURCESCUDA:.$(SRCEXTCUDA)=.o))
LIBS=-ladolc -lgsl -lgslcblas -lfftw3
OBJECTSCUDA := 
########################################################################

IEXT := h
DEPS:= $(shell find $(IDIR) -type f -name *.$(IEXT))

#CUDA specific header files -- not meant for external use
IEXTCUDA := hu
DEPSCUDA:= $(shell find $(IDIR) -type f -name *.$(IEXTCUDA))

CC=$(CXX)
#CC=g++
CCCUDA=nvcc
#CC=nvcc

.PHONY: all 
all:  Doxyfile $(PROJ_LIB) $(PROJ_PYLIB) $(PROJ_SHARED_LIB) $(MCMC_TOOL)

$(ODIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(DEPS) $(CONFIGFILE)
	$(CC) $(CFLAGS) -c -o $@ $<

#Write pertinent information to config file -- ie filepaths
$(CONFIGFILE): 
	@cp data/Config_template.txt $(CONFIGFILE)
	@pwd | awk '{print "#define GWAT_ROOT_DIRECTORY \""$$1"/\""}' >> $(CONFIGFILE)
	@echo "">>$(CONFIGFILE)
	@echo "#endif">>$(CONFIGFILE)

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

$(PROJ_SHARED_LIB) : $(OBJECTS) $(OBJECTSCUDA) | $(LDIR_LOCAL)
	$(CC) -shared -fopenmp $(OBJECTS) $(OBJECTSCUDA) -o $(LDIR_LOCAL)/$(LOCAL_SHARED_LIB) $(LIBS)

$(LDIR_LOCAL):
	mkdir -p $(LDIR_LOCAL)

$(PROJ_PYLIB): $(PROJ_LIB) $(PROJ_SHARED_LIB) $(PROJ_PYSRC)
	make -C $(PYDIR) 
	
$(TESTFISHER) : $(OBJECTS) $(TESTFISHEROBJ) | $(TESTDIR)
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

$(TEST) : $(OBJECTS) $(OBJECTSCUDA) $(TESTOBJ) | $(TESTDIR)
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

$(TESTDIR):
	mkdir -p $(TESTDIR)

Doxyfile: $(OBJECTS) $(OBJECTSCUDA) $(PROJ_PYLIB)
	doxygen Doxyfile
	#make -C docs/latex
	#cp docs/latex/refman.pdf ./
	#make -C docs/latex clean
	 
	
 
.PHONY: c
c: $(PROJ_LIB) $(PROJ_SHARED_LIB) $(MCMC_TOOL)

.PHONY: test
test: $(TEST) $(PROJ_LIB) $(PROJ_SHARED_LIB) $(PROJ_PYLIB) 

.PHONY: testfisher
testfisher: $(TESTFISHER) $(PROJ_LIB) $(PROJ_SHARED_LIB)

.PHONY: testc
testc: $(TEST) $(PROJ_LIB) $(PROJ_SHARED_LIB)

.PHONY: clean
clean:
	-rm build/*.o 
	-rm build_cuda/*.o 
.PHONY: remove
remove:
	-rm build/*.o build_cuda/* $(TEST) lib/*.a lib/*.so
	-rm include/gwat/GWATConfig.h
	-make -C $(PYDIR) remove

$(MCMC_TOOL): $(OBJECTS)
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)
