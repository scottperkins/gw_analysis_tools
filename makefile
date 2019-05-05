
SRCDIR=src
ODIR=build
IDIR=include
LDIR_LOCAL=lib

PYDIR=gw_analysis_tools_py
PYSRC=mcmc_routines_ext.pyx waveform_generator_ext.pyx
PROJ_PYSRC=$(addprefix $(PYDIR)/src/,$(PYSRC))

PYLIB=mcmc_routines_ext.cpp waveform_generator_ext.cpp
PROJ_PYLIB=$(addprefix $(PYDIR)/,$(PYLIB))

LIBS=-ladolc -lgsl -lgslcblas -lfftw3 -llal
LOCAL_LIB=libgwanalysistools.a
PROJ_LIB=$(addprefix $(LDIR_LOCAL)/,$(LOCAL_LIB))

TESTSRC=testing/test.cpp
TESTOBJ=testing/test.o
TESTDIR=bin
TEST=$(addprefix $(TESTDIR)/,exe.a)



#CFLAGS=-I$(IDIR) -I/opt/lalsuite/lalsimulation/src -I/opt/lalsuite/include -Wall -fPIC -g -O3 -std=c++11
CFLAGS=-I$(IDIR) -fopenmp  -fPIC -g -O2 -std=c++11
LFLAGS= -fopenmp
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(ODIR)/%,$(SOURCES:.$(SRCEXT)=.o))

IEXT := h
DEPS:= $(shell find $(IDIR) -type f -name *.$(IEXT))

#CC=g++-7
CC=g++
#CC=nvcc

.PHONY: all
all:  Doxyfile $(PROJ_LIB) $(PROJ_PYLIB)

$(ODIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(DEPS) 
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJECTS): | $(ODIR)

$(ODIR):
	mkdir -p $(ODIR)

$(TESTOBJ): $(TESTSRC)
	$(CC) $(CFLAGS) -c -o $@ $<

$(PROJ_LIB) : $(OBJECTS) | $(LDIR_LOCAL)
	ar rcs $(LDIR_LOCAL)/$(LOCAL_LIB) $(OBJECTS)

$(LDIR_LOCAL):
	mkdir -p $(LDIR_LOCAL)

$(PROJ_PYLIB): $(PROJ_LIB) $(PROJ_PYSRC)
	make -C $(PYDIR) 
	
$(TEST) : $(OBJECTS) $(TESTOBJ) | $(TESTDIR)
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

$(TESTDIR):
	mkdir -p $(TESTDIR)

Doxyfile: $(OBJECTS)
	doxygen Doxyfile
 
.PHONY: c
c: $(PROJ_LIB)

.PHONY: test
test: $(TEST) $(PROJ_LIB) $(PROJ_PYLIB) 

.PHONY: testc
testc: $(TEST) $(PROJ_LIB)  

.PHONY: clean
clean:
	rm build/*.o 
.PHONY: remove
remove:
	rm build/*.o $(TEST) lib/*.a
	make -C $(PYDIR) remove
