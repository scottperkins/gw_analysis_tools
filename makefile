SYSIDIR=${CDIR}
LDIR=${LD_LIBRARY_PATH}


SRCDIR=src
ODIR=build
IDIR=include
LDIR_LOCAL=lib

PYDIR=gw_analysis_tools_py
PYSRC=mcmc_routines_ext.pyx waveform_generator_ext.pyx
#PROJ_PYSRC=$(PYDIR)/$(PYSRC)
PROJ_PYSRC=$(addprefix $(PYDIR)/src/,$(PYSRC))

PYLIB=mcmc_routines_ext.cpp waveform_generator_ext.cpp
PROJ_PYLIB=$(addprefix $(PYDIR)/,$(PYLIB))

LIBS=-ladolc -lgsl -lgslcblas -lfftw3
LOCAL_LIB=libgwanalysistools.a
PROJ_LIB=$(addprefix $(LDIR_LOCAL)/,$(LOCAL_LIB))

TESTSRC=testing/test.cpp
TESTOBJ=testing/test.o
TEST=bin/exe.a



CFLAGS=-I$(IDIR) -I$(SYSIDIR) -Wall -fPIC -g -O3
LFLAGS=-L$(LDIR) 

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(ODIR)/%,$(SOURCES:.$(SRCEXT)=.o))

IEXT := h
DEPS:= $(shell find $(IDIR) -type f -name *.$(IEXT))

CC=g++

.PHONY: all
all:  Doxyfile $(PROJ_LIB) $(PROJ_PYLIB)

$(ODIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(DEPS) 
	$(CC) $(CFLAGS) -c -o $@ $<

$(TESTOBJ): $(TESTSRC)
	$(CC) $(CFLAGS) -c -o $@ $<

$(PROJ_LIB) : $(OBJECTS)
	ar rcs $(LDIR_LOCAL)/$(LOCAL_LIB) $(OBJECTS)

$(PROJ_PYLIB): $(PROJ_LIB) $(PROJ_PYSRC)
	make -C $(PYDIR) 
	
$(TEST) : $(OBJECTS) $(TESTOBJ)
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

Doxyfile: $(OBJECTS)
	doxygen Doxyfile
 
.PHONY: test
test: $(TEST) $(PROJ_LIB) $(PROJ_PYLIB) 

.PHONY: clean
clean:
	rm build/*.o 
.PHONY: remove
remove:
	rm build/*.o $(TEST) lib/*.a
	make -C $(PYDIR) remove
#.PHONY: docs
#docs:
#	doxygen Doxyfile
