
#PYSRC=mcmc_routines_ext.pyx waveform_generator_ext.pyx
PYSRC=./src/*.pyx
#PYLIB = mcmc_routines_ext.cpp waveform_generator_ext.cpp
PYLIB=./*.so
PROJ_LIB=../lib/libgwanalysistools.a

.PHONY: all
all: $(PYLIB)

$(PYLIB) : $(PYSRC) $(PROJ_LIB)
	python setup.py build_ext --inplace	
.PHONY: remove
remove:
	rm ./*.so

