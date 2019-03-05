PYSRC=./src/*.pyx
PYLIB=./*.so
PROJ_LIB=../lib/libgwanalysistools.a

.PHONY: all
all: $(PYLIB)

$(PYLIB) : $(PYSRC) $(PROJ_LIB)
	python setup.py build_ext --inplace	
.PHONY: remove
remove:
	rm ./*.so

