PYSRC=./src/*.pyx
PYLIB=./*.so
SETUPSCRIPT=setup.py
PROJ_LIB=../lib/libgwat.so

.PHONY: all
all: $(PYLIB) $(SETUPSCRIPT)

$(PYLIB) : $(PYSRC) $(PROJ_LIB) $(SETUPSCRIPT)
	python3 $(SETUPSCRIPT) build_ext -i
.PHONY: remove
remove:
	rm ./*.so

.PHONY: clean
clean:
	rm ./*.so
