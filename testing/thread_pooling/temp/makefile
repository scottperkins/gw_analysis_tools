SYSIDIR=${CDIR}
LDIR=${LD_LIBRARY_PATH}


SRCDIR=src
ODIR=build
IDIR=include
LIBS=-ladolc -lgsl -lgslcblas -lfftw3 -lpthread

OUT=bin/exe.a

CFLAGS=-I$(IDIR) -I$(SYSIDIR) -Wall -fPIC -g 
LFLAGS=-L$(LDIR) 

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(ODIR)/%,$(SOURCES:.$(SRCEXT)=.o))

IEXT := h
DEPS:= $(shell find $(IDIR) -type f -name *.$(IEXT))

#CC=g++
CC=g++-7

$(ODIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(DEPS) 
	$(CC) $(CFLAGS) -c -o $@ $<

#$(OUT) : $(OBJECTS)
#	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)
$(OUT) : $(OBJECTS)
	$(CC)  -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	rm build/*.o 
.PHONY: remove
remove:
	rm build/*.o $(OUT)
.PHONY: docs
docs:
	doxygen Doxyfile

 
