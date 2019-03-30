#!/bin/bash

gcc-7 -o test.o -c test.c
g++-7 -L. -o test test.o -lgwanalysistools -ladolc -lgsl -lm -lfftw3 -lgslcblas -llapack
