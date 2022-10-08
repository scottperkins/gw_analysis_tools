---
author: Scott Perkins
date: 2022-08-17
title: Installation Process
---

# Overview

To install this software library, you'll need to do the following steps:

## Install dependencies [^1] 

1. CMake \*
2. FFTW3 \*
3. GSL \*
4. OpenMP \*
5. Eigen3 \*
6. adol-c - Not as common, so probably will need to install from source. Use the appropriate openmp flag if that's desired 
7. HDF5 - Not as common. Might be available on certain systems, but will need to compile with gcc (not clang). This is a headache, but things won't compile properly if HDF5 is not compiled with the same compiler as GWAT
8. BayesShip - Another package from me. Download from git and follow install instructions

\* means the library is very common and is probably available through system package managers (apt/yum/etc) or through HPC infrastructure

[^1]: On systems where you have admin privileges, I would install everything in "/usr/local" unless there's good reason not to. For systems where you don't have admin privileges or if you want to keep the software isolated, I typically put a directory called ".local" into my home directory, and install everything there.

## Install library

1. Download the source code from github.
2. Make a directory called ``build'' in the source directory. 
3. Move into ``build/'' and run:

```bash
cmake .. 
```

4. If you need to modify compile settings (like turning on debugger flags or changing the install prefix), run:

```bash
ccmake .. 
```

Save it and rerun 

```bash
cmake .. 
```

5. Finally, run 

```bash
make 
make install
```






