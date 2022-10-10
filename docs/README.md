
A suite of tools useful for doing statistical studies  on gravitational wave science, including routines useful in MCMC studies, wave template generation, Fisher analysis, etc. Written in C++.

# Installation
You can manually install this software directly on the computer of interest, or you can run this through docker containers.

## Manual Installation

### Install dependencies [^1] 

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

### Install library

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

## Docker Usage

Several public images with this software already installed are maintained on DockerHub:

1. scottperkins/mcmc_dev 
	
	- An image with the mcmc_dev branch of GWAT installed. *This should be the default option for new users*

2. scottperkins/parityviolation
	
	- An image with the parityViolation branch of GWAT installed

3. scottperkins/gwat.jupyter

	- An image with mcmc_dev branch of GWAT installed and automatically launches a jupyter server on port 8989.

	- Useful for doing post-MCMC analysis, as there are python wrappers that unpack the data installed on this image

	- For usage notes on this, see [jupyter notes](extraDocs/jupyterDocker/index.md)


