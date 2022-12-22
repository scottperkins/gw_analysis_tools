---
author: Scott Perkins
date: 2022-12-23
title: Running with Singularity 
---

# Running with Singularity 

There are a few ways to run this code through a singularity container (which are much more common on scientific HPC clusters than Docker), and below is outlined an example workflow one might use (but definitely not the only way one can accomplish this).
This example will be done in the context of running on a cluster implementing SLURM, but that can be easily modified to suit the desired environment.

## via Docker
The method built around Docker will utilize the maintained images on DockerHub (e.g. scottperkins/gwat and scottperkins/gw.base.deps primarily).
If one creates one's own [docker image with custom development](custom_docker_standalone), simply replace all the instances of 
```bash
docker://scottperkins/gwat
```
with 
```bash
docker-archive://Docker-image.tgz
```

### Run the pre-packaged executables 

In this section, the steps to run the pre-packaged executables will be outlined.
For example, we will show how to use the "mcmc_gw_tool_v2" executable to run on pre-downloaded (and pre-processed, if applicable) data from GWOSC.
Assuming one is on the remote server, there should be all the necessary data, 
```
/project-dir
|------param_file.dat 
|------slurm_submission.sh
|------extra_data.dat
|------data/
|------------PSD.dat
|------------STRAIN.dat
```

In the param_file.dat, there should be all the necessary information to run in the container (appropriate file paths, etc.). 
A template can be found in "data/sample_config_files".

In the slurm_submission.sh, there should be the following steps:
```bash
#!/bin/sh

#SBATCH --time=<TIME>
#SBATCH --ntasks=<TASKS>
#SBATCH --job-name=<NAME>
#SBATCH --partition=<PARTITION>
#SBATCH --output=%x-%j.output
#SBATCH --error=%x-%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<EMAIL>

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load singularity
cd /scratch/users/$USER                         # go to scratch directory (WILL change for each user)
mkdir -p $SLURM_JOB_NAME/$SLURM_JOB_ID                # create job unique directory
cd $SLURM_JOB_NAME/$SLURM_JOB_ID                     # switch to it
mkdir data


# Copy in the data
cp $HOME/path/to/project-dir/data/* ./data/
cp $HOME/path/to/project-dir/param_file.dat ./
cp $HOME/path/to/project-dir/extra_data.dat ./
export SINGULARITY_CACHEDIR=$HOME/scratch/.singularity/ # Some people may need to put the cache in a specific place, as home directories don't always have enough space


# Launch the command. Note the bind option to mount the current directory (scratch project folder) to the container. The "pwd" option tells the container to run the following command from that directory in the container.
singularity exec  --bind $(pwd):/opt/project --pwd /opt/project docker://scottperkins/mcmc_dev mcmc_gw_tool_v2 param_file.dat
```

### Compile custom software linked to the compiled libraries

If one wants to run custom code *linked* to the pre-compiled (or custom, self-compiled) libraries, one can follow the following steps.
Working in the following directory structure
```
/project-dir
|------source-file.cpp 
|------runscript.cpp 
|------makefile
|------slurm_submission.sh
|------extra_data.dat
```

The source-file.cpp file has the custom software that will link to the existing libraries.
The file runscript.sh has the commands to compile and run the software, for example
```bash
#!/bin/bash

make 
./source-file.exe
```
The "make" utilities are a common configuration/compilation tool, but this could be anything from manually compiling to "cmake".
In the slurm_submission.sh, there should be the following information:
```bash
#!/bin/sh

#SBATCH --time=<TIME>
#SBATCH --ntasks=<TASKS>
#SBATCH --job-name=<NAME>
#SBATCH --partition=<PARTITION>
#SBATCH --output=%x-%j.output
#SBATCH --error=%x-%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<EMAIL>

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load singularity
cd /scratch/users/$USER                         # go to scratch directory (WILL change for each user)
mkdir -p $SLURM_JOB_NAME/$SLURM_JOB_ID                # create job unique directory
cd $SLURM_JOB_NAME/$SLURM_JOB_ID                     # switch to it
mkdir data


# Copy in the data
cp $HOME/path/to/project-dir/runscript.sh ./
cp $HOME/path/to/project-dir/source-file.cpp ./
cp $HOME/path/to/project-dir/makefile ./
cp $HOME/path/to/project-dir/extra_data.dat ./
export SINGULARITY_CACHEDIR=$HOME/scratch/.singularity/ # Some people may need to put the cache in a specific place, as home directories don't always have enough space


# Launch the command. Note the bind option to mount the current directory (scratch project folder) to the container. The "pwd" option tells the container to run the following command from that directory in the container.
singularity exec --no-home  --bind $(pwd):/opt/project --pwd /opt/project docker://scottperkins/mcmc_dev /bin/bash ./runscript.sh 
```
Note, the option "--no-home" stops singularity from binding the user's home directory in the container. This is *optional*, but if the user is experiencing problems with multiple installations of software, this helps isolate the container from the host system (hopefully resolving these issues).


## via Docker and Singularity Sandbox

TBD


