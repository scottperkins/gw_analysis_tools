
# Running with docker: two minimal working examples 

In these two examples, we download the image
[scottperkins/gwat](https://hub.docker.com/r/scottperkins/gwat),
and run an example usage of `gwat` that performs parameter
estimation on a gravitational wave injection.
The first example shows how to do this on your local computer
using docker,
the second shows how to do this with `singularity`, which the
container application that runs on clusters (instead of docker).

## Running on your local computer

1. Assuming you already have docker desktop installed,
start docker desktop, e.g. in the terminal type
```
systemctl --user start docker-desktop
``` 

2. Make a directory, e.g. `GW-Analysis`, and cd to that directory. 

3. Clone [gw_analysis_tools](https://github.com/scottperkins/gw_analysis_tools)
to the directory. cd to `gw_analysis_tools/examples/injection_IMRPhenomD`.

4. There should be bash script called `docker_launch.sh`.

Running this file downloads and runs the docker image 
[scottperkins/gwat](https://hub.docker.com/r/scottperkins/gwat)
from Docker Hub. 

5. Run the script 
```
bash docker_launch.sh
```
You should now be in the docker image, and there should be a directory
called `injection`. cd to that directory, and type `make` to build
the binary.

6. You can run the binary by
```
./bin/injection.exe
```
The output of the code will be saved
under `/data` (which also contains the injection parameters).
It may take some time to run on your computer. 

## Running on a cluster

1. Make sure singularity is available, e.g. through `module load singularity`

2. As above, make and cd to directory, e.g. `GW-Analysis`.

3. Clone [gw_analysis_tools](https://github.com/scottperkins/gw_analysis_tools)
to the directory. cd to `gw_analysis_tools/examples/injection_IMRPhenomD`.

4. There should be bash script called `singularity_launch.slurm`.
You may need to change some of the SBATCH commands for your local cluster.

Running this file downloads and runs the docker image 
[scottperkins/gwat](https://hub.docker.com/r/scottperkins/gwat)
from Docker Hub, and then runs the `runscript_for_slurm.sh`
script, which makes and runs the executable within the image.

5. Launch with 
```
slurm singularity_launch.slurm
```

Output should be saved in the relevant `data` directory
in your scratch folder (or elsewhere if you changed the .slurm file). 
