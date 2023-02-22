
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
You should now be in the docker image, and you should be in a directory called `/injection/`.

6. To make and run the binary, execute the following commands
```
make 
./bin/injection.exe
```
The output of the code will be saved
under `/data` (which also contains the injection parameters).
It may take some time to run on your computer. 

7. To examine the output, there are routines already packaged in `gw_analysis_tools` and `BayesShip` to unpack and process the raw data coming from the sampler. Included in the `python/` directory of this example is an example notebook that loads the data and plots a corner plot, trace plots, and some useful convergence diagnostics. To launch it, you can run 
```
./jupyter_launch.sh
```
from within the running container used for the sampling. This will host a jupyter kernel that can be accessed through the browser on your local computer at localhost:8989. The token is printed in the output of running the script, but can also be accessed through 
```
jupyter server list
```
Alternatively, you can run the `jupyter_docker_launch.sh` script from a terminal directly on your laptop using 
```
./jupyter_docker_launch.sh
```
which will launch a separate docker container hosted in a separate container. The server than can accessed just like the above method, but at localhost:8990 (note the different port).

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
