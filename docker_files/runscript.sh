#!/bin/bash

############################################
#Example commands to build the docker image, name it docker-gwat-dev, 
#and launch a container with a volume linked to the current directory.
#Must be run in /${PATH_TO_GWAT_DIR}/
#Note -- the image file only needs to be run once, and the container only has to be launched once per session
############################################


#build
docker build --tag docker-gwat-dev ./docker_files/
#run
docker run -v $(pwd):/gwat -it docker-gwat-dev /bin/bash
