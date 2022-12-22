---
author: Scott Perkins
date: 2022-12-22
title: Writing Custom Docker Images: An Example Workflow
---

# Writing Custom Docker Images: Example Workflows

## Developing in a Docker container


In the process of writing this C++ library, the authors have embraced Docker containers as a flexible way to portably modify and run this library, both for production on scientific HPC clusters and on local computers for development.
This is by no means a necessity, but an example of this workflow is shown below for the benefit of those running into compatibility/portability issues.
This is, of course, no long term solution and no replacement for developing reproducible workflows, but can be a helpful mitigation effort to not slow down the science.

The goal of this example is to show how one might develop the software library (modify the source code, compile it, and run tests) by using a Docker container.
To do this, we will use the publicly maintained base image, "scottperkins/gw.base.deps", hosted on DockerHub by the principle developer.
Note, there are publicly maintained pre-compiled images to use on DockerHub, but we will be making our own image to allow for developing the source code.

For this example, we will assume the user is in a directory with the following structure, where "BayesShip" represents any necessary dependency not included in the image "scottperkins/gw.base.deps".
```
/root-directory/
|-------BayesShip/
|-------other-project-code/
```
We will assume the user is in the "root-directory" to begin with.

1. Download the source code from the repository
```bash
git pull <git-repo-name>
```
2. Make any desired changes to the source code.
3. Write a docker-compose.yml file, following the following template and save it in "root-directory",
```bash
version: "3"
services:
  project:
    container_name: <name>
    image: scottperkins/gw.base.deps
    ports:
        - "8081:8081"
    tty: true
    stdin_open: true
    volumes:
        - $PWD/BayesShip:/opt/BayesShip
        - $PWD/gw_analysis_tools/:/opt/gw_analysis_tools
        - $PWD/other-project-code/:/opt/other-project-code
```
4. We now expect the following file structure, and that the user is in "root-directory":
```
/root-directory/
|-------BayesShip/
|-------gw_analysis_tools/
|-------other-project-code/
|-------docker-compose.yml
```
5. Build the image with the following command
```bash
docker-compose up -d
```
6. Now, one can login to the running container with
```bash
docker exec -it <name> bash
```
7. Inside the container, one can install any additional dependencies one needs, like python packages and other software libraries through the usual managers, like pip/apt/yum/etc.
8. For clarity, the install location for any compiled software should generally be on non-mounted volumes, but this is up to the discretion of the developer. 
9. If the necessary libraries and executables have been installed in the container (not in mounted volumes), then one can save the image file with 
```bash
docker save myimage:latest | gzip > myimage_latest.tar.gz
```
(custom_docker_standalone)=
## Writing a standalone Docker image 

This example is focused on the use case where one wants to use this software on a scientific HPC cluster where singularity might be used. 
One method of doing this is to first create a Docker image, outlined below.
For this example, we will assume the user is in a directory with the following structure, where "BayesShip" represents any necessary dependency not included in the image "scottperkins/gw.base.deps".
```
/root-directory/
|-------BayesShip/
```
We will assume the user is in the "root-directory" to begin with.

1. Download the source code from the repository
```bash
git pull <git-repo-name>
```
2. Make any desired changes to the source code.
3. Write a Dockerfile, following the following template and save it in "root-directory",
```bash
FROM scottperkins/gw.base.deps:latest # Pull the image from the online repository with all the dependencies pre-installed


# Install any other dependencies (like custom repositories)
RUN mkdir /opt/BayesShip # Make the directory for the source code *in the container*
RUN mkdir /opt/BayesShip/build # Make a build directory
COPY BayesShip/ /opt/BayesShip/ # Copy all source files into the container
WORKDIR /opt/BayesShip/build/ # Make the current working directory the build directory
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr/local ../ # Configure the installation
RUN make -j 4 install # Compile and install the library

# Install the main source code
RUN mkdir /opt/gw_analysis_tools
RUN mkdir /opt/gw_analysis_tools/build
COPY gw_analysis_tools/ /opt/gw_analysis_tools/
WORKDIR /opt/gw_analysis_tools/build/
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr/local ../
RUN make -j 4 
RUN make install 
WORKDIR /

# Remove the source files *from the container* in case one wants to host the image on DockerHub (removing the source files allows for maintaining privacy of the source code until public)
RUN rm -r /opt/gw_analysis_tools/
RUN rm -r /opt/BayesShip/
```
4. We now expect the following file structure, and that the user is in "root-directory":
```
/root-directory/
|-------BayesShip/
|-------gw_analysis_tools/
|-------<Dockerfile-name>
```
5. Build the image with the following command
```bash
docker build -f <Dockerfile-name> -t <tag> .
```
If the user wants to host on DockerHub, the \<tag\> name should be \<DockerHub-ID\>/\<tag-name\>.
6. If desired, the user can now push the image to DockerHub (which is *generally* public):
```bash
docker push <tag>
```
7. Now one can save the image file with 
```bash
docker save myimage:latest | gzip > myimage_latest.tar.gz
```

