#!/bin/bash
docker run  \
	-it \
        --platform linux/amd64 \
	-p 8990:8989\
	--rm \
	--name gw-injection-example-jupyter \
	-w /injection/ \
        -v $(pwd)/src:/injection/src \
        -v $(pwd)/data:/injection/data \
        -v $(pwd)/bin:/injection/bin \
        -v $(pwd)/build:/injection/build \
        -v $(pwd)/python:/injection/python \
        -v $(pwd)/makefile:/injection/makefile \
        -v $(pwd)/jupyter_launch.sh:/injection/jupyter_launch.sh \
	scottperkins/gwat.jupyter ./jupyter_launch.sh

echo | docker logs gw-injection-example-jupyter
