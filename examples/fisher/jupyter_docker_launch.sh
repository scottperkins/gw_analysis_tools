#!/bin/bash
docker run  \
        --platform linux/amd64 \
	-p 8990:8989\
	--rm \
	--name gw-injection-example-jupyter \
	-w /fisher/ \
        -v $(pwd)/src:/fisher/src \
        -v $(pwd)/data:/fisher/data \
        -v $(pwd)/python:/fisher/python \
        -v $(pwd)/makefile:/fisher/makefile \
	scottperkins/gwat.jupyter 

