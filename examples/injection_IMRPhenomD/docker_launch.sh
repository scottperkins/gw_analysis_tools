#!/bin/bash
docker run -it \
        --platform linux/amd64 \
        -v $(pwd)/src/:/injection/src \
        -v $(pwd)/data/:/injection/data \
        -v $(pwd)/bin/:/injection/bin \
        -v $(pwd)/build/:/injection/build \
        -v $(pwd)/makefile/:/injection/makefile \
	scottperkins/gwat
