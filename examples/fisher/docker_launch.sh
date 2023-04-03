#!/bin/bash
docker run -it \
        --platform linux/amd64 \
	--rm \
	-w /fisher/ \
	-p 8989:8989 \
        -v $(pwd)/src/:/fisher/src \
        -v $(pwd)/data/:/fisher/data \
        -v $(pwd)/makefile:/fisher/makefile \
	scottperkins/gwat
