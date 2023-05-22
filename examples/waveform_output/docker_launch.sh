#!/bin/bash
docker run -it \
        --platform linux/amd64 \
	--rm \
	-w /waveform/ \
	-p 8989:8989 \
        -v $(pwd)/src/:/waveform/src \
        -v $(pwd)/data/:/waveform/data \
        -v $(pwd)/makefile:/waveform/makefile \
	scottperkins/gwat
