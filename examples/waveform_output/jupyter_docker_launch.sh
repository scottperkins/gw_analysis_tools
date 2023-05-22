#!/bin/bash
docker run  \
        --platform linux/amd64 \
	-p 8990:8989\
	--rm \
	--name waveform-example-jupyter \
	-w /waveform/ \
        -v $(pwd)/src:/waveform/src \
        -v $(pwd)/data:/waveform/data \
        -v $(pwd)/python:/waveform/python \
        -v $(pwd)/makefile:/waveform/makefile \
	scottperkins/gwat.jupyter 

