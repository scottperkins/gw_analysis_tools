#!/bin/bash

cd /opt/BayesShip/build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make 
make install
cd /opt/gw_analysis_tools/build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make
make install
cd /opt/gw_analysis_tools/examples/

