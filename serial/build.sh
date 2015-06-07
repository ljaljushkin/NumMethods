#!/usr/bin/env bash

if [ -d "__cmake" ]; then
    echo "------------------->cleaning temp cmake folder<-------------------"
    rm -rf __cmake/*
else
    echo "------------------->creating temp cmake folder<-------------------"
    mkdir __cmake
fi

cd __cmake

echo "------------------->starting cmake build<-------------------------"
export CXX=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
cmake ../
make

echo "------------------->running solution<-----------------------------"
cd ../bin/
./solution input output time
