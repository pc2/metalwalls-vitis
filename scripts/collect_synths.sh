#!/usr/bin/bash

make clean
mkdir -p build
wd=`pwd`
for kernel in k0 sr lr cg
do
    path="${wd}_${kernel}"
    cp -r $path/build/*xclbin build
done
