#!/usr/bin/bash

make clean
wd=`pwd`
for kernel in k0 sr lr cg
do
    cd ..
    path="${wd}_${kernel}"
    rm -rf $path
    cp -r ${wd} $path
    cd $path
    sbatch -J "${kernel}" ./scripts/synth.sh ${kernel}_acc
    cd ${wd}
done
