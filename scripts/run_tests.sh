#!/bin/bash

#SBATCH -t 00:30:00
#SBATCH -n 3
#SBATCH -N 1
#SBATCH -J "metalwalls"
#SBATCH -p fpga
#SBATCH -A hpc-lco-kenter
#SBATCH --constraint xilinx_u280_xrt2.15

## Load environment modules
source env.sh
make build/metalwalls -j128

mkdir -p tmp && cd tmp
for device in "0000:a1:00.1" "0000:81:00.1" "0000:01:00.1"
do
    srun -n $SLURM_NNODES --spread-job xbutil reset --force --device $device
done
cd ../ && rmdir tmp

make run_test
