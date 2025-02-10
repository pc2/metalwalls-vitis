#!/bin/bash

#SBATCH -t 00:30:00
#SBATCH -n 15
#SBATCH -N 5
#SBATCH -J "metalwalls"
#SBATCH -p fpga
#SBATCH -A hpc-lco-kenter
#SBATCH --constraint xilinx_u280_xrt2.16

## Load environment modules
source env.sh

EXPERIMENTS="supercap graphene"

for e in $EXPERIMENTS
do
    if [ ! -e "scaling_${e}.json" ]; then
        echo "[]" > scaling_${e}.json
    fi
done

mkdir -p tmp && cd tmp
for device in "0000:a1:00.1" "0000:81:00.1" "0000:01:00.1"
do
    srun -n $SLURM_NNODES --spread-job xbutil reset --force --device $device
done
cd ../ && rmdir tmp

function run_metalwalls {
    export MW_BINARY_DIR=./build/

    srun -l -n $2 ./build/metalwalls ./experiments/$1/input_f90.dat || exit 1
    cp charges.out charges_$1.out
    cp metrics.json metrics_$1.json
    ./scripts/add_benchmark_run.py scaling_$1.json $2
}

for n_ranks in `seq 13 15`
do
    for e in $EXPERIMENTS
    do
        run_metalwalls $e $n_ranks || exit 1
    done
done
