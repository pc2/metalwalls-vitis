#!/bin/bash

#SBATCH -t 00:30:00
#SBATCH -n 39
#SBATCH -N 13
#SBATCH -J "metalwalls"
#SBATCH -p fpga
#SBATCH -A hpc-lco-kenter
#SBATCH --constraint xilinx_u280_xrt2.15

## Load environment modules
source env.sh
make build/metalwalls -j128

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

    rm -rf visualisation
    mkdir visualisation

    srun -n $2 ./build/metalwalls ./experiments/$1/input_f90.dat ./experiments/$1/data.inpt ./experiments/$1/runtime.inpt || exit 1

    cp charges.out charges_$1_$2.out
    cp metrics.json metrics_$1_$2.json
    rm -rf visualisation_$1_$2 
    mv visualisation visualisation_$1_$2
    ./scripts/add_benchmark_run.py scaling_$1.json $2
}

for n_ranks in `seq 37 39`
do
    for e in $EXPERIMENTS
    do
        run_metalwalls $e $n_ranks || exit 1
    done
done
