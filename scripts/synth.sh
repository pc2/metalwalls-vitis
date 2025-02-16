#!/usr/bin/bash
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -q fpgasynthesis
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --mail-type=ALL

source env.sh

make $@
