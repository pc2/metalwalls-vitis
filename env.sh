# Load software
module reset
module load lang fpga compiler toolchain vis
module load Julia/1.10.0-linux-x86_64 Python/3.10.4-GCCcore-11.3.0-bare xilinx/xrt/2.15 GCC/11.3.0 gompi/2022a VTK/9.2.2-foss-2022a

# Set environment variables
export XRT_INI_PATH=./xrt.ini

#increase stack size
ulimit -s 8182
