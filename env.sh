# Load software
module reset
module load fpga toolchain vis
module load xilinx/xrt/2.16 gompi/2023b VTK/9.3.0-foss-2023b

# Set environment variables
export XRT_INI_PATH=./xrt.ini

#increase stack size
ulimit -s 8182
