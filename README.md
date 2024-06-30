# MetalWalls Vitis

Metalwalls is an approach for performing molecular dynamics (MD) simulations in electrochemical systems. While there are a lot of applications available for molecular dynamics, where atoms are being investigated for their interactions and behaviour, there are not much for electrochemical systems, which Metalwalls tries to change. Fundamentally, it is an attempt to calculate the forces on atoms in order to get the updated positions of the atoms in relation to these forces. It is different from regular MD applications since it is using electric potentials and not electric charges. This makes it useful for simulating systems, where the energy is stored in the electrostatic form, for example supercapacitors.

The application was first developed with Fortran in 2007 [[1]](References) and has subsequently undergone constant development to improve functionality and performance. [[2]](References) [[3]](References) Also new accelerated approaches were developed using the Maxeler and the OneAPI toolchains, to accelerate the calculations of a small part of the original code with multiple FPGAs. [[4]](References) This work is based on the existing work and follows the same approach, but utilizes Vitis HLS. We were able to demonstrate the scaling on up to 45 FPGAs.

## Noctua 2

As we developed MetalWalls Vitis with the [Noctua 2 HPC cluster](https://pc2.uni-paderborn.de/systems-and-services/noctua-2), of the [Paderborn Center for Parallel Computing](https://pc2.uni-paderborn.de/) (PC2) in mind, we describe two ways if getting started with the application. The first one describes what we did, so anybody with access to Noctua 2 may replicate our steps by using prepared scripts and our makefile. If you would like to use the application on another cluster please refer to the [Other Cluster](Other-Cluster) sections. To get introduced to the Noctua 2 cluster we suggest referring to the [Getting started](https://upb-pc2.atlassian.net/wiki/spaces/PC2DOK/pages/1903652/Getting+Started) documentation of PC2, as it provides ample information of how to access and operate the cluster. More specific information on how to launch the application on Noctua 2 is provided in the sections below.

### Getting Started

#### Setting up the Environment

Before building or running anything, you need to set up your working environment. This is done by the `env.sh` script. Simply execute:
```
source env.sh
```
to load all of the required modules, download the required Python modules, and set the necessary environment variables.

### Usage

The application has five important targets: The host application, the three FPGA accelerator designs `k0_acc`, `lr_acc`, and `sr_acc`, and the unit testing application. The code of these targets is found in subdirectories of the `src` folder. The build system is implemented as a `Makefile`, which places all artifacts in the `build` directory. This directory is ignored by git and thus, code and artifacts are separated.

#### Running Unit Tests

The commands to build and run the unit tests in emulation mode are:
``` bash
make k0_acc lr_acc sr_acc TARGET=sw_emu # Build the three accelerator designs for emulation in serial
make run_test -j16 TARGET=sw_emu # Build the unit testing executable and run it.
```
The `run_test` target doesn't automatically build the accelerator designs since `make` might otherwise try to rebuild a hardware image that hasn't changed. This way, building the accelerator design is an explicit request. This also makes it possible to build and test only a single design, for example:
``` bash
make sr_acc TARGET=sw_emu # Only build the sr_acc design
make run_test -j128 TARGET=sw_emu TEST_PARAMS="[sr_acc]" # Tell catch2 to only run tests with the [sr_acc] tag.
```
It is also not recommended to build the accelerator designs in parallel since the build commands seem to interact and block each other, leading to failed compilations.

#### Generating reports

Building the synthesis reports is simple. Just run:
``` bash
make reports -j16
```
Unlike the final accelerator designs, building the reports in parallel doesn't seem to have side effects. You will find the reports in the folders `./build/k0_acc/`, `./build/lr_acc`, and`./build/sr_acc`.

#### Synthesizing Designs

Hardware synthesis is done the same way as building emulation binaries by running `make sr_acc`, but without the `TARGET=sw_emu` argument.

#### Running the Application

We introduced the "run"-script ```run_miniapp.sh```, which can be passed to the ```sbatch``` command to allocate a SLURM job, that executes the application for two experiments: ```supercap``` and ```graphene``` and aggregates the metrics in a single file. To execute the application for these two experiments, you can therefore execute:
```
sbatch run_miniapp.sh
```

There is also the possibility to enable the output of the charges in the vtk format for visualisation with ParaView for example. For this the application needs more input files, otherwise it will not run. This can be executed with:

```
sbatch run_miniapp_vtk.sh
```

## Other Cluster

If you would like to use MetalWalls Vitis on a cluster other than Noctua 2, the following sections should give you all instructions necessary to do so. Please keep in mind, that we just can give rough descriptions of what software you need to run the application, but won't be able to go into much detail, as a cluster could have a special way to provide software to their users. By cross referencing the target cluster's documentation on provided software, you should be able to run the appllication though.

### Getting started

#### Setting up the Environment

In order to run the application you need to set up your environment, such that the following software packages are loaded:

- Python3 (Version 3.10.8)
- GNU Compiler Collection (GCC) (Version 12.2.0)
- Xilinx Runtime Library (xrt) (Version 2.15)
- Visualization Toolkit (VTK) (Version 9.2)

When running the test in software emulation, it can be necessary to increase the stack size.

```
ulimit -s 8182
```

## Documentation

###  Assumptions in our code

As we provide three accelerator designs, we assume that the application is launched with at least three MPI ranks, to provide one accelerator for one kernel design each. The rest of our workload distribution handles rank distribution to devices, by usage of MPI node communicators, so the workload gets distributed according to the provided number of ranks. The number of ranks per node must match the number of FPGA devices per node.
### Runscript

As Noctua 2 uses SLURM, the following table provides brief descriptions on the SLURM command flags set via the ```run_miniapp.sh``` script, so you can edit them to suit your individual situation:

For more information about SLURM and flags of ```sbatch``` refer to the [PC2 documentation](https://upb-pc2.atlassian.net/wiki/spaces/PC2DOK/pages/1902952/Running+Compute+Jobs), regarding job submission or the [SLURM documentation](https://slurm.schedmd.com/).

## References
[1] https://gitlab.com/ampere2/metalwalls

[2] Abel Marin-Laflèche, Matthieu Haefele, Laura Scalfi, Alessandro Coretti, Thomas Dufils, Guillaume Jeanmairet, Stewart K. Reed, Alessandra Serva, Roxanne Berthin, Camille Bacon, Sara Bonella, Benjamin Rotenberg, Paul A. Madden, and Mathieu Salanne. MetalWalls: A classical molecular dynamics software dedicated to the simulation of electrochemical systems. Journal of Open Source Software, 5(53):2373, 2020, https://doi.org/10.21105/joss.02373

[3] Alessandro Coretti, Camille Bacon, Roxanne Berthin, Alessandra Serva, Laura Scalfi, Iurii Chubak, Kateryna Goloviznina, Matthieu Haefele, Abel Marin-Laflèche, Benjamin Rotenberg, Sara Bonella, and Mathieu Salanne. MetalWalls: Simulating electrochemical interfaces between polarizable electrolytes and metallic electrodes. 157(18):184801, 2022, https://doi.org/10.1063/5.0101777

[4] Charles Prouveur, Matthieu Haefele, Tobias Kenter, and Nils Voss. FPGA acceleration for HPC supercapacitor simulations. In Proceedings of the Platform for Advanced Scientific Computing Conference, 2023, https://dl.acm.org/doi/10.1145/3592979.3593419
