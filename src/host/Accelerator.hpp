/**
 * Copyright 2023-2024
 *
 * Akanksha Ramakant Choudante <akanksha@mail.upb.de>
 * Amir Kouhpayehzadeh Esfahani <amirk@mail.upb.de>
 * Nils Horsmann <horsmann@mail.upb.de>
 * René Lammert <rlammert@mail.upb.de>
 * Jan-Oliver Opdenhövel <joo@mail.upb.de>
 * Gerrit Pape <papeg@mail.upb.de>
 *
 * This file is part of metalwalls-vitis
 *
 * metalwalls-vitis is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General
 * Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * metalwalls-vitis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with metalwalls-vitis. If not, see
 * <https://www.gnu.org/licenses/>.
 */
#pragma once
#include "xcl2.hpp"
#include <optional>
#include <string>

#include "ComputeUnit.hpp"
#include "Experiment.hpp"

typedef struct kernel_run_times
{
    double start;
    double end;
} kernel_run_times_t;

enum class AcceleratorDesign
{
    K0_ACC,
    LR_ACC,
    SR_ACC,
};

/**
 * @brief Wrapper for the FPGA accelerator.
 *
 * An object of this class is a representation of physical, programmed FPGA
 * accelerator that can execute computations. This is done by wrapping the
 * internal context and command queue, as well as the FPGA program and the
 * individual kernels.
 */
class Accelerator
{
  public:
    /**
     * @brief Set up an FPGA accelerator
     *
     * This constructor sets up the FPGA accelerator with the given index with
     * the requested design and sets up the context and queue for computations.
     *
     * The path to the design's binary path is constructed from the environment
     * variables `MW_BINARY_DIR` and `XCL_EMULATION_MODE`. In general, the
     * expected path of the binaries is
     * `${MW_BINARY_DIR}/${design}_${XCL_EMULATION_MODE}.xclbin`. If
     * `MW_BINARY_DIR` is not set, the current working directory (i.e. `.`) is
     * used, and if `XCL_EMULATION_MODE` is not set, the hardware design will be
     * used (i.e. `hw`). For example, if you want to load the SR design in
     * `./build/` for software emulation, you would set the variables
     * `MW_BINARY_DIR=./build/` and `XCL_EMULATION_MODE=sw_emu`. Then, the design
     * in
     * `./build/sr_acc_sw_emu.xclbin` is loaded onto the FPGA.
     *
     * @param design The design to load onto the FPGA.
     * @param device_index The index of the accelerator to wrap.
     */
    Accelerator(Experiment experiment, AcceleratorDesign design, int unsigned device_index);

    /**
     * @brief Run the k0 kernel.
     *
     * This function will set the scalar parameters, copy the vector parameters
     * to the FPGA, and enqueue the k0 kernel. After the kernel is finished, the
     * result is copied to `k0`.
     *
     * This function may crash if the program loaded onto the FPGA does not
     * contain the k0 kernel.
     *
     * @param num_atoms The total number of atoms in the experiment.
     * @param n_start The index of the first atom to process.
     * @param n_end The index of the last atom to process.
     * @param alpha The kernel parameter alpha.
     * @param alphasq The kernel parameter alphasq.
     * @param sqrpialpha The kernel parameter sqrpialpha.
     * @param k0PotFactor The kernel parameter k0PotFactor.
     * @param z The z values for every atom. This has to be an array with at
     * least `num_atoms` values.
     * @param q The q values for every atom. This has to be an array with at
     * least `num_atoms` values.
     * @param k0 Output vector. This has to be an array with at least `num_atoms`
     * values.
     */
    kernel_run_times_t k0(int n_start, int n_end, double const *q, double *k0);

    /**
     * @brief Run the short-range kernel.
     *
     * This function will set the scalar parameters, copy the vector parameters
     * to the FPGA, and enqueue the short-range kernel. After the kernel is
     * finished, the result is copied to `lr`.
     *
     * The position coordinates are expected to be stored in an interleaved
     * fashion: For the i-th atom, the x coordinate is stored in `xyz[4*i]`, the
     * y coordinate is stored in `xyz[4*i+1]`, and the z coordinate is stored in
     * `xyz[4*i+2]`; `xyz[4*i+3]` is a padding element. Therefore, the `xyz`
     * array needs to contain at least `4*numWall1` elements.
     *
     * This function may crash if the program loaded onto the FPGA does not
     * contain the short-range kernel.
     *
     * @param a The kernel parameter a.
     * @param alpha The kernel parameter alpha.
     * @param b The kernel parameter b.
     * @param param_eta The kernel parameter eta.
     * @param numWall1 The total number of atoms in the experiment.
     * @param rcutsq The kernel parameter rcutsq.
     * @param xyz The position coordinates of every atom. This has to be an array
     * with at least `4*numWall1` values.
     * @param q The q values for every atom. This has to be an array with at
     * least `numWall1` values.
     * @param sr Output vector. This has to be an array with at least `numWall1`
     * values.
     */
    kernel_run_times_t sr(const int n_start, const int n_end, double const *q, double *sr);

    /**
     * @brief Run the long-range kernel.
     *
     * This function will set the scalar parameters, copy the vector parameters
     * to the FPGA, and enqueue the long-range kernel. After the kernel is
     * finished, the result is copied to `lr`.
     *
     * This function may crash if the program loaded onto the FPGA does not
     * contain the long-range kernel.
     *
     * @param num_atoms The number of atoms.
     * @param mode_start
     * @param mode_end
     * @param kx_start
     * @param ky_start
     * @param kz_start
     * @param kyMax
     * @param kzMax
     * @param alphaconst
     * @param rksqmax
     * @param twopiarec
     * @param twopibrec
     * @param twopicrec
     * @param lrPotFactor
     * @param x Input vector This has to be an array with at least `num_atoms`
     * elements.
     * @param y Input vector This has to be an array with at least `num_atoms`
     * elements.
     * @param z Input vector This has to be an array with at least `num_atoms`
     * elements.
     * @param q Input vector This has to be an array with at least `num_atoms`
     * elements.
     * @param lr Output vector. This has to beSR, I had to raise the precision requirements to 1e-8 so that the Graphene
     * tests pass. However, The LR test cases fail entirely. There are even instances where the LR kernel produces zeros
     * and the reference data contains something like -27. That's not good and this needs to be investigated. an array
     * with at least `num_atoms` elements.
     */
    kernel_run_times_t lr(int mode_start, int mode_end, double *q, double *lr);

    /**
     * @brief Return the wrapped command queue.
     * @return The wrapped command queue.
     */
    inline cl::CommandQueue get_queue()
    {
        return queue;
    }

    /**
     * @brief Return the wrapped context.
     * @return The wrapped context.
     */
    inline cl::Context get_context()
    {
        return context;
    }

    /**
     * @brief Return two boundries for a workload, which is devisible by the given chunk_size.
     * @return The wrapped context.
     */
    static std::vector<size_t> get_boundaries(size_t workload, size_t number, size_t chunk_size, size_t rank);

    static std::vector<size_t> get_distribution(size_t total_rows, size_t number, size_t chunk_size);

  private:
    Experiment experiment;
    AcceleratorDesign design;
    std::optional<ComputeUnit> k0_cu;
    std::optional<ComputeUnit> sr_cu;
    std::optional<ComputeUnit> lr_cu;
    cl::CommandQueue queue;
    cl::Context context;
    cl::Program program;
};
