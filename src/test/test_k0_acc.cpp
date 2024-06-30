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
#include "../host/Accelerator.hpp"
#include "../host/Experiment.hpp"
#include "../host/xcl2.hpp"
#include "KernelRun.hpp"
#include "catch_amalgamated.hpp"
#include "cpu_reference.hpp"

#define BEST_PRECISION 1e-9

const int vector_size = 8;

void test_k0_accelerator_with_extracted_data(const char *experiment_name, const char *target, int start, int end,
                                             double precision)
{
    Accelerator acc(Experiment(experiment_name), AcceleratorDesign::K0_ACC, 0);

    KernelRuns runs(experiment_name, "k0", target);
    runs.deserialize();

    if (end == -1)
    {
        end = runs.iterations;
    }
    for (int i = start; i < end; i++)
    {
        double mean_error = 0.0;
        double min_error = std::numeric_limits<double>::infinity();
        double max_error = -std::numeric_limits<double>::infinity();
        double l2_distance = 0.0;

        std::vector<double, aligned_allocator<double>> kernel_out(runs.num);

        kernel_run_times_t cl_runtime = acc.k0(runs.start, runs.end, runs.q[i].data(), kernel_out.data());

        for (int j = 0; j < runs.num; j++)
        {
            double error = std::abs(kernel_out[j] - runs.result[i][j]);
            mean_error += error / runs.num;
            min_error = std::min(min_error, error);
            max_error = std::max(max_error, error);
            l2_distance += error * error;
        }
        l2_distance = std::sqrt(l2_distance);

        std::cout << "mean_error: " << mean_error << ", min_error: " << min_error << ", max_error: " << max_error
                  << ", l2_distance: " << l2_distance << std::endl;
        std::cout << "measured CL event runtime = " << cl_runtime.end - cl_runtime.start << "ns" << std::endl;
        CHECK(max_error <= precision);
    }
}

TEST_CASE("test k0 accelerator with extracted fpga supercap data", "[k0_acc]")
{
    if (std::getenv("XCL_EMULATION_MODE") != nullptr)
    {
        test_k0_accelerator_with_extracted_data("supercap", "fpga", 0, 5, BEST_PRECISION);
    }
    else
    {
        test_k0_accelerator_with_extracted_data("supercap", "fpga", 0, -1, BEST_PRECISION);
    }
}

TEST_CASE("test k0 accelerator with extracted fpga graphene data", "[k0_acc]")
{
    if (std::getenv("XCL_EMULATION_MODE") != nullptr)
    {
        test_k0_accelerator_with_extracted_data("graphene", "fpga", 0, 5, BEST_PRECISION);
    }
    else
    {
        test_k0_accelerator_with_extracted_data("graphene", "fpga", 0, -1, BEST_PRECISION);
    }
}

TEST_CASE("Benchmark k0 accelerator on supercap", "[k0_acc][!benchmark]")
{
    Accelerator acc(Experiment("supercap"), AcceleratorDesign::K0_ACC, 0);
    KernelRuns runs("supercap", "k0", "fpga");
    runs.deserialize();
    std::vector<double, aligned_allocator<double>> kernel_out(runs.num);

    BENCHMARK("supercap k0 iteration 42")
    {
        acc.k0(runs.start, runs.end, runs.q[42].data(), kernel_out.data());
    };
}

TEST_CASE("Benchmark k0 accelerator on graphene", "[k0_acc][!benchmark]")
{
    Accelerator acc(Experiment("graphene"), AcceleratorDesign::K0_ACC, 0);
    KernelRuns runs("graphene", "k0", "fpga");
    runs.deserialize();
    std::vector<double, aligned_allocator<double>> kernel_out(runs.num);

    BENCHMARK("graphene k0 iteration 42")
    {
        acc.k0(runs.start, runs.end, runs.q[42].data(), kernel_out.data());
    };
}

void test_k0_accelerator_with_applied_distribution(const char *experiment_name, const char *target)
{
    // test 15 variations, primes up to and including maximum of 46
    std::vector<int> test_splits = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 46};

    Experiment experiment(experiment_name, vector_size);
    Accelerator acc(experiment, AcceleratorDesign::K0_ACC, 0);
    KernelRuns runs(experiment_name, "k0", target);
    runs.deserialize();

    std::vector<double, aligned_allocator<double>> single_kernel_out(experiment.num_pad, 0.0);
    std::vector<double, aligned_allocator<double>> multiple_kernel_out(experiment.num_pad, 0.0);
    std::vector<double, aligned_allocator<double>> multiple_kernel_out_part(experiment.num_pad, 0.0);

    std::vector<size_t> n_boundaries = Accelerator::get_boundaries(experiment.num, 1, 1, 0);

    acc.k0(n_boundaries[0], n_boundaries[1], runs.q[42].data(), single_kernel_out.data());

    for (unsigned int ite = 0; ite < test_splits.size(); ite++)
    {
        for (unsigned int elem = 0; elem < multiple_kernel_out.size(); elem++)
        {
            multiple_kernel_out[elem] = 0.0;
        }

        for (int k0_rank = 0; k0_rank < test_splits[ite]; k0_rank++)
        {
            for (unsigned int elem = 0; elem < multiple_kernel_out_part.size(); elem++)
            {
                multiple_kernel_out_part[elem] = 0.0;
            }

            n_boundaries = acc.get_boundaries(experiment.num, test_splits[ite], 1, k0_rank);
            acc.k0(n_boundaries[0], n_boundaries[1], runs.q[42].data(), multiple_kernel_out_part.data());

            for (unsigned int i = n_boundaries[0]; i < n_boundaries[1]; i++)
            {
                multiple_kernel_out[i] += multiple_kernel_out_part[i];
            }
        }

        for (int j = 0; j < experiment.num; j++)
        {
            CHECK(single_kernel_out[j] == multiple_kernel_out[j]);
        }
    }
}

TEST_CASE("Compare k0 single fpga results with multiple fpga version (testing 15 variations) on supercap",
          "[k0_acc][dist][supercap]")
{
    test_k0_accelerator_with_applied_distribution("supercap", "fpga");
}

TEST_CASE("Compare k0 single fpga results with multiple fpga version (testing 15 variations) on graphene",
          "[k0_acc][dist][graphene]")
{
    test_k0_accelerator_with_applied_distribution("graphene", "fpga");
}
