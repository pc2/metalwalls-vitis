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
#include "../host/xcl2.hpp"
#include "../src/host/Accelerator.hpp"
#include "KernelRun.hpp"
#include "catch_amalgamated.hpp"
#include "cpu_reference.hpp"

#define BEST_PRECISION 1e-9

const int vector_size = 8;

void test_lr_accelerator_with_extracted_data(const char *experiment_name, const char *target, int start, int end,
                                             double precision)
{
    Accelerator acc(Experiment(experiment_name), AcceleratorDesign::LR_ACC, 0);

    KernelRuns runs = KernelRuns(experiment_name, "lr", target);
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

        std::vector<double, aligned_allocator<double>> out_fpga(runs.num);

        kernel_run_times_t cl_runtime = acc.lr(runs.start, runs.end, runs.q[i].data(), out_fpga.data());

        for (int j = 0; j < runs.num; j++)
        {
            double error = std::abs(out_fpga[j] - runs.result[i][j]);
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

TEST_CASE("test lr accelerator with extracted fpga supercap data", "[lr_acc]")
{
    if (std::getenv("XCL_EMULATION_MODE") != nullptr)
    {
        test_lr_accelerator_with_extracted_data("supercap", "fpga", 0, 5, BEST_PRECISION);
    }
    else
    {
        test_lr_accelerator_with_extracted_data("supercap", "fpga", 0, -1, BEST_PRECISION);
    }
}

TEST_CASE("test lr accelerator with extracted fpga graphene data", "[lr_acc]")
{
    if (std::getenv("XCL_EMULATION_MODE") != nullptr)
    {
        test_lr_accelerator_with_extracted_data("graphene", "fpga", 0, 5, BEST_PRECISION);
    }
    else
    {
        test_lr_accelerator_with_extracted_data("graphene", "fpga", 0, -1, BEST_PRECISION);
    }
}

TEST_CASE("Benchmark lr accelerator on supercap", "[lr_acc][!benchmark]")
{
    Accelerator acc(Experiment("supercap"), AcceleratorDesign::LR_ACC, 0);
    KernelRuns runs = KernelRuns("supercap", "lr", "fpga");
    runs.deserialize();
    std::vector<double, aligned_allocator<double>> out_fpga(runs.num);

    BENCHMARK("supercap lr iteration 20")
    {
        acc.lr(runs.start, runs.end, runs.q[42].data(), out_fpga.data());
    };
}

TEST_CASE("Benchmark lr accelerator on graphene", "[lr_acc][!benchmark]")
{
    Accelerator acc(Experiment("graphene"), AcceleratorDesign::LR_ACC, 0);
    KernelRuns runs = KernelRuns("graphene", "lr", "fpga");
    runs.deserialize();
    std::vector<double, aligned_allocator<double>> out_fpga(runs.num);

    BENCHMARK("graphene lr iteration 42")
    {
        acc.lr(runs.start, runs.end, runs.q[42].data(), out_fpga.data());
    };
}

void test_lr_accelerator_with_applied_distribution(const char *experiment_name, const char *target, double precision)
{
    // test 15 variations, primes up to and including maximum of 46
    std::vector<int> test_splits = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 46};
    const int i_iteration = 1;

    Experiment experiment(experiment_name, vector_size);
    Accelerator acc(experiment, AcceleratorDesign::LR_ACC, 0);
    KernelRuns runs(experiment_name, "lr", target);
    runs.deserialize();

    std::vector<double, aligned_allocator<double>> multiple_kernel_out(experiment.num_pad, 0.0);
    std::vector<double, aligned_allocator<double>> multiple_kernel_sum(experiment.num_pad, 0.0);

    for (unsigned int ite = 0; ite < test_splits.size(); ite++)
    {
        for (unsigned int elem = 0; elem < multiple_kernel_sum.size(); elem++)
        {
            multiple_kernel_sum[elem] = 0.0;
        }

        for (int lr_rank = 0; lr_rank < test_splits[ite]; lr_rank++)
        {
            for (unsigned int elem = 0; elem < multiple_kernel_out.size(); elem++)
            {
                multiple_kernel_out[elem] = 0.0;
            }

            std::vector<size_t> k_boundaries = acc.get_boundaries(experiment.numKModes, test_splits[ite], 1, lr_rank);

            acc.lr(k_boundaries[0], k_boundaries[1], runs.q[i_iteration].data(), multiple_kernel_out.data());

            for (int i = 0; i < experiment.num; i++)
            {
                multiple_kernel_sum[i] += multiple_kernel_out[i];
            }
        }

        double mean_error = 0.0;
        double min_error = std::numeric_limits<double>::infinity();
        double max_error = -std::numeric_limits<double>::infinity();

        for (int j = 0; j < experiment.num; j++)
        {
            double error = std::abs(runs.result[i_iteration][j] - multiple_kernel_sum[j]);
            mean_error += error / runs.num;
            min_error = std::min(min_error, error);
            max_error = std::max(max_error, error);
        }

        std::cout << "no. of ranks: " << test_splits[ite] << ", mean error: " << mean_error
                  << ", min_error: " << min_error << ", max_error: " << max_error << std::endl;
        CHECK(max_error <= precision);
    }
}

TEST_CASE("Compare lr single fpga results with multiple fpga version (testing 15 variations) on supercap",
          "[lr_acc][dist][supercap]")
{
    test_lr_accelerator_with_applied_distribution("supercap", "fpga", BEST_PRECISION);
}

TEST_CASE("Compare lr single fpga results with multiple fpga version (testing 15 variations) on graphene",
          "[lr_acc][dist][graphene]")
{
    test_lr_accelerator_with_applied_distribution("graphene", "fpga", BEST_PRECISION);
}
