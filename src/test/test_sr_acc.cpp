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
#include "../host/xcl2.hpp"
#include "KernelRun.hpp"
#include "catch_amalgamated.hpp"
#include "cpu_reference.hpp"

#define BEST_PRECISION 1e-9

const int vector_size = 8;

void test_sr_accelerator_with_extracted_data(const char *experiment_name, const char *target, int iteration_start,
                                             int iteration_end, double precision)
{
    Experiment experiment(experiment_name);
    Accelerator acc(experiment, AcceleratorDesign::SR_ACC, 0);
    KernelRuns runs[2] = {KernelRuns(experiment_name, "sr1", target), KernelRuns(experiment_name, "sr2", target)};
    runs[0].deserialize();
    runs[1].deserialize();

    int num = experiment.num;

    if (iteration_end == -1)
    {
        iteration_end = runs[0].iterations;
    }

    for (int i = iteration_start; i < iteration_end; i++)
    {
        double mean_error = 0.0;
        double min_error = std::numeric_limits<double>::infinity();
        double max_error = -std::numeric_limits<double>::infinity();
        double l2_distance = 0.0;

        std::vector<double, aligned_allocator<double>> q(num);
        for (int i_run = 0; i_run < 2; i_run++)
        {
            for (int j = 0; j < num / 2; j++)
            {
                q[j + i_run * num / 2] = runs[i_run].q[i][j];
            }
        }

        std::vector<double, aligned_allocator<double>> kernel_out(num);

        kernel_run_times_t cl_runtime = acc.sr(0, num, q.data(), kernel_out.data());

        for (int i_run = 0; i_run < 2; i_run++)
        {
            for (int j = 0; j < num / 2; j++)
            {
                double error = std::abs(kernel_out[j + i_run * num / 2] - runs[i_run].result[i][j]);
                mean_error += error / num;
                min_error = std::min(min_error, error);
                max_error = std::max(max_error, error);
                l2_distance += error * error;
            }
        }
        l2_distance = std::sqrt(l2_distance);

        std::cout << "mean_error: " << mean_error << ", min_error: " << min_error << ", max_error: " << max_error
                  << ", l2_distance: " << l2_distance << std::endl;
        std::cout << "measured CL event runtime = " << cl_runtime.end - cl_runtime.start << "ns" << std::endl;
        CHECK(max_error <= precision);
    }
}

TEST_CASE("test sr accelerator with extracted fpga supercap data", "[sr_acc]")
{
    if (std::getenv("XCL_EMULATION_MODE") != nullptr)
    {
        test_sr_accelerator_with_extracted_data("supercap", "fpga", 0, 5, BEST_PRECISION);
    }
    else
    {
        test_sr_accelerator_with_extracted_data("supercap", "fpga", 0, -1, BEST_PRECISION);
    }
}

TEST_CASE("test sr accelerator with extracted fpga graphene data", "[sr_acc]")
{
    if (std::getenv("XCL_EMULATION_MODE") != nullptr)
    {
        test_sr_accelerator_with_extracted_data("graphene", "fpga", 0, 5, BEST_PRECISION);
    }
    else
    {
        test_sr_accelerator_with_extracted_data("graphene", "fpga", 0, -1, BEST_PRECISION);
    }
}

TEST_CASE("Benchmark sr accelerator on supercap", "[sr_acc][!benchmark]")
{
    Experiment experiment("supercap");
    Accelerator acc(experiment, AcceleratorDesign::SR_ACC, 0);
    KernelRuns runs[2] = {KernelRuns("supercap", "sr1", "fpga"), KernelRuns("supercap", "sr2", "fpga")};
    runs[0].deserialize();
    runs[1].deserialize();

    int num = experiment.num;

    std::vector<double, aligned_allocator<double>> q(num);
    for (int i_run = 0; i_run < 2; i_run++)
    {
        for (int j = 0; j < num / 2; j++)
        {
            q[j + i_run * num / 2] = runs[i_run].q[2][j];
        }
    }

    std::vector<double, aligned_allocator<double>> kernel_out(num);

    BENCHMARK("iteration 2")
    {
        acc.sr(0, num, q.data(), kernel_out.data());
    };
}

TEST_CASE("Benchmark sr accelerator on graphene", "[sr_acc][!benchmark]")
{
    Experiment experiment("graphene");
    Accelerator acc(experiment, AcceleratorDesign::SR_ACC, 0);
    KernelRuns runs[2] = {KernelRuns("graphene", "sr1", "fpga"), KernelRuns("graphene", "sr2", "fpga")};
    runs[0].deserialize();
    runs[1].deserialize();

    int num = experiment.num;

    std::vector<double, aligned_allocator<double>> q(num);
    for (int i_run = 0; i_run < 2; i_run++)
    {
        for (int j = 0; j < num / 2; j++)
        {
            q[j + i_run * num / 2] = runs[i_run].q[2][j];
        }
    }

    std::vector<double, aligned_allocator<double>> kernel_out(num);

    BENCHMARK("iteration 2")
    {
        acc.sr(0, num, q.data(), kernel_out.data());
    };
}

void test_sr_accelerator_with_applied_distribution(const char *experiment_name, const char *target)
{
    // test 15 variations, primes up to and including maximum of 46
    std::vector<int> test_splits = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 46};

    Experiment experiment(experiment_name);
    Accelerator acc(experiment, AcceleratorDesign::SR_ACC, 0);
    KernelRuns runs[2] = {KernelRuns(experiment_name, "sr1", target), KernelRuns(experiment_name, "sr2", target)};
    runs[0].deserialize();
    runs[1].deserialize();

    int num = experiment.num;

    std::vector<double, aligned_allocator<double>> q(num);
    for (int i_run = 0; i_run < 2; i_run++)
    {
        for (int j = 0; j < num / 2; j++)
        {
            q[j + i_run * num / 2] = runs[i_run].q[2][j];
        }
    }

    std::vector<double, aligned_allocator<double>> single_kernel_out(experiment.num_pad, 0.0);
    std::vector<double, aligned_allocator<double>> multiple_kernel_out(experiment.num_pad, 0.0);
    std::vector<double, aligned_allocator<double>> multiple_kernel_out_part(experiment.num_pad, 0.0);

    std::vector<size_t> n_boundaries = Accelerator::get_boundaries(experiment.num, 1, 1, 0);

    acc.sr(n_boundaries[0], n_boundaries[1], q.data(), single_kernel_out.data());

    for (unsigned int ite = 0; ite < test_splits.size(); ite++)
    {
        for (unsigned int elem = 0; elem < multiple_kernel_out.size(); elem++)
        {
            multiple_kernel_out[elem] = 0.0;
        }

        for (int sr_rank = 0; sr_rank < test_splits[ite]; sr_rank++)
        {
            for (unsigned int elem = 0; elem < multiple_kernel_out_part.size(); elem++)
            {
                multiple_kernel_out_part[elem] = 0.0;
            }

            n_boundaries = acc.get_boundaries(experiment.num, test_splits[ite], 1, sr_rank);
            acc.sr(n_boundaries[0], n_boundaries[1], q.data(), multiple_kernel_out_part.data());

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

TEST_CASE("Compare sr single fpga results with multiple fpga version (testing 15 variations) on supercap",
          "[sr_acc][dist][supercap]")
{
    test_sr_accelerator_with_applied_distribution("supercap", "fpga");
}

TEST_CASE("Compare sr single fpga results with multiple fpga version (testing 15 variations) on graphene",
          "[sr_acc][dist][graphene]")
{
    test_sr_accelerator_with_applied_distribution("graphene", "fpga");
}
