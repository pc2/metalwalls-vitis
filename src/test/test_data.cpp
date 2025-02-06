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
#include "../host/Experiment.hpp"
#include "../host/xcl2.hpp"
#include "KernelRun.hpp"
#include "catch_amalgamated.hpp"
#include "cpu_reference.hpp"

TEST_CASE("read experiment", "[data]")
{
    Experiment supercap_experiment("supercap");
    Experiment random_experiment(42, 1024);
    std::cout << "created random experiment" << std::endl;
    REQUIRE(supercap_experiment == supercap_experiment);
    REQUIRE(random_experiment == random_experiment);
    REQUIRE_FALSE(supercap_experiment == random_experiment);

    Experiment another_random_experiment(43, 1024);
    REQUIRE(another_random_experiment == another_random_experiment);
    REQUIRE_FALSE(another_random_experiment == random_experiment);
    REQUIRE_FALSE(another_random_experiment == supercap_experiment);
}

TEST_CASE("random serialization", "[data]")
{
    Experiment random_experiment(42, 1024);
    KernelRuns random_runs_0(random_experiment, 7, 10);
    KernelRuns random_runs_1(random_experiment, 8, 10);

    REQUIRE(random_runs_0 == random_runs_0);
    REQUIRE(random_runs_1 == random_runs_1);

    REQUIRE_FALSE(random_runs_0 == random_runs_1);

    SKIP();
    random_runs_0.serialize("./random_data.csv");

    KernelRuns copy = KernelRuns(random_runs_0.experiment, "test", "random", "cpu");
    copy.deserialize("./random_data.csv");

    REQUIRE(random_runs_0 == copy);
}

void verify_k0_data(const char *experiment, const char *target, double precision)
{
    KernelRuns runs(experiment, "k0", target);
    runs.deserialize();
    std::vector<std::vector<double, aligned_allocator<double>>> results(runs.iterations);

#pragma omp parallel for shared(runs, results)
    for (int i = 0; i < runs.iterations; i++)
    {
        results[i] = std::vector<double, aligned_allocator<double>>(runs.experiment.num);
        double *result = results[i].data();

        MWPot_k0PotCPU(runs.experiment.a, runs.experiment.b, runs.experiment.alpha, 0, runs.experiment.num,
                       runs.experiment.z.data(), runs.q[i].data(), result);
    }

    for (int i = 0; i < runs.iterations; i++)
    {
        double *result = results[i].data();
        for (int j = 0; j < runs.experiment.num; j++)
        {
            REQUIRE_THAT(runs.result[i][j], Catch::Matchers::WithinRel(result[j], precision));
        }
    }
}

TEST_CASE("verify supercap fpga k0 data", "[data]")
{
    verify_k0_data("supercap", "fpga", 1e-10);
}

TEST_CASE("verify supercap femu k0 data", "[data]")
{
    verify_k0_data("supercap", "femu", 1e-8);
}

TEST_CASE("verify graphene fpga k0 data", "[data]")
{
    verify_k0_data("graphene", "fpga", 1e-10);
}

TEST_CASE("verify graphene femu k0 data", "[data]")
{
    verify_k0_data("graphene", "femu", 1e-8);
}

void verify_sr_data(int instance, const char *experiment, const char *target, double precision)
{
    char name_arr[16];
    snprintf(name_arr, 16, "sr%i", instance);
    KernelRuns runs(experiment, std::string(name_arr), target);
    runs.deserialize();
    std::vector<std::vector<double, aligned_allocator<double>>> results(runs.iterations);

    std::vector<double, aligned_allocator<double>> xyz(4 * runs.num);
    init_xyz(runs.num, xyz.data(), runs.experiment.x.data(), runs.experiment.y.data(), runs.experiment.z.data());

#pragma omp parallel for shared(runs, results, xyz)
    for (int i = 0; i < runs.iterations; i++)
    {
        results[i] = std::vector<double, aligned_allocator<double>>(runs.num);
        double *result = results[i].data();

        MWPot_srPotWallCPU(runs.experiment.a, runs.experiment.alpha, runs.experiment.b, runs.experiment.eta, runs.num,
                           runs.experiment.rcutsq, xyz.data(), runs.q[i].data(), result);
    }

    for (int i = 0; i < runs.iterations; i++)
    {
        double *result = results[i].data();
        for (int j = 0; j < runs.num; j++)
        {
            REQUIRE_THAT(runs.result[i][j], Catch::Matchers::WithinRel(result[j], precision));
        }
    }
}

TEST_CASE("verify supercap fpga sr1 data", "[data]")
{
    verify_sr_data(1, "supercap", "fpga", 1e-8);
}

TEST_CASE("verify supercap femu sr1 data", "[data]")
{
    verify_sr_data(1, "supercap", "femu", 1e-8);
}

TEST_CASE("verify supercap fpga sr2 data", "[data]")
{
    verify_sr_data(1, "supercap", "fpga", 1e-8);
}

TEST_CASE("verify supercap femu sr2 data", "[data]")
{
    verify_sr_data(1, "supercap", "femu", 1e-8);
}

TEST_CASE("verify graphene fpga sr1 data", "[data]")
{
    verify_sr_data(1, "graphene", "fpga", 1e-6);
}

TEST_CASE("verify graphene femu sr1 data", "[data]")
{
    verify_sr_data(1, "graphene", "femu", 1e-8);
}

TEST_CASE("verify graphene fpga sr2 data", "[data]")
{
    verify_sr_data(1, "graphene", "fpga", 1e-6);
}

TEST_CASE("verify graphene femu sr2 data", "[data]")
{
    verify_sr_data(1, "graphene", "femu", 1e-8);
}

void verify_lr_data(const char *experiment, const char *target, double precision)
{
    KernelRuns runs(experiment, "lr", target);
    runs.deserialize();
    std::vector<std::vector<double, aligned_allocator<double>>> results(runs.iterations);

    std::vector<double, aligned_allocator<double>> xyz(4 * runs.num);
    init_xyz(runs.num, xyz.data(), runs.experiment.x.data(), runs.experiment.y.data(), runs.experiment.z.data());

    std::vector<double, aligned_allocator<double>> cossinSkx(2 * runs.num * (runs.experiment.kxMax + 1));
    init_instream_cossinSk(runs.num, runs.experiment.kxMax, 0, runs.experiment.twopiarec, xyz.data(), cossinSkx.data());
    std::vector<double, aligned_allocator<double>> cossinSky(2 * runs.num * (runs.experiment.kyMax + 1));
    init_instream_cossinSk(runs.num, runs.experiment.kyMax, 1, runs.experiment.twopibrec, xyz.data(), cossinSky.data());
    std::vector<double, aligned_allocator<double>> cossinSkz(2 * runs.num * (runs.experiment.kzMax + 1));
    init_instream_cossinSk(runs.num, runs.experiment.kzMax, 2, runs.experiment.twopicrec, xyz.data(), cossinSkz.data());

#pragma omp parallel for shared(runs, results, xyz, cossinSkx, cossinSky, cossinSkz)
    for (int i = 0; i < runs.iterations; i++)
    {
        results[i] = std::vector<double, aligned_allocator<double>>(runs.num);
        double *result = results[i].data();

        MWPot_lrPotCPU(runs.experiment.a, runs.experiment.alpha, runs.experiment.b, runs.experiment.c,
                       runs.experiment.kxMax, runs.experiment.kyMax, runs.experiment.kzMax, 0, runs.num,
                       runs.experiment.rksqmax, cossinSkx.data(), cossinSky.data(), cossinSkz.data(), runs.q[i].data(),
                       result);
    }

    for (int i = 0; i < runs.iterations; i++)
    {
        double *result = results[i].data();
        for (int j = 0; j < runs.num; j++)
        {
            CHECK_THAT(runs.result[i][j], Catch::Matchers::WithinAbs(result[j], precision));
        }
    }
}

TEST_CASE("verify supercap fpga lr data", "[data]")
{
    verify_lr_data("supercap", "fpga", 1e-10);
}

TEST_CASE("verify supercap femu lr data", "[data]")
{
    verify_lr_data("supercap", "femu", 2e-13);
}
