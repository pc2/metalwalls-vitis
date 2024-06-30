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

void test_distribution_algorithm(const char *experiment_name, const int num_type_units)
{
    Experiment experiment(experiment_name);
    int target_num_atoms = experiment.num;
    int target_numKModes = experiment.numKModes;

    std::vector<size_t> n_boundaries;
    std::vector<size_t> mode_boundaries;

    int prev_n_end = 0;
    int cur_n_start = 0;
    int cur_n_end = 0;

    int prev_mode_end = 0;
    int cur_mode_start = 0;
    int cur_mode_end = 0;

    for (int type_rank = 0; type_rank < num_type_units; type_rank++)
    {
        n_boundaries = Accelerator::get_boundaries(experiment.num, num_type_units, 1, type_rank);
        cur_n_start = n_boundaries[0];
        cur_n_end = n_boundaries[1];
        CHECK_THAT(prev_n_end, Catch::Matchers::WithinAbs(cur_n_start, 1));
        prev_n_end = cur_n_end;

        mode_boundaries = Accelerator::get_boundaries(experiment.numKModes, num_type_units, 1, type_rank);
        cur_mode_start = mode_boundaries[0];
        cur_mode_end = mode_boundaries[1];
        CHECK(prev_mode_end == cur_mode_start);
        prev_mode_end = cur_mode_end;
    }
    CHECK(prev_n_end == target_num_atoms);
    CHECK(prev_mode_end == target_numKModes);
}

TEST_CASE("test distribution algorithm for supercap for up to 46 units", "[dist][supercap]")
{
    for (int i = 0; i < 46; i++)
    {
        test_distribution_algorithm("supercap", i + 1);
    }
}

TEST_CASE("test distribution algorithm for graphene for up to 46 units", "[dist][graphene]")
{
    for (int i = 0; i < 46; i++)
    {
        test_distribution_algorithm("graphene", i + 1);
    }
}
