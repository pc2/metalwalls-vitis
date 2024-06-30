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
#include <cmath>
#include <random>

/**
 * @brief Class holding all the constants of an experiment
 *
 */
class Experiment
{
  public:
    int num, num1, kxMax, kyMax, kzMax, max_iterations, num_pad;
    double alpha, eta, a, b, c, res_tol, rcutsq, rksqmax;
    int num2;
    double sqrpialpha;
    double k0PotFactor;
    double alphasq;
    double arec;
    double brec;
    double etasqrt2;
    double selfPotFactor;
    double twopiarec;
    double twopibrec;
    double twopicrec;
    double alphaconst;
    double zweight;
    double lrPotFactor;
    int numKModes;
    std::vector<double, aligned_allocator<double>> p;
    std::vector<double, aligned_allocator<double>> x_cg;
    std::vector<double, aligned_allocator<double>> b_cg;
    std::vector<double, aligned_allocator<double>> res;
    std::vector<double, aligned_allocator<double>> Ap;
    std::vector<double, aligned_allocator<double>> x;
    std::vector<double, aligned_allocator<double>> y;
    std::vector<double, aligned_allocator<double>> z;
    std::vector<double, aligned_allocator<double>> xyz;

    /**
     * @brief Calculates all constants which are derived from the input constants
     *
     * Used by constructor
     *
     */
    void postprocess();

    /**
     * @brief Resizes the vector to fit all the data
     *
     * Used by constructor
     *
     */
    void resize();

    /**
     * @brief Construct a new Experiment object and pads the data to chunk_size
     *
     * @param path_to_inputfile
     * @param chunk_size
     */
    Experiment(std::string path_to_inputfile, int chunk_size = 1);

    /**
     * @brief Construct a new Experiment object assuming a default path
     *
     * @param name Name of just the experiment
     * @param chunk_size
     */
    Experiment(const char *name, int chunk_size = 1);

    /**
     * @brief Construct a random experiment
     *
     * Can be used for testing.
     *
     * @param seed
     * @param num
     */
    Experiment(size_t seed, int num);

    /**
     * @brief Default constructor
     *
     */
    Experiment();

    /**
     * @brief Helper function for comparing two vectors of doubles
     *
     * Elements are considered equal if closer than 1e-10
     *
     * @param lhs
     * @param rhs
     * @return true if vectors are not equal
     */
    inline bool unequal_vector(const std::vector<double, aligned_allocator<double>> &lhs,
                               const std::vector<double, aligned_allocator<double>> &rhs) const;

    /**
     * @brief Compares two experiments by all individual values
     *
     * Elements are considered equal if closer than 1e-10
     *
     * @return true if Experiments are equal
     */
    bool operator==(const Experiment &rhs) const;
};
