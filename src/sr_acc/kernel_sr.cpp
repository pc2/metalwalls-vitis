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
#include "stdint.h"
#include <cmath>

#define SR_INNER_U 4
#define SR_READ_U 1
#define II_CYCLES 32

const int num_max = 53760 + SR_INNER_U;

const int supercap_trips = 2496;
const int trips = supercap_trips;

extern "C"
{
    void kernel_sr(const int n_start, const int n_end, double a, double alpha, double b, double eta, int32_t num,
                   double rcutsq, double *xyz, double *q, double *sr)
    {

#pragma HLS INTERFACE m_axi port = xyz offset = slave bundle = gmem1
#pragma HLS INTERFACE m_axi port = q offset = slave bundle = gmem2
#pragma HLS INTERFACE m_axi port = sr offset = slave bundle = gmem3

        double arec = 1.0 / a;
        double brec = 1.0 / b;
        double etasqrt2 = eta / sqrt(2.0);
        const int write_end = n_end < num ? n_end : num;

        // array declaration
        double x_loc[2][num_max / 2];
        double y_loc[2][num_max / 2];
        double z_loc[2][num_max / 2];
        double q_loc[2][num_max / 2];

    Read_inp:
        // we need to separate the xyz to 3 different vectors and read them in
        // parallel from 3 different HBM banks later
        for (int i_half = 0; i_half < 2; i_half++)
        {
#pragma HLS loop_tripcount min = 2 max = 2
            for (int i = 0; i < num / 2; i++)
            {
#pragma HLS loop_tripcount min = num_max / 2 max = num_max / 2
#pragma HLS loop_flatten
                int global_i = i + i_half * num / 2;
                x_loc[i_half][i] = xyz[4 * global_i];
                y_loc[i_half][i] = xyz[4 * global_i + 1];
                z_loc[i_half][i] = xyz[4 * global_i + 2];
                q_loc[i_half][i] = q[global_i];
            }
        }

    main_loop:
        for (int global_i = n_start; global_i < n_end; global_i++)
        {
#pragma HLS loop_tripcount min = num_max max = num_max
#pragma HLS pipeline off
            int i_half, i;

            if (global_i < num / 2)
            {
                i_half = 0;
                i = global_i;
            }
            else
            {
                i_half = 1;
                i = global_i - num / 2;
            }

            double xi = x_loc[i_half][i];
            double yi = y_loc[i_half][i];
            double zi = z_loc[i_half][i];

            // Shift register for accumulating values
            double s_reg[II_CYCLES];

            for (int i_reg = 0; i_reg < II_CYCLES; i_reg++)
            {
#pragma HLS unroll
                s_reg[i_reg] = 0.0;
            }

            // unroll the inner loop
        inner_loop:
            for (int jb = 0; jb < num / 2; jb += SR_INNER_U)
            {
#pragma HLS loop_tripcount min = num_max / SR_INNER_U / 2 max = num_max / SR_INNER_U / 2
#pragma HLS pipeline II = 1
#pragma HLS loop_flatten off
                double vi_block = 0.0;

            SR_INNER_U_loop:
                for (int jj = 0; jj < SR_INNER_U; jj++)
                {
#pragma HLS unroll
                    int j = jb + jj;

                    //  Compute minimum image distance with 2D periodic boundary
                    //  conditions
                    double xij = x_loc[i_half][j] - xi;
                    double xscale = floor(xij * arec + 0.5);
                    xij = xij - a * xscale;

                    double yij = y_loc[i_half][j] - yi;
                    double yscale = floor(yij * brec + 0.5);
                    yij = yij - b * yscale;

                    double zij = z_loc[i_half][j] - zi;

                    double drnorm2 = xij * xij + yij * yij + zij * zij;
                    bool cut = (drnorm2 < rcutsq) & (j != i);
                    double drnorm = sqrt(drnorm2);
                    // using standard erfc not the customized implemented
                    double erfdiff = (erfc(alpha * drnorm) - erfc(etasqrt2 * drnorm)) / drnorm;
                    double update = cut ? q_loc[i_half][j] * erfdiff : 0.0;
                    if (j < num / 2)
                    {
                        vi_block += update;
                    }
                }

                s_reg[II_CYCLES - 1] = vi_block + s_reg[0];
                // parallel shift
                for (char s = 0; s < II_CYCLES - 1; s++)
                {
#pragma HLS unroll
                    s_reg[s] = s_reg[s + 1];
                }
            }

            double s_sum = 0;
            // parallel reduction of partial sums from shift register
            for (char s = 0; s < II_CYCLES - 1; ++s)
            {
#pragma HLS unroll
                s_sum += s_reg[s];
            }
            sr[global_i] = s_sum;
        }
    }
}
