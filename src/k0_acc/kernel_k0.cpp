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
#include <cmath>

#if AP_FLOAT
#define II_CYCLES 32
#else
#define II_CYCLES 20
#endif

#ifndef K0_INNER_U
#define K0_INNER_U 8
#endif

#ifndef K0_READ_U
#define K0_READ_U 8
#endif

const int num_max = 53760 + K0_INNER_U;

const int supercap_atoms = 2496;

extern "C"
{
    void kernel_k0(const int num_atoms, const int n_start, const int n_end, const double alpha, const double alphasq,
                   const double sqrpialpha, const double k0PotFactor, double *z, double *q, double *k0)
    {
#pragma HLS INTERFACE m_axi port = z offset = slave bundle = gmem1
#pragma HLS INTERFACE m_axi port = q offset = slave bundle = gmem2
#pragma HLS INTERFACE m_axi port = k0 offset = slave bundle = gmem3
        double z_loc[num_max];
        double q_loc[num_max];
        double k0_loc[num_max];
        const int write_end = n_end < num_atoms ? n_end : num_atoms;

    read:
        for (int word = 0; word < num_atoms; word += K0_READ_U)
        {
#pragma HLS loop_tripcount min = supercap_atoms / K0_READ_U max = supercap_atoms / K0_READ_U
#pragma HLS PIPELINE II = 1
            for (int element = 0; element < K0_READ_U; element++)
            {
#pragma HLS loop_tripcount min = K0_READ_U max = K0_READ_U
#pragma HLS UNROLL
                const int i = word + element;
                z_loc[i] = z[i];
                q_loc[i] = q[i];
            }
        }

    compute:
        for (int ii = n_start; ii < n_end; ii++)
        {
#pragma HLS loop_tripcount min = supercap_atoms max = supercap_atoms
            double s_reg[II_CYCLES];
            double minus_zi = -z_loc[ii];

        all_atoms_outer:
            for (int jb = 0; jb < num_atoms; jb += K0_INNER_U)
            {
#pragma HLS PIPELINE
#pragma HLS LOOP_FLATTEN OFF
#pragma HLS loop_tripcount min = supercap_atoms / K0_INNER_U max = supercap_atoms / K0_INNER_U
                double s_block = 0.0;
            all_atoms_inner:
                for (int jj = 0; jj < K0_INNER_U; jj++)
                {
#pragma HLS loop_tripcount min = K0_INNER_U max = K0_INNER_U
#pragma HLS UNROLL
                    if (jb + jj < num_atoms)
                    {
                        double zij = z_loc[jb + jj] + minus_zi;
                        double zijsq = zij * zij;
                        double minus_zijsq = -zijsq;
                        double expval = exp(minus_zijsq * alphasq);
                        double erfval = erf(zij * alpha);
                        double c = q_loc[jb + jj] * (sqrpialpha * expval + M_PI * zij * erfval);
                        s_block += c;
                    }
                }
                // feedback to right end of shift register
                s_reg[II_CYCLES - 1] = s_block + (jb >= (II_CYCLES - 1) * K0_INNER_U ? s_reg[0] : 0.0);

            shift:
                for (int s = 0; s < II_CYCLES - 1; s++)
                {
                    // defining shift register
                    s_reg[s] = s_reg[s + 1];
                }
            }
            double s_sum = 0.0;

        sum_shift_register:
            for (int s = 0; s < II_CYCLES - 1; ++s)
            {
#pragma HLS UNROLL
                s_sum += s_reg[s];
            }

            k0_loc[ii] = -(s_sum * k0PotFactor);
        }

        for (int ii = (n_start / K0_READ_U) * K0_READ_U; ii < ((write_end + K0_READ_U - 1) / K0_READ_U) * K0_READ_U;
             ii++)
        {
#pragma HLS loop_tripcount min = supercap_atoms / K0_INNER_U max = supercap_atoms / K0_INNER_U
            k0[ii] = k0_loc[ii];
        }

        for (int ii = write_end; ii < n_end; ii++)
        {

            k0[ii] = 0.0;
        }
    } // end of the kernel
} // end of the extern C
