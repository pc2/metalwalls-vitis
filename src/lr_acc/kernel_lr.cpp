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
#include "hls_stream.h"
#include "hls_vector.h"
#include <math.h>
#include <stdint.h>

const int LR_U = 8;
const int II_CYCLES = 32;

const int num_max = 53760 + LR_U;

const int mode_const = 146598;

union LrPacket {
    std::array<std::array<double, LR_U>, 2> kxyz;
    std::array<double, 2> skread;
};

using LrStream = hls::stream<LrPacket>;

void compute(const int num_atoms, const int mode_start, const int mode_end, const int kx_start, const int ky_start,
             const int kz_start, double kyMax, double kzMax, double twopiarec, double twopibrec, double twopicrec,
             const double *x, const double *y, const double *z, const double *q, LrStream &out_stream)
{
    double x_local[num_max];
    double y_local[num_max];
    double z_local[num_max];
    double q_local[num_max];
#pragma HLS array_reshape variable = x_local type = cyclic factor = LR_U
#pragma HLS array_reshape variable = y_local type = cyclic factor = LR_U
#pragma HLS array_reshape variable = z_local type = cyclic factor = LR_U
#pragma HLS array_reshape variable = q_local type = cyclic factor = LR_U

read_loop:
    for (int ib = 0; ib < num_atoms; ib += LR_U)
    {
#pragma HLS loop_tripcount min = num_max / LR_U max = num_max / LR_U
#pragma HLS pipeline II = 1
        for (int ie = 0; ie < LR_U; ie++)
        {
#pragma HLS unroll
            int i = ib + ie;
            x_local[i] = x[i];
            y_local[i] = y[i];
            z_local[i] = z[i];
            q_local[i] = q[i];
        }
    }

    int kx = kx_start;
    int ky = ky_start;
    int kz = kz_start;

mode_loop:
    for (int imode = mode_start; imode < mode_end; ++imode)
    {
#pragma HLS loop_tripcount min = mode_const max = mode_const
        double rkx = kx * twopiarec;
        double rky = ky * twopibrec;
        double rkz = kz * twopicrec;

        double cosSk[II_CYCLES];
        double sinSk[II_CYCLES];

        for (int ii = 0; ii < II_CYCLES; ii++)
        {
#pragma HLS unroll
            cosSk[ii] = sinSk[ii] = 0.0;
        }

    atom_loop:
        for (int ib = 0; ib < num_atoms; ib += LR_U)
        {
#pragma HLS loop_tripcount min = num_max / LR_U max = num_max / LR_U

            double cosSkUpdate = 0.0;
            double sinSkUpdate = 0.0;

            LrPacket packet;
            for (int ie = 0; ie < LR_U; ie++)
            {
#pragma HLS unroll
                int i = ib + ie;

                double prodx = x_local[i] * rkx;
                double prody = y_local[i] * rky;
                double prodz = z_local[i] * rkz;
                double cosKx = cos(prodx);
                double sinKx = sin(prodx);
                double cosKy = cos(prody);
                double sinKy = sin(prody);
                double cosKz = cos(prodz);
                double sinKz = sin(prodz);

                double cosKxKy = cosKx * cosKy - sinKx * sinKy;
                double sinKxKy = sinKx * cosKy + cosKx * sinKy;
                double cosKxKyKz = cosKxKy * cosKz - sinKxKy * sinKz;
                double sinKxKyKz = sinKxKy * cosKz + cosKxKy * sinKz;

                packet.kxyz[0][ie] = cosKxKyKz;
                packet.kxyz[1][ie] = sinKxKyKz;

                if (i < num_atoms)
                {
                    cosSkUpdate += q_local[i] * cosKxKyKz;
                    sinSkUpdate += q_local[i] * sinKxKyKz;
                }
            }
            out_stream.write(packet);

            cosSk[II_CYCLES - 1] = cosSk[0] + cosSkUpdate;
            sinSk[II_CYCLES - 1] = sinSk[0] + sinSkUpdate;

            for (char s = 0; s < II_CYCLES - 1; s++)
            {
#pragma HLS unroll
                cosSk[s] = cosSk[s + 1];
                sinSk[s] = sinSk[s + 1];
            }
        }

        double cosSkread = 0.0;
        double sinSkread = 0.0;

        for (char s = 0; s < II_CYCLES - 1; ++s)
        {
#pragma HLS unroll
            cosSkread += cosSk[s];
            sinSkread += sinSk[s];
        }

        LrPacket packet;
        packet.skread = {cosSkread, sinSkread};
        out_stream.write(packet);

        kx = (kz == kzMax & ky == kyMax) ? ++kx : kx;
        ky = kz == kzMax ? ky == kyMax ? -kyMax : ++ky : ky;
        kz = kz == kzMax ? -kzMax : ++kz;
    }
}

void invert_stream(const int num_atoms, const int mode_start, const int mode_end, LrStream &in_stream,
                   LrStream &out_stream)
{
    LrPacket packet_buffer[2][num_max / LR_U + 1];
#pragma HLS array_partition variable = packet_buffer dim = 1 type = block factor = 2

    int n_packets_per_mode = num_atoms / LR_U;
    if (n_packets_per_mode % LR_U != 0)
    {
        n_packets_per_mode += 1;
    }

mode_loop:
    for (int imode = mode_start; imode < mode_end + 1; imode++)
    {
#pragma HLS loop_tripcount min = mode_const max = mode_const
    atom_loop:
        for (int i_packet = 0; i_packet < n_packets_per_mode; i_packet++)
        {
#pragma HLS loop_tripcount min = num_max / LR_U max = num_max / LR_U
            if (imode < mode_end)
            {
                packet_buffer[imode & 0b1][i_packet] = in_stream.read();
            }

            if (imode >= mode_start + 1)
            {
                out_stream.write(packet_buffer[(~imode) & 0b1][i_packet]);
            }
        }
        if (imode < mode_end)
        {
            out_stream.write(in_stream.read());
        }
    }
}

void accumulate(const int num_atoms, const int mode_start, const int mode_end, const int kx_start, const int ky_start,
                const int kz_start, double kyMax, double kzMax, double twopiarec, double twopibrec, double twopicrec,
                double rksqmax, double alphaconst, double lrPotFactor, LrStream &in_stream, double *lr)
{
    double lr_local[2][num_max];
#pragma HLS array_partition variable = lr_local dim = 1 type = block factor = 2
#pragma HLS array_partition variable = lr_local dim = 2 type = cyclic factor = LR_U

init_local_buffer:
    for (int ib = 0; ib < num_atoms; ib += LR_U)
    {
#pragma HLS loop_tripcount min = num_max / LR_U max = num_max / LR_U
#pragma HLS pipeline II = 1
        for (int ie = 0; ie < LR_U; ie++)
        {
#pragma HLS unroll
            int i = ib + ie;
            lr_local[0][i] = 0.0;
            lr_local[1][i] = 0.0;
        }
    }

    int kx = kx_start;
    int ky = ky_start;
    int kz = kz_start;

mode_loop:
    for (int imode = mode_start; imode < mode_end; imode++)
    {
#pragma HLS loop_tripcount min = mode_const max = mode_const
        double rkx = kx * twopiarec;
        double rky = ky * twopibrec;
        double rkz = kz * twopicrec;
        double rkk = rkx * rkx + rky * rky + rkz * rkz;
        double rkkValid = rkk < rksqmax;
        double alphaSk = rkkValid ? exp(alphaconst * rkk) / rkk : 0.0;

        LrPacket packet = in_stream.read();
        double cosSkread = packet.skread[0];
        double sinSkread = packet.skread[1];

    atom_loop:
        for (int ib = 0; ib < num_atoms; ib += LR_U)
        {
#pragma HLS loop_tripcount min = num_max / LR_U max = num_max / LR_U

            packet = in_stream.read();
            for (int ie = 0; ie < LR_U; ie++)
            {
#pragma HLS unroll
                int i = ib + ie;
                double coskxyzread = packet.kxyz[0][ie];
                double sinkxyzread = packet.kxyz[1][ie];

                if (i < num_atoms)
                {
                    lr_local[(~imode) & 0b1][i] =
                        lr_local[imode & 0b1][i] +
                        lrPotFactor * alphaSk * (coskxyzread * cosSkread + sinkxyzread * sinSkread);
                }
            }
        }

        kx = (kz == kzMax & ky == kyMax) ? ++kx : kx;
        ky = kz == kzMax ? ky == kyMax ? -kyMax : ++ky : ky;
        kz = kz == kzMax ? -kzMax : ++kz;
    }

write_results:
    for (int ib = 0; ib < num_atoms; ib += LR_U)
    {
#pragma HLS loop_tripcount min = num_max / LR_U max = num_max / LR_U
        for (int ie = 0; ie < LR_U; ie++)
        {
            int i = ib + ie;
            lr[i] = lr_local[mode_end & 0b1][i];
        }
    }
}

extern "C"
{
    void kernel_lr(const int num_atoms, const int mode_start, const int mode_end, const int kx_start,
                   const int ky_start, const int kz_start, const int kyMax, const int kzMax, const double alphaconst,
                   const double rksqmax, const double twopiarec, const double twopibrec, const double twopicrec,
                   const double lrPotFactor, const double *x, const double *y, const double *z, const double *q,
                   double *lr)
    {
#pragma HLS INTERFACE m_axi port = x offset = slave bundle = gmem1
#pragma HLS INTERFACE m_axi port = y offset = slave bundle = gmem2
#pragma HLS INTERFACE m_axi port = z offset = slave bundle = gmem3
#pragma HLS INTERFACE m_axi port = q offset = slave bundle = gmem4
#pragma HLS INTERFACE m_axi port = lr offset = slave bundle = gmem5
#pragma HLS dataflow

        LrStream stream0, stream1;

        compute(num_atoms, mode_start, mode_end, kx_start, ky_start, kz_start, kyMax, kzMax, twopiarec, twopibrec,
                twopicrec, x, y, z, q, stream0);
        invert_stream(num_atoms, mode_start, mode_end, stream0, stream1);
        accumulate(num_atoms, mode_start, mode_end, kx_start, ky_start, kz_start, kyMax, kzMax, twopiarec, twopibrec,
                   twopicrec, rksqmax, alphaconst, lrPotFactor, stream1, lr);
    }
}
