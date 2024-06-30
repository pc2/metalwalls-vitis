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
#include "Experiment.hpp"

std::string construct_path(std::string name)
{
    char path[100];
    snprintf(path, 100, "./experiments/%s/input_f90.dat", name.c_str());
    return std::string(path);
}

Experiment::Experiment()
{
}

void Experiment::postprocess()
{
    num2 = num - num1;
    alphasq = alpha * alpha;
    sqrpialpha = sqrt(M_PI) / alpha;
    k0PotFactor = 2.0 / (a * b);
    arec = 1.0 / a;
    brec = 1.0 / b;
    etasqrt2 = eta / sqrt(2.0);
    selfPotFactor = (2.0 * alpha - sqrt(2.0) * eta) / sqrt(M_PI);
    twopiarec = 2 * M_PI / a;
    twopibrec = 2 * M_PI / b;
    twopicrec = 2 * M_PI / c;
    alphaconst = -1.0 / (4.0 * alphasq);
    zweight = 2.0 * M_PI / c;
    lrPotFactor = 4.0 / ((a) * (b)) * zweight;
    numKModes = (2 * kzMax + 1) * (2 * kxMax * kyMax + kxMax + kyMax);

    for (int32_t i = 0; i < num; ++i)
    {
        xyz[4 * i] = x[i];
        xyz[4 * i + 1] = y[i];
        xyz[4 * i + 2] = z[i];
        xyz[4 * i + 3] = 0.0;
    }
}

void Experiment::resize()
{
    p.resize(num_pad, 0.0);
    x_cg.resize(num_pad, 0.0);
    b_cg.resize(num_pad, 0.0);
    res.resize(num_pad, 0.0);
    Ap.resize(num_pad, 0.0);
    x.resize(num_pad, 0.0);
    y.resize(num_pad, 0.0);
    z.resize(num_pad, 0.0);
    xyz.resize(4 * num_pad, 0.0);
}

Experiment::Experiment(std::string path_to_inputfile, int chunk_size)
{
    std::ifstream infile(path_to_inputfile);

    infile >> res_tol >> num >> num1 >> a >> b >> c >> alpha >> eta >> rcutsq >> rksqmax >> kxMax >> kyMax >> kzMax >>
        max_iterations;

    num_pad = num % chunk_size == 0 ? num : num + chunk_size - num % chunk_size;
    resize();

    for (int i = 0; i < num; i++)
    {
        infile >> b_cg[i] >> x_cg[i] >> res[i] >> p[i] >> Ap[i] >> x[i] >> y[i] >> z[i];
    }

    postprocess();
}

Experiment::Experiment(const char *name, int chunk_size) : Experiment(construct_path(name), chunk_size)
{
}

Experiment::Experiment(size_t seed, int num) : num(num)
{
    std::seed_seq s{seed};
    std::mt19937 generator(s);
    std::uniform_real_distribution<double> dist(0, 1.0);
    auto gen_random = [&]() { return dist(generator); };
    num_pad = num;
    num1 = gen_random() * num;
    kxMax = gen_random() * 100;
    kyMax = gen_random() * 100;
    kzMax = gen_random() * 100;
    max_iterations = 1000 + num;

    alpha = gen_random();
    eta = gen_random();
    a = gen_random();
    b = gen_random();
    c = gen_random();
    res_tol = gen_random();
    rcutsq = gen_random();
    rksqmax = gen_random();

    resize();

    for (int i = 0; i < num; i++)
    {
        b_cg[i] = gen_random();
        x_cg[i] = gen_random();
        res[i] = gen_random();
        p[i] = gen_random();
        Ap[i] = gen_random();
        x[i] = gen_random();
        y[i] = gen_random();
        z[i] = gen_random();
    }

    postprocess();
}

inline bool Experiment::unequal_vector(const std::vector<double, aligned_allocator<double>> &lhs,
                                       const std::vector<double, aligned_allocator<double>> &rhs) const
{
    if (lhs.size() != rhs.size())
    {
        return true;
    }
    for (size_t i = 0; i < lhs.size(); i++)
    {
        if (abs(lhs[i] - rhs[i]) > 1e-10)
        {
            return true;
        }
    }
    return false;
}

bool Experiment::operator==(const Experiment &rhs) const
{

    if (abs(num2 - rhs.num2) > 1e-10 || abs(alphasq - rhs.alphasq) > 1e-10 ||
        abs(sqrpialpha - rhs.sqrpialpha) > 1e-10 || abs(k0PotFactor - rhs.k0PotFactor) > 1e-10 ||
        abs(arec - rhs.arec) > 1e-10 || abs(brec - rhs.brec) > 1e-10 || abs(etasqrt2 - rhs.etasqrt2) > 1e-10 ||
        abs(selfPotFactor - rhs.selfPotFactor) > 1e-10 || abs(twopiarec - rhs.twopiarec) > 1e-10 ||
        abs(twopibrec - rhs.twopibrec) > 1e-10 || abs(twopicrec - rhs.twopicrec) > 1e-10 ||
        abs(alphaconst - rhs.alphaconst) > 1e-10 || abs(zweight - rhs.zweight) > 1e-10 ||
        abs(lrPotFactor - rhs.lrPotFactor) > 1e-10 || abs(numKModes - rhs.numKModes) > 1e-10)
    {
        return false;
    }

    if (unequal_vector(p, rhs.p) || unequal_vector(x_cg, rhs.x_cg) || unequal_vector(b_cg, rhs.b_cg) ||
        unequal_vector(res, rhs.res) || unequal_vector(Ap, rhs.Ap) || unequal_vector(x, rhs.x) ||
        unequal_vector(y, rhs.y) || unequal_vector(z, rhs.z))
    {

        return false;
    }

    return true;
}
