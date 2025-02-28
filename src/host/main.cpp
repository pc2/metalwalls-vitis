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
#include "../test/cpu_reference.cpp"
#include "../test/cpu_reference.hpp"
#include "Accelerator.hpp"
#include "Visualisation.hpp"
#include "json.hpp"
#include <mpi.h>

using json = nlohmann::json;

const int num_max = 53760;
const int num_max2 = 26880;

const int vector_size = 8;

#define DEBUG

double sum(const double *__restrict x, const int N_elem)
{
    double total_sum = 0.0;
    for (int i = 0; i < N_elem; i++)
    {
        total_sum += x[i];
    }
    return (total_sum);
}

double dot_product(const double *__restrict x, const double *__restrict y, const int N_elem)
{
    double product = 0.0;
    for (int i = 0; i < N_elem; i++)
    {
        product = product + x[i] * y[i];
    }
    return (product);
}

double calc_k0_time(double num_splits, double num)
{
    double unroll_factor = 8;
    double freq = 300.0e6;

    double n_read_cycles = 3 + std::ceil(num / unroll_factor);
    double n_compute_cycles = 166 + (num / num_splits) * (328.0 + std::ceil(num / unroll_factor));
    double n_write_cycles = 4 + std::ceil(num / unroll_factor);
    return (n_read_cycles + n_compute_cycles + n_write_cycles) / freq;
}

double calc_sr_time(double num_splits, double num)
{
    double unroll_factor = 4;
    double freq = 238.2e6;

    double n_read_cycles = 76 + num;
    double n_compute_cycles = 288 + (num / num_splits) * (464 + std::ceil((num / 2) / unroll_factor));
    return (n_read_cycles + n_compute_cycles) / freq;
}

double calc_lr_time(double num_splits, double num, double modes)
{
    double unroll_factor = 8;
    double freq = 190.1e6;

    double n_read_cycles = 3 + num;
    double n_compute_cycles = 264 + (modes / num_splits + 1) * (155 + std::ceil(num / unroll_factor) + 1);
    double n_write_cycles = 4 + num;
    return (n_read_cycles + n_compute_cycles + n_write_cycles) / freq;
}

double calculate_residual(Experiment &experiment, int iter, double *res, double *b_cg, double *Ap, double *p,
                          double *x_cg, double residual)
{
    if (iter == 0)
    {
        for (int i = 0; i < experiment.num; i++)
        {
            res[i] = b_cg[i] - Ap[i];
        }

        // p0 = z0
        for (int i = 0; i < experiment.num; i++)
        {
            p[i] = res[i];
        }

        return dot_product(res, res, experiment.num);
    }
    else
    {
        const double pAp = dot_product(p, Ap, experiment.num);
        double alpha_cg = residual / pAp;
        for (int i = 0; i < experiment.num; i++)
        {
            x_cg[i] = x_cg[i] + alpha_cg * p[i];
        }
        for (int i = 0; i < experiment.num; i++)
        {
            res[i] = res[i] - alpha_cg * Ap[i];
        }

        double residual_new = dot_product(res, res, experiment.num);

        /// Setup for next iteration
        double beta = residual_new / residual;
        for (int i = 0; i < experiment.num; i++)
        {
            p[i] = res[i] + beta * p[i];
        }
        return residual_new;
    }
}

class Metrics
{
  public:
    Metrics(Experiment &experiment) : experiment(experiment)
    {
    }

    void write()
    {
        std::ofstream metrics_file("metrics.json");
        metrics_file << metrics << std::endl;
    }

    void setup(int num_k0_units, int num_sr_units, int num_lr_units)
    {
        metrics["setup"] = {{"loadout", {{"k0", num_k0_units}, {"sr", num_sr_units}, {"lr", num_lr_units}}},
                            {"predicted_iter_time",
                             {{"k0", calc_k0_time(num_k0_units, experiment.num)},
                              {"sr", calc_sr_time(num_sr_units, experiment.num)},
                              {"lr", calc_lr_time(num_lr_units, experiment.num, experiment.numKModes)}}}};
        metrics["iteration"] = json::array();
    }

    void add_iteration(std::vector<double, aligned_allocator<double>> &accelerator_durations,
                       std::vector<double, aligned_allocator<double>> &accelerator_run_starts,
                       std::vector<double, aligned_allocator<double>> &accelerator_run_ends, double total_duration,
                       double residual_norm)
    {
        metrics["iteration"].push_back({{"rank_durations", accelerator_durations},
                                        {"run_starts", accelerator_run_starts},
                                        {"run_ends", accelerator_run_ends},
                                        {"total_duration", total_duration},
                                        {"residual_norm", residual_norm}});
    }

    void finalize(int iter, double runtime, double residual_norm)
    {
        metrics["final"] = {{"n_iterations", iter}, {"runtime", runtime}, {"residual_norm", residual_norm}};
    }

    json metrics;
    Experiment experiment;
};

int main(int argc, char **argv)
{
    int mpi_size;
    int mpi_rank = 0;

    if (MPI_Init(NULL, NULL) != MPI_SUCCESS)
    {
        std::cout << "Failed to initialize MPI\n";
        exit(-1);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (mpi_size < 3)
    {
        if (mpi_rank == 0)
        {
            std::cerr << "Illegal number of ranks. Metalwalls Vitis needs at least three ranks!" << std::endl;
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    MPI_Comm node_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, OMPI_COMM_TYPE_NODE, 0, MPI_INFO_NULL, &node_comm);

    int node_rank, node_size;
    MPI_Comm_rank(node_comm, &node_rank);
    MPI_Comm_size(node_comm, &node_size);

    std::string path_to_inputfile;
    std::string path_to_datafile;
    std::string path_to_runtimefile;

    if ((argc != 2) && (argc != 4))
    {
        if (mpi_rank == 0)
        {
            std::cerr << "Incorrect number of command line arguments given!" << std::endl;
            std::cerr << "Usage: " << argv[0] << " <input_f90.dat> [<data.inpt> <runtime.inpt>]" << std::endl;
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    bool emulation = std::getenv("XCL_EMULATION_MODE") != nullptr;

    bool host_cg = std::getenv("HOST_CG") != nullptr;

    bool visualise = (argc == 4) && host_cg;

    path_to_inputfile = argv[1];
    if (visualise)
    {
        path_to_datafile = argv[2];
        path_to_runtimefile = argv[3];
    }

    if (mpi_rank == 0)
    {
        printf("Launched with:\n");
        printf("Inputfile: %s\n", path_to_inputfile.c_str());
        printf("Visualisation %s\n", visualise ? "enabled" : "disabled");
        printf("CG Step on %s\n", host_cg ? "Host" : "Accelerator");
    }

    Experiment experiment(path_to_inputfile, vector_size);
    Metrics metrics(experiment);

    int num_k0_units = 0;
    int num_sr_units = 0;
    int num_lr_units = 0;
    double time_k0, time_sr, time_lr;
    double t_max = std::numeric_limits<double>::infinity();
    int acc_size = mpi_size - (host_cg ? 0 : 1);
    for (int i = 1; i <= acc_size - 2; i++)
    {
        time_k0 = calc_k0_time(i, experiment.num);
        for (int j = 1; j <= acc_size - i - 1; j++)
        {
            time_sr = calc_sr_time(j, experiment.num);
            int k = acc_size - i - j;
            time_lr = calc_lr_time(k, experiment.num, experiment.numKModes);
            if (t_max > fmax(time_k0, fmax(time_sr, time_lr)))
            {
                t_max = fmax(time_k0, fmax(time_sr, time_lr));
                if (mpi_rank == 0)
                    std::cout << "Time for (" << i << ", " << j << ", " << k
                              << ") = " << fmax(time_k0, fmax(time_sr, time_lr)) * 1000 << "ms" << std::endl;
                num_k0_units = i;
                num_sr_units = j;
                num_lr_units = k;
            };
        }
    }

    assert(acc_size == num_k0_units + num_sr_units + num_lr_units);

    metrics.setup(num_k0_units, num_sr_units, num_lr_units);
    if (mpi_rank == 0)
    {
        metrics.write();
    }

    Visualisation visualisation;
    if ((mpi_rank == 0) && visualise)
    {
        visualisation = Visualisation(experiment, path_to_datafile, path_to_runtimefile);
        visualisation.write_species();
    }

    std::vector<double, aligned_allocator<double>> accelerator_durations(mpi_size);
    std::vector<double, aligned_allocator<double>> accelerator_run_starts(mpi_size);
    std::vector<double, aligned_allocator<double>> accelerator_run_ends(mpi_size);

    AcceleratorDesign design;
    if (mpi_rank < num_k0_units)
    {
        design = AcceleratorDesign::K0_ACC;
    }
    else if (mpi_rank < num_k0_units + num_sr_units)
    {
        design = AcceleratorDesign::SR_ACC;
    }
    else if (mpi_rank < num_k0_units + num_sr_units + num_lr_units)
    {
        design = AcceleratorDesign::LR_ACC;
    }
    else
    {
        design = AcceleratorDesign::CG_ACC;
    }

    Accelerator acc = Accelerator(experiment, design, emulation ? 0 : node_rank);

    int iter = 0;
    double residual = 1.0;
    double residual_norm = 1.0;
    std::vector<double, aligned_allocator<double>> p(experiment.num_pad);
    std::vector<double, aligned_allocator<double>> x_cg(experiment.num_pad);
    std::vector<double, aligned_allocator<double>> b_cg(experiment.num_pad);
    std::vector<double, aligned_allocator<double>> res(experiment.num_pad);
    std::vector<double, aligned_allocator<double>> Ap(experiment.num_pad);

    for (int32_t i = 0; i < experiment.num_pad; ++i)
    {
        p[i] = experiment.p[i];
        x_cg[i] = experiment.x_cg[i];
        b_cg[i] = experiment.b_cg[i];
        res[i] = experiment.res[i];
        Ap[i] = experiment.Ap[i];
    }

    std::vector<size_t> boundaries;
    if (design == AcceleratorDesign::K0_ACC)
    {
        int k0_part_rank = mpi_rank;
        boundaries = acc.get_boundaries(experiment.num, num_k0_units, 1, k0_part_rank);
    }

    else if (design == AcceleratorDesign::SR_ACC)
    {
        int sr_part_rank = mpi_rank - num_k0_units;
        boundaries = acc.get_boundaries(experiment.num, num_sr_units, 1, sr_part_rank);
    }
    else if (design == AcceleratorDesign::LR_ACC)
    {
        int lr_part_rank = mpi_rank - num_k0_units - num_sr_units;
        boundaries = acc.get_boundaries(experiment.numKModes, num_lr_units, 1, lr_part_rank);
    }

    // Initial iteration
    auto iteration_start_instant = std::chrono::high_resolution_clock::now();

    if (design == AcceleratorDesign::K0_ACC)
    {
        acc.k0(boundaries[0], boundaries[1], p.data(), Ap.data());
    }
    else if (design == AcceleratorDesign::SR_ACC)
    {
        acc.sr(boundaries[0], boundaries[1], p.data(), Ap.data());
    }
    else if (design == AcceleratorDesign::LR_ACC)
    {
        acc.lr(boundaries[0], boundaries[1], p.data(), Ap.data());
    }

    if (host_cg && mpi_rank == 0)
    {
        for (int32_t i = 0; i < experiment.num; ++i)
        {
            Ap[i] -= experiment.selfPotFactor * p[i];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, Ap.data(), experiment.num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (!host_cg)
    {
        if (mpi_rank == (mpi_size - 1))
        {
            acc.cg(0, residual, b_cg.data(), p.data(), res.data(), x_cg.data(), p.data(), res.data(), &residual,
                   x_cg.data(), Ap.data());
        }
        MPI_Bcast(p.data(), experiment.num, MPI_DOUBLE, mpi_size - 1, MPI_COMM_WORLD);    // needed by k0/sr/lr/cg
        MPI_Bcast(res.data(), experiment.num, MPI_DOUBLE, mpi_size - 1, MPI_COMM_WORLD);  // needed by cg
        MPI_Bcast(&residual, 1, MPI_DOUBLE, mpi_size - 1, MPI_COMM_WORLD);                // needed to check convergence
        MPI_Bcast(x_cg.data(), experiment.num, MPI_DOUBLE, mpi_size - 1, MPI_COMM_WORLD); // needed by cg and as result
    }
    else
    {
        // Setup initial residual
        residual =
            calculate_residual(experiment, iter, res.data(), b_cg.data(), Ap.data(), p.data(), nullptr, residual);
    }

    if (mpi_rank == 0)
    {
        auto iteration_end_instant = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> iteration_duration = iteration_end_instant - iteration_start_instant;
        std::cout << "Initial Iteration: " << iteration_duration.count() * 1000.0 << " ms" << std::endl;
    }

    kernel_run_times_t kernel_run_time;

    auto start_instant = std::chrono::high_resolution_clock::now();

    /************************************** Start Conjugate Gradient (CG) ****************************/

    for (iter = 1; iter < experiment.max_iterations; ++iter)
    {
        if (mpi_rank == 0)
        {
            printf("we are starting iteration number %d with residual = %+4.15e, residual_norm =  %+4.15e \n", iter,
                   residual, residual_norm);
        }

        auto start_iteration_instant = std::chrono::high_resolution_clock::now();

        if (design == AcceleratorDesign::K0_ACC)
        {
            kernel_run_time = acc.k0(boundaries[0], boundaries[1], p.data(), Ap.data());
        }
        else if (design == AcceleratorDesign::SR_ACC)
        {
            kernel_run_time = acc.sr(boundaries[0], boundaries[1], p.data(), Ap.data());
        }
        else if (design == AcceleratorDesign::LR_ACC)
        {
            kernel_run_time = acc.lr(boundaries[0], boundaries[1], p.data(), Ap.data());
        }
        else if (design == AcceleratorDesign::CG_ACC)
        {
            for (int32_t i = 0; i < experiment.num; ++i)
            {
                Ap[i] = 0.0;
            }
        }

        auto end_accelerator_instant = std::chrono::high_resolution_clock::now();

        if (host_cg && mpi_rank == 0)
        {
            for (int32_t i = 0; i < experiment.num; ++i)
            {
                Ap[i] -= experiment.selfPotFactor * p[i];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, Ap.data(), experiment.num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (!host_cg)
        {
            if (mpi_rank == (mpi_size - 1))
            {
                kernel_run_time = acc.cg(iter, residual, b_cg.data(), p.data(), res.data(), x_cg.data(), p.data(),
                                         res.data(), &residual, x_cg.data(), Ap.data());
            }
            MPI_Bcast(p.data(), experiment.num, MPI_DOUBLE, mpi_size - 1, MPI_COMM_WORLD);   // needed by k0/sr/lr/cg
            MPI_Bcast(res.data(), experiment.num, MPI_DOUBLE, mpi_size - 1, MPI_COMM_WORLD); // needed by cg
            MPI_Bcast(&residual, 1, MPI_DOUBLE, mpi_size - 1, MPI_COMM_WORLD); // needed to check convergence
            MPI_Bcast(x_cg.data(), experiment.num, MPI_DOUBLE, mpi_size - 1,
                      MPI_COMM_WORLD); // needed by cg and as result

            residual_norm = sqrt(residual);
            if (residual_norm < experiment.res_tol)
            {
                break;
            }
        }
        else
        {
            residual = calculate_residual(experiment, iter, res.data(), b_cg.data(), Ap.data(), p.data(), x_cg.data(),
                                          residual);
            if ((mpi_rank == 0) && visualise)
            {
                visualisation.electrodes.set_x_cg(iter, x_cg.data());
            }

            residual_norm = sqrt(residual);
            if (residual_norm < experiment.res_tol)
            {
                break;
            }
        }

        auto end_iteration_instant = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> iteration_duration = end_iteration_instant - start_iteration_instant;

        std::chrono::duration<double> local_accelerator_duration = end_accelerator_instant - start_iteration_instant;
        double local_accelerator_seconds = local_accelerator_duration.count();
        MPI_Gather(&local_accelerator_seconds, 1, MPI_DOUBLE, accelerator_durations.data(), 1, MPI_DOUBLE, 0,
                   MPI_COMM_WORLD);
        MPI_Gather(&kernel_run_time.start, 1, MPI_DOUBLE, accelerator_run_starts.data(), 1, MPI_DOUBLE, 0,
                   MPI_COMM_WORLD);
        MPI_Gather(&kernel_run_time.end, 1, MPI_DOUBLE, accelerator_run_ends.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        metrics.add_iteration(accelerator_durations, accelerator_run_starts, accelerator_run_ends,
                              iteration_duration.count(), residual_norm);
        if (mpi_rank == 0)
        {
            metrics.write();
        }
    }

    if (mpi_rank == 0)
    {
        // Write Output
        std::ofstream outfile;
        outfile.precision(16);
        outfile.open("charges.out"); // TODO: change output directory or something
        outfile << "# Electrode atom charges in atomic unit\n";
        outfile << "# -------------------------------------\n";
        outfile << "# Charge:   1 au = 1.602176565e-19 C = 1 e\n";
        outfile << "# step      " << iter << "\n";
        for (int i = 0; i < experiment.num; i++)
        {
            outfile << x_cg[i] << "\n";
        }
        outfile.close();

        if (visualise)
        {
            visualisation.electrodes.write(iter);
        }

        // Report on overall time for CG
        auto end_instant = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end_instant - start_instant;
        std::cout << "Overall Time for CG:            " << duration.count() * 1000.0 << " ms";
        std::cout << "number of iterations: " << iter << std::endl;

        metrics.finalize(iter, duration.count(), residual_norm);
        if (mpi_rank == 0)
        {
            metrics.write();
        }
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
