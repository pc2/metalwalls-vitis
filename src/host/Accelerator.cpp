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
#include "Accelerator.hpp"
#include <cassert>
#include <iomanip>
#include <sstream>

kernel_run_times_t execution_time_ms(cl::Event &exec_event)
{
    cl_int err;

    uint64_t nstimestart, nstimeend;
    OCL_CHECK(err, err = exec_event.getProfilingInfo<uint64_t>(CL_PROFILING_COMMAND_START, &nstimestart));
    OCL_CHECK(err, err = exec_event.getProfilingInfo<uint64_t>(CL_PROFILING_COMMAND_END, &nstimeend));

    double ntimestart_d = (double)(nstimestart) / 1000000;
    double ntimeend_d = (double)(nstimeend) / 1000000;
    return kernel_run_times_t{ntimestart_d, ntimeend_d};
}

inline std::string get_binary_file_path(AcceleratorDesign design)
{
    std::stringstream path_stream;

    char *raw_binary_dir = getenv("MW_BINARY_DIR");
    if (raw_binary_dir == nullptr)
    {
        path_stream << "./";
    }
    else
    {
        path_stream << raw_binary_dir << "/";
    }

    switch (design)
    {
    case AcceleratorDesign::K0_ACC:
        path_stream << "k0_acc_";
        break;
    case AcceleratorDesign::LR_ACC:
        path_stream << "lr_acc_";
        break;
    case AcceleratorDesign::SR_ACC:
        path_stream << "sr_acc_";
        break;
    case AcceleratorDesign::CG_ACC:
        path_stream << "cg_acc_";
        break;
    default:
        std::cerr << "Illegal design identifier!" << std::endl;
        exit(1);
        break;
    }

    char *raw_xcl_mode = getenv("XCL_EMULATION_MODE");
    if (raw_xcl_mode == nullptr)
    {
        path_stream << "hw";
    }
    else
    {
        std::string xcl_mode(raw_xcl_mode);
        if (xcl_mode == "sw_emu")
        {
            path_stream << "sw_emu";
        }
        else if (xcl_mode == "hw_emu")
        {
            path_stream << "hw_emu";
        }
        else
        {
            std::cerr << "Illegal emulation mode '" << raw_xcl_mode << "'!" << std::endl;
            exit(1);
        }
    }

    path_stream << ".xclbin";
    return path_stream.str();
}

Accelerator::Accelerator(Experiment experiment, AcceleratorDesign design, int unsigned device_index)
    : experiment(experiment), design(design), k0_cu(std::nullopt), sr_cu(std::nullopt), lr_cu(std::nullopt)
{
    cl_int err;

    auto devices = xcl::get_xil_devices();
    if (devices.size() <= device_index)
    {
        throw std::out_of_range("Not enough devices available");
    }
    auto device = devices[device_index];

    OCL_CHECK(err, context = cl::Context(device, nullptr, nullptr, nullptr, &err));
    OCL_CHECK(err, queue = cl::CommandQueue(context, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE,
                                            &err));

    auto binary_file_path = get_binary_file_path(design);
    auto fileBuf = xcl::read_binary_file(binary_file_path);
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};

    std::cout << "Programming device " << device_index << std::endl;
    OCL_CHECK(err, program = cl::Program(context, {device}, bins, nullptr, &err));

    if (design == AcceleratorDesign::K0_ACC)
    {
        k0_cu = ComputeUnit("kernel_k0", program, queue, context);

        k0_cu->set_scalar(0, experiment.num);
        k0_cu->set_scalar(3, experiment.alpha);
        k0_cu->set_scalar(4, experiment.alphasq);
        k0_cu->set_scalar(5, experiment.sqrpialpha);
        k0_cu->set_scalar(6, experiment.k0PotFactor);
        k0_cu->write_buffer(7, experiment.z.data(), experiment.num_pad).wait();
    }
    else if (design == AcceleratorDesign::SR_ACC)
    {
        sr_cu = ComputeUnit("kernel_sr", program, queue, context);

        sr_cu->set_scalar(2, experiment.a);
        sr_cu->set_scalar(3, experiment.alpha);
        sr_cu->set_scalar(4, experiment.b);
        sr_cu->set_scalar(5, experiment.eta);
        sr_cu->set_scalar(6, experiment.num);
        sr_cu->set_scalar(7, experiment.rcutsq);
        sr_cu->write_buffer(8, experiment.xyz.data(), experiment.num_pad * 4).wait();
    }
    else if (design == AcceleratorDesign::LR_ACC)
    {
        lr_cu = ComputeUnit("kernel_lr", program, queue, context);

        lr_cu->set_scalar(0, experiment.num);
        lr_cu->set_scalar(6, experiment.kyMax);
        lr_cu->set_scalar(7, experiment.kzMax);
        lr_cu->set_scalar(8, experiment.alphaconst);
        lr_cu->set_scalar(9, experiment.rksqmax);
        lr_cu->set_scalar(10, experiment.twopiarec);
        lr_cu->set_scalar(11, experiment.twopibrec);
        lr_cu->set_scalar(12, experiment.twopicrec);
        lr_cu->set_scalar(13, experiment.lrPotFactor);

        std::vector<cl::Event> write_events;

        write_events.push_back(lr_cu->write_buffer(14, experiment.x.data(), experiment.num_pad));
        write_events.push_back(lr_cu->write_buffer(15, experiment.y.data(), experiment.num_pad));
        write_events.push_back(lr_cu->write_buffer(16, experiment.z.data(), experiment.num_pad));

        for (auto &event : write_events)
        {
            event.wait();
        }
    }
    else if (design == AcceleratorDesign::CG_ACC)
    {
        if (!cg_cu.has_value())
        {
            cg_cu = ComputeUnit("kernel_cg", program, queue, context);
        }

        cg_cu->set_scalar(0, experiment.num);
        cg_cu->set_scalar(2, experiment.selfPotFactor);

        cg_cu->write_buffer(4, experiment.b_cg.data(), experiment.num).wait();
    }
    else
    {
        std::cerr << "unsupported AcceleratorDesign" << std::endl;
        exit(1);
    }
}

kernel_run_times_t Accelerator::k0(int n_start, int n_end, double const *q, double *k0)
{
    if (!k0_cu.has_value())
    {
        std::cerr << "k0 ComputeUnit uninitialized" << std::endl;
        exit(1);
    }

    std::vector<cl::Event> write_events;
    std::vector<cl::Event> execute_events;

    k0_cu->set_scalar(1, n_start);
    k0_cu->set_scalar(2, n_end);

    write_events.push_back(k0_cu->write_buffer(8, q, experiment.num_pad));
    write_events.push_back(k0_cu->write_buffer(9, k0, experiment.num_pad));

    execute_events.push_back(k0_cu->run(write_events));

    cl::Event read_event = k0_cu->read_buffer(9, k0, experiment.num_pad, execute_events);

    cl_int err;
    OCL_CHECK(err, err = read_event.wait());

    assert(execute_events.size() == 1);
    return execution_time_ms(execute_events[0]);
}

kernel_run_times_t Accelerator::sr(const int n_start, const int n_end, double const *q, double *sr)
{
    if (!sr_cu.has_value())
    {
        std::cerr << "sr ComputeUnit uninitialized" << std::endl;
        exit(1);
    }

    std::vector<cl::Event> write_events;
    std::vector<cl::Event> execute_events;

    sr_cu->set_scalar(0, n_start);
    sr_cu->set_scalar(1, n_end);

    write_events.push_back(sr_cu->write_buffer(9, q, experiment.num_pad));
    write_events.push_back(sr_cu->write_buffer(10, sr, experiment.num_pad));

    execute_events.push_back(sr_cu->run(write_events));

    cl::Event read_event = sr_cu->read_buffer(10, sr, experiment.num_pad, execute_events);

    cl_int err;
    OCL_CHECK(err, err = read_event.wait());

    assert(execute_events.size() == 1);
    return execution_time_ms(execute_events[0]);
}

kernel_run_times_t Accelerator::lr(int mode_start, int mode_end, double *q, double *lr)
{
    if (!lr_cu.has_value())
    {
        std::cerr << "lr ComputeUnit uninitialized" << std::endl;
        exit(1);
    }

    std::vector<cl::Event> write_events;
    std::vector<cl::Event> execute_events;

    int kx_start = 0, ky_start = 0, kz_start = 0;
    auto compute_xyz_start = [&]() {
        int mode = 0;
        for (kx_start = 0; kx_start <= experiment.kxMax; ++kx_start)
        {
            for (ky_start = (kx_start == 0) ? 1 : -experiment.kyMax; ky_start <= experiment.kyMax; ++ky_start)
            {
                for (kz_start = -experiment.kzMax; kz_start <= +experiment.kzMax; ++kz_start)
                {
                    if (mode == mode_start)
                        return;
                    else
                        mode++;
                }
            }
        }
    };
    compute_xyz_start();

    lr_cu->set_scalar(1, mode_start);
    lr_cu->set_scalar(2, mode_end);
    lr_cu->set_scalar(3, kx_start);
    lr_cu->set_scalar(4, ky_start);
    lr_cu->set_scalar(5, kz_start);

    write_events.push_back(lr_cu->write_buffer(17, q, experiment.num_pad));
    write_events.push_back(lr_cu->write_buffer(18, lr, experiment.num_pad));

    execute_events.push_back(lr_cu->run(write_events));

    cl::Event read_event = lr_cu->read_buffer(18, lr, experiment.num_pad, execute_events);

    cl_int err;
    OCL_CHECK(err, err = read_event.wait());

    assert(execute_events.size() == 1);
    return execution_time_ms(execute_events[0]);
}

kernel_run_times_t Accelerator::cg(const int iter, const double rsold, const double *b_cg, const double *q_in,
                                   const double *res_in, const double *x_cg_in, double *q_out, double *res_out,
                                   double *rsnew, double *x_cg_out, double *Ap)
{
    std::vector<cl::Event> write_events;
    std::vector<cl::Event> execute_events;
    std::vector<cl::Event> read_events;

    cg_cu->set_scalar(1, iter);
    cg_cu->set_scalar(3, rsold);

    write_events.push_back(cg_cu->write_buffer(5, q_in, experiment.num));
    write_events.push_back(cg_cu->write_buffer(6, res_in, experiment.num));
    write_events.push_back(cg_cu->write_buffer(7, x_cg_in, experiment.num));
    write_events.push_back(cg_cu->write_buffer(8, x_cg_out, experiment.num));

    write_events.push_back(cg_cu->write_buffer(9, q_out, experiment.num));
    write_events.push_back(cg_cu->write_buffer(10, res_out, experiment.num));
    write_events.push_back(cg_cu->write_buffer(11, rsnew, 1));
    write_events.push_back(cg_cu->write_buffer(12, Ap, experiment.num));

    execute_events.push_back(cg_cu->run(write_events));

    read_events.push_back(cg_cu->read_buffer(8, x_cg_out, experiment.num, execute_events));
    read_events.push_back(cg_cu->read_buffer(9, q_out, experiment.num, execute_events));
    read_events.push_back(cg_cu->read_buffer(10, res_out, experiment.num, execute_events));
    read_events.push_back(cg_cu->read_buffer(11, rsnew, 1, execute_events));

    cl_int err;
    for (auto &read_event : read_events)
    {
        OCL_CHECK(err, err = read_event.wait());
    }
    assert(execute_events.size() == 1);
    return execution_time_ms(execute_events[0]);
}

std::vector<size_t> Accelerator::get_boundaries(size_t workload, size_t number, size_t chunk_size, size_t rank)
{
    std::vector<size_t> boundaries;
    size_t lower_bound, upper_bound;
    size_t base_workload = ((workload / chunk_size) / number) * chunk_size;
    if ((workload / chunk_size) % number > rank)
    {
        base_workload += chunk_size;
        lower_bound = base_workload * rank;
        upper_bound = lower_bound + base_workload;
    }
    else
    {
        size_t num_bigger_workloads = (workload / chunk_size) % number;
        lower_bound = (base_workload + chunk_size) * num_bigger_workloads;
        lower_bound += (rank - num_bigger_workloads) * base_workload;
        upper_bound = lower_bound + base_workload;
    }
    boundaries.push_back(lower_bound);
    boundaries.push_back(upper_bound);
    return boundaries;
}

std::vector<size_t> Accelerator::get_distribution(size_t total_rows, size_t number, size_t chunk_size)
{
    std::vector<size_t> distribution;
    size_t base_rows = ((total_rows / chunk_size) / number) * chunk_size;
    size_t rest = ((total_rows / chunk_size) % number) * chunk_size;
    bool rest_flag = true;
    for (size_t i = 0; i < number; i++)
    {
        size_t rows_per_entry = base_rows;
        if (rest > chunk_size)
        {
            rows_per_entry += chunk_size;
            rest -= chunk_size;
        }
        else if (rest_flag && rest > 0)
        {
            rest_flag = false;
            rows_per_entry += rest;
        }
        if (rows_per_entry > 0)
        {
            distribution.push_back(rows_per_entry);
        }
    }
    return distribution;
}
