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
#include <cassert>
#include <map>
#include <string>
#include <vector>

/**
 * @brief Representation of a physical compute unit.
 *
 * This class represents and wraps an arbitrary compute unit on an FPGA that
 * has scalar arguments and is attached to possibly multiple HBM banks.
 *
 * Scalar arguments can be set using the @ref set_scalar method. It directly
 * calls @ref cl::Kernel::setArg, which means it does not check whether the
 * requested argument actually exists or whether the used type is correct.
 * Buffer arguments have to be set using the @ref write_buffer method. This
 * will allocate an OpenCL buffer with the required size, copy the data from
 * the host onto the FPGA, and set the requested kernel argument. Then, you can
 * start the compute unit using the @ref run method. After that, you can read
 * out the results of the computation by using the @ref read_buffer method.
 *
 * It is intentionally not possible to let a compute unit work on a buffer that
 * already exists since OpenCL buffers are supposed to represent HBM banks.
 * Letting two compute units simultanously work on the same HBM bank will
 * decrease their performance. Dependencies between reads, writes, and
 * computations also have to be handled manually using the return OpenCL events
 * and `wait_for_events` vectors. Lastly, you have to make sure that all
 * arguments were set at least once before running the compute unit.
 *
 * Compute unit is a term introduced by AMD to describe an instance of a
 * kernel. Following their definitions, a "kernel" appears to be the
 * description of a compute unit that has to be connected to memory or other
 * compute units, while compute units are the actual instances of these
 * kernels.
 */
class ComputeUnit
{
  public:
    /**
     * @brief Extract a compute unit from an FPGA program.
     *
     * @param name The name of the compute unit to extract.
     * @param program The program to extact the compute unit from.
     * @param queue The queue associated with the context and program.
     * @param context The context associated with the queue.
     */
    ComputeUnit(std::string name, cl::Program program, cl::CommandQueue queue, cl::Context context);

    /**
     * @brief Set a scalar parameter of the compute unit.
     *
     * Once set, the parameter retains its value until it's overriden.
     *
     * @tparam T The type of the scalar parameter to set.
     * @param index The parameter index.
     * @param param The parameter value.
     */
    template <typename T> void set_scalar(cl_int index, T param)
    {
        cl_int err;
        OCL_CHECK(err, err = kernel.setArg(index, param));

        // Remove any previously allocated buffers.
        device_buffer_sizes.erase(index);
        device_buffers.erase(index);
    }

    /**
     * @brief Allocate a new CL buffer, copy the data from the host to the
     * buffer, and pass it to the compute unit as an argument.
     *
     * This will only enqueue the copy operation and return immediately. Use the
     * returned event object to await the completion or handle dependencies.
     *
     * @tparam T The type of the buffer elements.
     * @param index The index of the parameter to set.
     * @param input_buffer A pointer to the host data to write. This has to be an
     * array with at least `length` elements.
     * @param length The number of elements to allocate on the device and to
     * copy.
     * @return The event of the copy operation.
     */
    template <typename T> cl::Event write_buffer(cl_int index, T const *input_buffer, size_t length)
    {
        cl_int err;

        size_t buffer_size = length * sizeof(T);

        if (device_buffer_sizes.count(index) == 0 || device_buffer_sizes[index] != buffer_size)
        {
            device_buffer_sizes[index] = buffer_size;

            OCL_CHECK(err, device_buffers[index] = cl::Buffer(context, CL_MEM_READ_WRITE, buffer_size, nullptr, &err));
            OCL_CHECK(err, err = kernel.setArg(index, device_buffers[index]));
        }

        cl::Event event;
        OCL_CHECK(err, err = queue.enqueueWriteBuffer(device_buffers[index], false, 0, buffer_size, input_buffer,
                                                      nullptr, &event));

        return event;
    }

    /**
     * @brief Copy the contents of a device buffer to the host.
     *
     * This will only enqueue the copy operation and return immediately. Use the
     * returned event object to await the completion or handle dependencies.
     *
     * @tparam T The type of the buffer elements.
     * @param index The index of the buffer to read.
     * @param output_buffer A pointer to the host buffer to write to. This has to
     * be an array with at least `length` elements.
     * @param length The number of elements in the device buffer.
     * @param wait_for_events A list of events that have to finish before the
     * read operation can start.
     * @return The event of the copy operation.
     */
    template <typename T>
    cl::Event read_buffer(cl_int index, T *output_buffer, size_t length, std::vector<cl::Event> wait_for_events)
    {
        size_t output_buffer_size = length * sizeof(T);
        assert(output_buffer_size == device_buffer_sizes.at(index));

        cl_int err;
        cl::Event event;

        OCL_CHECK(err, err = queue.enqueueReadBuffer(device_buffers[index], false, 0, device_buffer_sizes.at(index),
                                                     output_buffer, &wait_for_events, &event));

        return event;
    }

    /**
     * @brief Let the compute unit run.
     *
     * This will only enqueue the execution and return immediately. Use the
     * `wait_for_events` parameter and the return `cl::Event` object to await
     * completion or handle dependencies.
     *
     * @param wait_for_events A list of events that have to finish before the
     * computation can start.
     * @return The event of the computation.
     */
    cl::Event run(std::vector<cl::Event> wait_for_events);

  protected:
    cl::Context context;
    cl::CommandQueue queue;
    cl::Program program;
    std::string name;
    cl::Kernel kernel;
    std::map<size_t, size_t> device_buffer_sizes;
    std::map<size_t, cl::Buffer> device_buffers;
};
