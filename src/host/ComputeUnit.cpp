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
#include "ComputeUnit.hpp"

ComputeUnit::ComputeUnit(std::string name, cl::Program program, cl::CommandQueue queue, cl::Context context)
    : context(context), queue(queue), program(program), name(name), kernel(), device_buffer_sizes(), device_buffers()
{
    cl_int err;
    OCL_CHECK(err, kernel = cl::Kernel(program, name.c_str(), &err));
}

cl::Event ComputeUnit::run(std::vector<cl::Event> wait_for_events)
{
    cl_int err;

    cl::Event event;
    OCL_CHECK(err, err = queue.enqueueTask(kernel, &wait_for_events, &event));
    queue.finish();
    return event;
}
