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
#include "Visualisation.hpp"

Species::Species(std::ifstream &data, std::ifstream &runtime)
{
    std::string ignore;
    runtime >> ignore >> name;
    std::getline(runtime, ignore);
    runtime >> ignore >> count;
    std::getline(runtime, ignore);
    for (int i = 0; i < 3; i++)
    {
        std::string line;
        std::getline(runtime, line);
        std::stringstream lines(line);
        if (line.find("charge") != std::string::npos)
        {
            std::string word;
            lines >> word >> word;
            if (word == "neutral")
            {
                charge = 0.0;
                break;
            }
            else if (word == "point" || word == "gaussian")
            {
                lines >> word;
                charge = stod(word);
                break;
            }
            else
            {
                charge = stod(word);
                break;
            }
        }
    }
    coordinates.resize(count);
    for (uint32_t i = 0; i < count; i++)
    {
        coordinates[i].resize(3);
        std::string data_name;
        data >> data_name;
        if (data_name != name)
        {
            std::cerr << "parsing the wrong coordinate: " << data_name << ", expecting: " << name << std::endl;
            return;
        }
        data >> coordinates[i][0] >> coordinates[i][1] >> coordinates[i][2];
    }
}

void Species::write()
{
    std::stringstream filename;
    filename << "./visualisation/" << name << ".vtu" << std::endl;

    vtkNew<vtkUnstructuredGrid> grid;

    vtkNew<vtkDoubleArray> scalars;
    scalars->SetName("charge");

    vtkNew<vtkPoints> points;
    for (uint32_t i = 0; i < count; i++)
    {
        points->InsertNextPoint(coordinates[i][0], coordinates[i][1], coordinates[i][2]);
        scalars->InsertNextValue(charge);
    }

    grid->SetPoints(points);
    grid->GetPointData()->AddArray(scalars);

    vtkNew<vtkXMLDataSetWriter> writer;
    writer->SetFileName(filename.str().c_str());
    writer->SetInputData(grid);
    writer->Write();
}

void Species::print()
{
    std::cout << "name: " << name << std::endl;
    std::cout << "count: " << count << std::endl;
    std::cout << "charge: " << charge << std::endl;
    std::cout << "first coords: " << coordinates[0][0] << " " << coordinates[0][1] << " " << coordinates[0][2]
              << std::endl;
    std::cout << "last coords: " << coordinates[count - 1][0] << " " << coordinates[count - 1][1] << " "
              << coordinates[count - 1][2] << std::endl;
}

Electrodes::Electrodes(Experiment &experiment) : num(experiment.num), max_iterations(experiment.max_iterations)
{
    x_cg_data.resize(max_iterations);
    for (int i = 0; i < max_iterations; i++)
    {
        x_cg_data[i].resize(num);
    }

    for (int i = 0; i < num; i++)
    {
        points->InsertNextPoint(experiment.x[i], experiment.y[i], experiment.z[i]);
        x_cg_data[0][i] = experiment.x_cg[i];
    }
}

void Electrodes::set_x_cg(int iter, double *x_cg)
{
    for (int i = 0; i < num; i++)
    {
        x_cg_data[iter][i] = x_cg[i];
    }
}

void Electrodes::write_step(int iter)
{
    std::stringstream filename;
    filename << "./visualisation/electrodes." << std::setfill('0') << std::setw(3) << iter << ".vtu" << std::endl;

    vtkNew<vtkUnstructuredGrid> grid;
    grid->SetPoints(points);

    vtkNew<vtkDoubleArray> x_cg_scalars;
    x_cg_scalars->SetName("charge");

    for (int i = 0; i < num; i++)
    {
        x_cg_scalars->InsertNextValue(x_cg_data[iter][i]);
    }

    grid->GetPointData()->AddArray(x_cg_scalars);

    vtkNew<vtkXMLDataSetWriter> writer;
    writer->SetFileName(filename.str().c_str());
    writer->SetInputData(grid);
    writer->Write();
}

void Electrodes::write(int end_iter)
{
    if (end_iter < max_iterations)
    {
        // also write the current iteration
        end_iter++;
    }
    for (int i = 0; i < end_iter; i++)
    {
        write_step(i);
    }
}

Electrodes::Electrodes()
{
}

Visualisation::Visualisation(Experiment &experiment, std::string path_to_datafile, std::string path_to_runtimefile)
    : electrodes(experiment)
{
    std::ifstream datastream(path_to_datafile);
    std::ifstream runtimestream(path_to_runtimefile);
    // read to start of coordinates
    for (std::string line; line.find("# coordinates") == std::string::npos; std::getline(datastream, line))
    {
    }

    for (std::string line; std::getline(runtimestream, line);)
    {
        if (line.find("species_type") != std::string::npos)
        {
            Species parsed_species(datastream, runtimestream);
            // dont save electrode data, as we already have it
            if (parsed_species.name != "C1" && parsed_species.name != "C2")
            {
                species.push_back(parsed_species);
            }
        }
    }
}

Visualisation::Visualisation()
{
}

void Visualisation::write_species()
{
    for (auto &elem : species)
    {
        elem.write();
    }
}
