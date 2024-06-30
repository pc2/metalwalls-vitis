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

#include "Experiment.hpp"

#include <fstream>
#include <sstream>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyVertex.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLDataSetWriter.h>

/**
 * @brief Helper class used for parsing the species
 *
 * The species are the model for all the molecules and atoms which are part of the overall experiment.
 *
 * Right now they are not needed for the overall application itself and are only used for visualisation.
 */
class Species
{
  public:
    /**
     * @brief Construct a new Species object
     *
     * The constructor assuemes the input streams are right before the relevant data and parses the information of the
     * species.
     *
     * Should not be used directly, as the Visualisation construcotr guarantess the correct usage.
     *
     * @param data Input stream of the data file parsed up to the Species
     * @param runtime Input stream of the runtime file parsed up to the species
     */
    Species(std::ifstream &data, std::ifstream &runtime);

    /**
     * @brief Writes a vtk file for this species
     *
     */
    void write();

    /**
     * @brief Prints the parsed data
     *
     * Only the first and the last coordinate are printed
     *
     */
    void print();

    std::string name;
    uint32_t count;
    double charge;
    std::vector<std::vector<double>> coordinates;
};

/**
 * @brief Representation of the electrodes
 *
 * This class is used for feeding the data of the charges for every timestep into it.
 *
 * The internal data is read from the already parsed Experiment.
 *
 * The write helper functions are used from the Visualisation class.
 *
 */
class Electrodes
{
  public:
    /**
     * @brief Construct a new Electrodes object
     *
     * @param experiment Already parsed Experiment
     */
    Electrodes(Experiment &experiment);

    /**
     * @brief Derault constructor
     *
     */
    Electrodes();

    /**
     * @brief Set the charges for one specific iteration
     *
     * @param iter CG iteration
     * @param x_cg pointer to the array of all charges for this iteration
     */
    void set_x_cg(int iter, double *x_cg);

    /**
     * @brief Writes the vtk file for one iteration
     *
     * @param iter CG iteration
     */
    void write_step(int iter);

    /**
     * @brief Writes the vtk files for all iterations
     *
     * @param end_iter Last iteration to write
     */
    void write(int end_iter);

    int num, max_iterations;
    std::vector<std::vector<double>> x_cg_data;
    std::vector<std::vector<double>> p_data;
    vtkNew<vtkPoints> points;
};

/**
 * @brief Class handling the visualisation
 *
 * Uses the data from the already parsed experiment and parses the needed data from the two input files.
 *
 * Creates a vtk file for every species, which contains the coordinates and the charge of each element.
 *
 * The data of the electrodes are written in a time series for every step of the CG iteration.
 *
 */
class Visualisation
{
  public:
    /**
     * @brief Construct a new Visualisation object
     *
     * @param experiment Reference to already parsed Experiment
     * @param path_to_datafile Path of the input file which contains the data
     * @param path_to_runtimefile Path of the input file which contains the runtime information
     */
    Visualisation(Experiment &experiment, std::string path_to_datafile, std::string path_to_runtimefile);

    /**
     * @brief Derault constructor
     *
     */
    Visualisation();

    /**
     * @brief Writes the vtk files
     *
     */
    void write_species();

    Electrodes electrodes;
    std::vector<Species> species;
};
