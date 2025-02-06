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
#include "../host/Experiment.hpp"
#include <filesystem>
#include <fstream>

#define ROOT_EXPERIMENT_DIR "experiments"
#define ROOT_DATA_DIR "/scratch/hpc-lco-kenter/papeg/mw/data/"

/**
 * @brief Reads an integer before a comma from a stringstream
 *
 * @param line
 * @return int
 */
inline int read_int(std::stringstream &line)
{
    std::string entry;
    std::getline(line, entry, ',');
    return stoi(entry);
}

/**
 * @brief Reads a double before a comma from a stringstream
 *
 * @param line
 * @return double
 */
inline double read_double(std::stringstream &line)
{
    std::string entry;
    std::getline(line, entry, ',');
    return stod(entry);
}

/**
 * \brief Base class for kernel runs
 *
 * A kernel run contains the input data for a certain kernel and the output
 * data it is supposed to produce. The usage pattern is that a kernel run
 * object is instantiated using the application or dataset name and the
 * iteration index. The constructor will then load the data associated with the
 * dataset and make it available in the fields of the subclass. To change the data
 * to another iteration an update method can be used.
 */
class KernelRuns
{
  public:
    /**
     * @brief The corresponding experiment object
     *
     */
    Experiment experiment;

    /**
     * \brief The name of the experiment that was the
     * source of this run.
     */
    std::string experiment_name;

    /**
     * \brief The name of the represented kernel.
     */
    std::string name;

    /**
     * \brief The execution mode of the kernel.
     */
    std::string mode;

    /**
     * \brief The path to the data store file
     */
    std::string path;

    /**
     * @brief The number of available atoms
     *
     *  can be different from experiment.num
     */
    int num;

    /**
     * \brief The total number of available iterations.
     */
    int iterations = 0;

    /**
     * @brief Scalar variable values
     *
     */
    int start, end, kx_start, ky_start, kz_start;

    /**
     * @brief Vector variable values
     *
     */
    std::vector<std::vector<double, aligned_allocator<double>>> q, result;

    /**
     * \brief Instantiate a new kernel run object.
     *
     * \param experiment Reference to the corresponding experiment object
     * \param experiment_name The name of the experiment that was the source of this run.
     * \param name The name of the represented kernel.
     * \param mode The execution mode of the kernel.
     */
    KernelRuns(Experiment experiment, std::string experiment_name, std::string name, std::string mode)
        : experiment(experiment), experiment_name(experiment_name), name(name), mode(mode), num(experiment.num)
    {
        char basepath_arr[1024];
        snprintf(basepath_arr, 1024, "%s/%s/%s/%s/", ROOT_DATA_DIR, experiment_name.c_str(), mode.c_str(),
                 name.c_str());

        char path_arr[1033];
        snprintf(path_arr, 1033, "%s/data.csv", basepath_arr);
        path = std::string(path_arr);

        q.resize(experiment.max_iterations);
        result.resize(experiment.max_iterations);
        for (int i = 0; i < experiment.max_iterations; i++)
        {
            q[i].resize(num);
            result[i].resize(num);
        }
    }

    std::string construct_inputfile(std::string &experiment_name)
    {
        char inputfile_arr[1024];
        snprintf(inputfile_arr, 1024, "%s/%s/input_f90.dat", ROOT_EXPERIMENT_DIR, experiment_name.c_str());
        return std::string(inputfile_arr);
    }

    KernelRuns(std::string experiment_name, std::string name, std::string mode)
        : KernelRuns(Experiment(construct_inputfile(experiment_name)), experiment_name, name, mode)
    {
    }

    KernelRuns(KernelRuns &copy) : KernelRuns(copy.experiment, copy.experiment_name, copy.name, copy.mode)
    {
    }

    /**
     * @brief Creates a random KernelRuns object
     *
     * @param experiment Reference to the corresponding experiment object
     * @param seed For generating the randomness
     * @param iterations How many random iterations to generate
     */
    KernelRuns(Experiment experiment, size_t seed, int gen_iterations) : KernelRuns(experiment, "test", "random", "cpu")
    {
        std::seed_seq s{seed};
        std::mt19937 generator(s);
        std::uniform_real_distribution<float> dist(0, 1.0);
        auto gen_random = [&]() { return dist(generator); };

        start = num * (gen_random() / 2);
        end = num * ((gen_random() / 2) + 0.5);
        kx_start = gen_random() * experiment.kxMax;
        ky_start = gen_random() * experiment.kyMax;
        kz_start = gen_random() * experiment.kzMax;
        iterations = gen_iterations;
        for (int i = 0; i < iterations; i++)
        {
            for (int n = 0; n < num; n++)
            {
                q[i][n] = gen_random();
                result[i][n] = gen_random();
            }
        }
    }

    void update_scalars(int num_u, int start_u, int end_u, int kx_start_u, int ky_start_u, int kz_start_u)
    {
        num = num_u;
        start = start_u;
        end = end_u;
        kx_start = kx_start_u;
        ky_start = ky_start_u;
        kz_start = kz_start_u;
    }

    /**
     * @brief Updates the variable data with the data of a new iteration
     *
     * @param iteration the new iteration
     */
    void update_iteration(int iteration, double *q_in, double *result_in)
    {
        iterations = iteration + 1;

        for (int n = 0; n < num; n++)
        {
            q[iteration][n] = q_in[n];
            result[iteration][n] = result_in[n];
        }
    }

    /**
     * \brief Load the dataset from disk, deserialize it, and store it in the
     * fields.
     */
    void deserialize()
    {
        iterations = 0;
        std::fstream fs(path, fs.in);
        std::string line;
        std::string entry;

        std::getline(fs, line);
        std::stringstream linestream(line);

        num = read_int(linestream);
        start = read_int(linestream);
        end = read_int(linestream);
        kx_start = read_int(linestream);
        ky_start = read_int(linestream);
        kz_start = read_int(linestream);
        while (std::getline(fs, line))
        {
            linestream = std::stringstream(line);
            for (int n = 0; n < num; n++)
            {
                q[iterations][n] = read_double(linestream);
            }
            std::getline(fs, line);
            linestream = std::stringstream(line);
            for (int n = 0; n < num; n++)
            {
                result[iterations][n] = read_double(linestream);
            }
            iterations++;
        }
    }

    /**
     * @brief Deserialize from custom path
     *
     */
    void deserialize(std::string custom_path)
    {
        std::string save_path = path;
        path = custom_path;
        deserialize();
        path = save_path;
    }

    /**
     * \brief Serialize the variable data of the iteration and write them to disk.
     */
    void serialize()
    {
        std::fstream fs(path, fs.out | fs.trunc);
        // write with full numeric precision
        fs << std::setprecision(std::numeric_limits<double>::digits10);
        fs << num << "," << start << "," << end << "," << kx_start << "," << ky_start << "," << kz_start << ","
           << "\n";
        for (int i = 0; i < iterations; i++)
        {
            for (int n = 0; n < num; n++)
            {
                fs << q[i][n] << ",";
            }
            fs << "\n";
            for (int n = 0; n < num; n++)
            {
                fs << result[i][n] << ",";
            }
            fs << "\n";
        }
    }

    /**
     * @brief Serialize from custom path
     *
     */
    void serialize(std::string custom_path)
    {
        std::string save_path = path;
        path = custom_path;
        serialize();
        path = save_path;
    }

    /**
     * \brief Compare the fields of the runs and test for elementwise equality.
     */
    inline bool operator==(const KernelRuns &rhs) const
    {
        if (!(experiment == rhs.experiment))
        {
            return false;
        }
        if (iterations != rhs.iterations)
        {
            return false;
        }
        if (q.size() != rhs.q.size())
        {
            return false;
        }
        if (result.size() != rhs.result.size())
        {
            return false;
        }
        for (int i = 0; i < iterations; i++)
        {

            for (int n = 0; n < num; n++)
            {
                if (abs(q[i][n] - rhs.q[i][n]) > 1e-10)
                {
                    return false;
                }
                if (abs(result[i][n] - rhs.result[i][n]) > 1e-10)
                {
                    return false;
                }
            }
        }
        return true;
    }
};
