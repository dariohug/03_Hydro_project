#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

// Define the initial profiles
void initialize(std::vector<double> &f, double dx, const std::string &profile)
{
    int N = f.size();
    // simple step profile
    if (profile == "step")
    {
        for (int i = 0; i < N; ++i)
        {
            double x = i * dx;
            if (x < 0.4)
                f[i] = 1.0;
            else if (x >= 0.4 && x < 0.6)
                f[i] = 2.0;
            else
                f[i] = 1.0;
        }
        // gaussian profile
    }
    else if (profile == "gaussian")
    {
        double sigma = 0.1;
        for (int i = 0; i < N; ++i)
        {
            double x = i * dx;
            f[i] = 1.0 + std::exp(-std::pow(x - 0.5, 2) / (2 * sigma * sigma));
        }
    }
    else
    {
        std::cerr << "Unknown profile" << std::endl;
    }
}

// Write data to a CSV file
void save_to_csv(const std::string &filename, const std::vector<std::vector<double>> &data, double dx, double dt)
{
    std::ofstream file(filename);
    // write all headers (x coordinate and every timestep)
    file << "x";
    for (size_t t = 0; t < data.size(); ++t)
        file << ",t=" << t * dt;
    file << "\n";

    // Write all data lines
    int N = data[0].size();
    for (int i = 0; i < N; ++i)
    {
        file << i * dx;
        for (size_t t = 0; t < data.size(); ++t)
            file << "," << data[t][i];
        file << "\n";
    }
    file.close();
}

int main(int argc, char *argv[])
{
    // A bit of defensive coding
    if (argc < 6)
    {
        std::cerr << "Usage: " << argv[0] << " <timesteps>" << "<filename>" << "<Resolution>" << "<starting profile>" << std::endl;
        return 1;
    }

    int timesteps = std::atoi(argv[1]);
    if (timesteps <= 0)
    {
        std::cerr << "Error: timesteps must be a positive integer." << std::endl;
        return 1;
    }
    int resolution = std::atoi(argv[3]);
    std::string filename = argv[2];
    std::string starting_profile = argv[4];
    std::string slope_lim_method = argv[5];

    std::cout << starting_profile << std::endl;

    // Parameters
    const double v0 = 1.0;
    const double x_min = 0.0, x_max = 1.0;
    int N = resolution;
    double dx = (x_max - x_min) / N;
    double dt = 0.5 * dx / v0; // given CFL condition

    std::vector<double> f(N, 0.0);     // Initialize the starting condition grid
    std::vector<double> f_new(N, 0.0); // Initialize the next timestep grid

    // use step starting condition or gaussian starting condition
    initialize(f, dx, starting_profile);

    // store data for output as a vector
    std::vector<std::vector<double>> results;
    results.push_back(f);

    // Time integration loop
    for (int t = 0; t < timesteps; ++t)
    {
        for (int i = 0; i < N; ++i)
        {

            // Periodic boundary conditions
            int i_upwind = (i == 0) ? N - 1 : i - 1;
            int i_downwind = (i == N - 1) ? 0 : i + 1;

            // Compute slope ratio r for slope limiter
            // Slope ratio
            double denominator_r_right = (f[i_downwind] - f[i] == 0.0) ? 1e-9 : f[i_downwind] - f[i];
            double r_right = (f[i] - f[i_upwind]) / denominator_r_right; // For F_{i+1/2} -> Ensure denominator != 0

            double denominator_r_left = (f[i] - f[i_upwind] == 0.0) ? 1e-9 : f[i] - f[i_upwind];
            double r_left = (f[i_upwind] - f[i]) / denominator_r_left; // For F_{i-1/2} -> Ensure denominator != 0

            double phi_right, phi_left;

            if (slope_lim_method == "vanleer")
            {
                // Van Leer slope limiter
                phi_right = (r_right + std::abs(r_right)) / (1.0 + std::abs(r_right));
                phi_left = (r_left + std::abs(r_left)) / (1.0 + std::abs(r_left));
            }
            else if (slope_lim_method == "minmod")
            {
                // MinMod slope limiter
                phi_right = std::max(0.0, std::min(r_right, 1.0));
                phi_left = std::max(0.0, std::min(r_left, 1.0));
            }
            else
            {
                std::cerr << "not a valid slope limiter method." << std::endl;
            }

            // compute value at interpoints f_i+1/2 and f_i-1/2
            double f_right = f[i] + 0.5 * phi_right * (f[i_downwind] - f[i]);
            double f_left = f[i_upwind] + 0.5 * phi_left * (f[i] - f[i_upwind]);

            // Compute fluxes at both interfaces
            double F_right = v0 * f_right;
            double F_left = v0 * f_left;

            // Update f_new[i] using flux difference
            f_new[i] = f[i] - (dt / dx) * (F_right - F_left);
        }

        f = f_new;            // Update the solution
        results.push_back(f); // Save current state
    }

    // Save results to a file
    save_to_csv(filename, results, dx, dt);

    std::string return_message = "Simulations completed, results saved in: ";
    std::cout << return_message << filename << std::endl;

    return 0;
}