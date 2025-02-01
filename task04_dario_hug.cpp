#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

// Define the initial profiles
void initialize(std::vector<std::vector<double>> &f, double dx, const std::string &profile)
{
    int N = f.size();
    // simple step profile
    if (profile == "step")
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                double x = i * dx;
                double y = j * dx;
                if (x < 0.4 && y < 0.4)
                    f[i][j] = 1.0;
                else if (x >= 0.4 && x < 0.6 && y >= 0.4 && y < 0.6)
                    f[i][j] = 2.0;
                else
                    f[i][j] = 1.0;
            }
        }
        // gaussian profile
    }
    else if (profile == "gaussian")
    {
        double sigma = 0.1;
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                double x = i * dx;
                double y = j * dx;
                f[i][j] = 1.0 + std::exp(-(std::pow(x - 0.5, 2) + std::pow(y - 0.5, 2)) / (2 * sigma * sigma));
            }
        }
    }
    else
    {
        std::cerr << "Unknown profile" << std::endl;
    }
}

// Write data to a CSV file
void save_to_csv(const std::string &filename, const std::vector<std::vector<std::vector<double>>> &data, double dx, double dt)
{
    std::ofstream file(filename);
    // write all headers (x coordinate and every timestep)
    file << "x,y";
    for (size_t t = 0; t < data.size(); ++t)
        file << ",t=" << t * dt;
    file << "\n";

    // Write all data lines
    int N = data[0].size();
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            file << i * dx << "," << j * dx;
            for (size_t t = 0; t < data.size(); ++t)
                file << "," << data[t][i][j];
            file << "\n";
        }
    }
    file.close();
}

// Solve the two-dimensional advection problem for a single timestep
void solve_advection_equation(std::vector<std::vector<double>> &f, int N, double dx, double dt, double v0, std::string slope_lim_method)
{

    std::vector<std::vector<double>> f_new(N, std::vector<double>(N, 0.0)); // Temporary array to store updated values

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            // Periodic boundary conditions
            int i_upwind = (i == 0) ? N - 1 : i - 1;
            int i_downwind = (i == N - 1) ? 0 : i + 1;
            int j_upwind = (j == 0) ? N - 1 : j - 1;
            int j_downwind = (j == N - 1) ? 0 : j + 1;

            // Compute slope ratio r for Van Leer limiter
            double r_right = (f[i][j] - f[i_upwind][j]) / (f[i_downwind][j] - f[i][j] + 1e-8); // For F_{i+1/2} -> Ensure denominator != 0
            double r_left = (f[i_upwind][j] - f[i][j]) / (f[i][j] - f[i_upwind][j] + 1e-8);    // For F_{i-1/2} -> Ensure denominator != 0
            double r_top = (f[i][j] - f[i][j_upwind]) / (f[i][j_downwind] - f[i][j] + 1e-8);   // For F_{j+1/2} -> Ensure denominator != 0
            double r_bottom = (f[i][j_upwind] - f[i][j]) / (f[i][j] - f[i][j_upwind] + 1e-8);  // For F_{j-1/2} -> Ensure denominator != 0

            double phi_right, phi_left, phi_top, phi_bottom;

            if (slope_lim_method == "vanleer")
            {
                // Van Leer slope limiter
                phi_right = (r_right + std::abs(r_right)) / (1.0 + std::abs(r_right));
                phi_left = (r_left + std::abs(r_left)) / (1.0 + std::abs(r_left));
                phi_top = (r_top + std::abs(r_top)) / (1.0 + std::abs(r_top));
                phi_bottom = (r_bottom + std::abs(r_bottom)) / (1.0 + std::abs(r_bottom));
            }
            else if (slope_lim_method == "minmod")
            {
                // MinMod slope limiter
                phi_right = std::max(0.0, std::min(2 * r_right, std::min((1 + 2 * r_right) / 3.0, 2.0)));
                phi_left = std::max(0.0, std::min(2 * r_left, std::min((1 + 2 * r_left) / 3.0, 2.0)));
                phi_top = std::max(0.0, std::min(2 * r_top, std::min((1 + 2 * r_top) / 3.0, 2.0)));
                phi_bottom = std::max(0.0, std::min(2 * r_bottom, std::min((1 + 2 * r_bottom) / 3.0, 2.0)));
            }
            else
            {
                std::cerr << "not a valid slope limiter method." << std::endl;
            }

            // Reconstruct left and right states for box boundary F(i+1/2)
            double f_L_right = f[i][j] + 0.5 * phi_right * (f[i_downwind][j] - f[i][j]);
            double f_R_right = f[i_downwind][j] - 0.5 * phi_right * (f[i_downwind][j] - f[i][j]);

            // Reconstruct left and right states for box boundary F_{i-1/2}
            double f_L_left = f[i_upwind][j] + 0.5 * phi_left * (f[i][j] - f[i_upwind][j]);
            double f_R_left = f[i][j] - 0.5 * phi_left * (f[i][j] - f[i_upwind][j]);

            // Reconstruct top and bottom states for box boundary F(j+1/2)
            double f_L_top = f[i][j] + 0.5 * phi_top * (f[i][j_downwind] - f[i][j]);
            double f_R_top = f[i][j_downwind] - 0.5 * phi_top * (f[i][j_downwind] - f[i][j]);

            // Reconstruct top and bottom states for box boundary F_{j-1/2}
            double f_L_bottom = f[i][j_upwind] + 0.5 * phi_bottom * (f[i][j] - f[i][j_upwind]);
            double f_R_bottom = f[i][j] - 0.5 * phi_bottom * (f[i][j] - f[i][j_upwind]);

            // Compute fluxes at both interfaces using Riemann Solver
            double F_right = 0.5 * (v0 * f_L_right + v0 * f_R_right) - 0.5 * v0 * (f_R_right - f_L_right);
            double F_left = 0.5 * (v0 * f_L_left + v0 * f_R_left) - 0.5 * v0 * (f_R_left - f_L_left);
            double F_top = 0.5 * (v0 * f_L_top + v0 * f_R_top) - 0.5 * v0 * (f_R_top - f_L_top);
            double F_bottom = 0.5 * (v0 * f_L_bottom + v0 * f_R_bottom) - 0.5 * v0 * (f_R_bottom - f_L_bottom);

            // Update f_new[i][j] using flux difference
            f_new[i][j] = f[i][j] - (dt / dx) * (F_right - F_left + F_top - F_bottom);
        }
    }

    f = f_new;
}

int main(int argc, char *argv[])
{
    // A bit of defensive coding
    if (argc < 6)
    {
        std::cerr << "Usage: " << argv[0] << " <timesteps>" << "<filename>" << "<Resolution>" << "<starting profile>" << "slope_limiter_metohd" << std::endl;
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

    std::vector<std::vector<double>> f(N, std::vector<double>(N, 0.0)); // Initialize the starting condition grid

    // use step starting condition or gaussian starting condition
    initialize(f, dx, starting_profile);

    // store data for output as a vector
    std::vector<std::vector<std::vector<double>>> results;
    results.push_back(f);

    // Time integration loop
    for (int t = 0; t < timesteps; ++t)
    {
        solve_advection_equation(f, N, dx, dt, v0, slope_lim_method);
        results.push_back(f); // Save current state
    }

    // Save results to a file
    save_to_csv(filename, results, dx, dt);

    std::string return_message = "Simulations completed, results saved in: ";
    std::cout << return_message << filename << std::endl;

    return 0;
}