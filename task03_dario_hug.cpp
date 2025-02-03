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
        std::cerr << "Unknown profile:" << profile;
        return;
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

// Function to solve the diffusion equation using RK2
void solve_diffusion_equation(std::vector<double> &f, int N, double dx, double dt, double diffusion_coefficient)
{
    std::vector<double> f_new(N, 0.0); // Temporary array to store updated values
    std::vector<double> k1(N, 0.0);
    std::vector<double> k2(N, 0.0);

    // Compute k1 (Euler step)
    for (int i = 0; i < N; ++i)
    {
        int i_upwind = (i == 0) ? N - 1 : i - 1;   // Periodic boundary for i-1
        int i_downwind = (i == N - 1) ? 0 : i + 1; // Periodic boundary for i+1

        double d2fdx2 = (f[i_downwind] - 2 * f[i] + f[i_upwind]) / (dx * dx); // Second spatial derivative
        k1[i] = dt * diffusion_coefficient * d2fdx2;                          // First step
    }

    // Compute k2 (use intermediate state)
    for (int i = 0; i < N; ++i)
    {
        int i_upwind = (i == 0) ? N - 1 : i - 1;   // Periodic boundary for i-1
        int i_downwind = (i == N - 1) ? 0 : i + 1; // Periodic boundary for i+1

        double d2fdx2 = ((f[i_downwind] + k1[i_downwind]) - 2 * (f[i] + k1[i]) + (f[i_upwind] + k1[i_upwind])) / (dx * dx);
        k2[i] = dt * diffusion_coefficient * d2fdx2; // Second step
    }

    // Update solution with RK2 formula
    for (int i = 0; i < N; ++i)
    {
        f_new[i] = f[i] + 0.5 * (k1[i] + k2[i]);
    }

    f = f_new; // Update f for the next step
}

// Function to solve the advection equation using RK II
void solve_advection_equation(std::vector<double> &f, int N, double dx, double dt, double v0)
{
    std::vector<double> f_new(N, 0.0); // Temporary array to store updated values

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

        // MinMod slope limiter
        double phi_right = std::max(0.0, std::min(r_right, 1.0));
        double phi_left = std::max(0.0, std::min(r_left, 1.0));

        // compute value at interpoints f_i+1/2 and f_i-1/2
        double f_right = f[i] + 0.5 * phi_right * (f[i_downwind] - f[i]);
        double f_left = f[i_upwind] + 0.5 * phi_left * (f[i] - f[i_upwind]);

        // Compute fluxes at both interfaces
        double F_right = v0 * f_right;
        double F_left = v0 * f_left;

        // Update f_new[i] using flux difference
        f_new[i] = f[i] - (dt / dx) * (F_right - F_left);
    }
    f = f_new; // Update f for the next step
}

int main(int argc, char *argv[])
{
    // A bit of defensive coding
    if (argc != 6)
    {
        std::cerr << "Usage: " << argv[0] << " <timesteps> <filename> <resolution> <starting profile> <diffusion_coefficient>" << std::endl;
        return 1;
    }

    int timesteps = std::atoi(argv[2]);
    if (timesteps <= 0)
    {
        std::cerr << "Error: timesteps must be a positive integer." << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    int resolution = std::atoi(argv[3]);
    std::string starting_profile = argv[4];
    double diffusion_coefficient = std::atof(argv[5]);

    // Parameters
    const double v0 = 1.0;
    const double x_min = 0.0, x_max = 1.0;
    int N = resolution;
    double dx = (x_max - x_min) / N;
    double dt = dx * dx / (4 * diffusion_coefficient); // dt based on D
    // double dt = dx / (2 * v0); // dt based on CFL rule

    std::vector<double> f(N, 0.0); // Initialize the starting condition grid

    // Use step or gaussian starting condition
    initialize(f, dx, starting_profile);

    // Store data for output as a vector
    std::vector<std::vector<double>> results;
    results.push_back(f);

    // Time integration loop
    for (int t = 0; t < timesteps; t++)
    {
        // Step 1: Solve diffusion equation for ∆t/2
        solve_diffusion_equation(f, N, dx, dt / 2, diffusion_coefficient);

        // Step 2: Solve advection equation for ∆t
        solve_advection_equation(f, N, dx, dt, v0);

        // Step 3: Solve diffusion equation for ∆t/2 again
        solve_diffusion_equation(f, N, dx, dt / 2, diffusion_coefficient);

        results.push_back(f); // Save current state
    }

    // Save results to a file
    save_to_csv(filename, results, dx, dt);

    std::string return_message = "Simulations completed, results saved in: ";
    std::cout << return_message << filename << std::endl;

    return 0;
}