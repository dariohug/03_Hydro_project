#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

// Define the initial conditions
void initialize(std::vector<double>& f, double dx, const std::string& profile) {
    int N = f.size();
    if (profile == "step") {
        for (int i = 0; i < N; ++i) {
            double x = i * dx;
            if (x < 0.4) f[i] = 1.0;
            else if (x >= 0.4 && x < 0.6) f[i] = 2.0;
            else f[i] = 1.0;
        }
    } else if (profile == "gaussian") {
        double sigma = 0.1;
        for (int i = 0; i < N; ++i) {
            double x = i * dx;
            f[i] = 1.0 + std::exp(-std::pow(x - 0.5, 2) / (2 * sigma * sigma));
        }
    }
}

// Write data to a CSV file
void save_to_csv(const std::string& filename, const std::vector<std::vector<double>>& data, double dx, double dt) {
    std::ofstream file(filename);
    file << "x";
    for (size_t t = 0; t < data.size(); ++t)
        file << ",t=" << t * dt;
    file << "\n";

    int N = data[0].size();
    for (int i = 0; i < N; ++i) {
        file << i * dx;
        for (size_t t = 0; t < data.size(); ++t)
            file << "," << data[t][i];
        file << "\n";
    }
    file.close();
}

int main() {

    std::cout << "Hello, World!" << std::endl;
    // Parameters
    const double v0 = 1.0;
    const double x_min = 0.0, x_max = 1.0;
    int N = 100;
    double dx = (x_max - x_min) / N;
    double dt = 0.5 * dx / v0; // CFL condition
    int timesteps = 500; // Number of timesteps to simulate

    // Initialize the grid and solution
    std::vector<double> f(N, 0.0), f_new(N, 0.0);
    initialize(f, dx, "gaussian");

    // Store data for output
    std::vector<std::vector<double>> results;
    results.push_back(f);

    // Time integration loop
    for (int t = 0; t < timesteps; ++t) {
        for (int i = 0; i < N; ++i) {
            // First-order upwind differencing
            int i_upwind = (i == 0) ? N - 1 : i - 1; // Periodic boundary
            f_new[i] = f[i] - v0 * dt / dx * (f[i] - f[i_upwind]);
        }

        f = f_new; // Update solution
        results.push_back(f); // Save current state
    }

    // Save results to a file
    save_to_csv("advection_results.csv", results, dx, dt);

    std::cout << "Simulation complete. Results saved to 'advection_results.csv'.\n";
    return 0;
}