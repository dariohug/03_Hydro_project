#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

// Define the initial profiles
void initialize(std::vector<double>& f, double dx, const std::string& profile) {
    int N = f.size();
    // simple step profile
    if (profile == "step") {
        for (int i = 0; i < N; ++i) {
            double x = i * dx;
            if (x < 0.4) f[i] = 1.0;
            else if (x >= 0.4 && x < 0.6) f[i] = 2.0;
            else f[i] = 1.0;
        }
    // gaussian profile
    } else if (profile == "gaussian") {
        double sigma = 0.1;
        for (int i = 0; i < N; ++i) {
            double x = i * dx;
            f[i] = 1.0 + std::exp(-std::pow(x - 0.5, 2) / (2 * sigma * sigma));
        }
    } else {
        std::cerr << "Unknown profile" << std::endl;
    }
}

// Write data to a CSV file
void save_to_csv(const std::string& filename, const std::vector<std::vector<double>>& data, double dx, double dt) {
    std::ofstream file(filename);
    // write all headers (x coordinate and every timestep)
    file << "x";
    for (size_t t = 0; t < data.size(); ++t)
        file << ",t=" << t * dt;
    file << "\n";

    // Write all data lines
    int N = data[0].size();
    for (int i = 0; i < N; ++i) {
        file << i * dx;
        for (size_t t = 0; t < data.size(); ++t)
            file << "," << data[t][i];
        file << "\n";
    }
    file.close();
}

int main(int argc, char* argv[]) {
    // A bit of defensive coding
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <timesteps>" << "filename" << "starting profile" << std::endl;
        return 1;
    }
    
    int timesteps = std::atoi(argv[1]);
    if (timesteps <= 0) {
        std::cerr << "Error: timesteps must be a positive integer." << std::endl;
        return 1;
    }
    std::string filename = argv[2];
    std::string starting_profile = argv[3];

    // Parameters
    const double v0 = 1.0;
    const double x_min = 0.0, x_max = 1.0;
    int N = 100;
    double dx = (x_max - x_min) / N;
    double dt = 0.5 * dx / v0; // given CFL condition

    std::vector<double> f(N, 0.0); // Initialize the starting condition grid
    std::vector<double> f_new(N, 0.0); // Initialize the next timestep grid

    // use step starting condition or gaussian starting condition
    initialize(f, dx, starting_profile);

    // store data for output as a vector
    std::vector<std::vector<double>> results;
    results.push_back(f);

    // time integration loop
    for (int t = 0; t < timesteps; ++t) {
        for (int i = 0; i < N; ++i) {
            // first-order upwind differencing
            int i_upwind = (i == 0) ? N - 1 : i - 1; // be sure to exclude i == -1
            f_new[i] = f[i] - v0 * dt / dx * (f[i] - f[i_upwind]); // compute gridpoint at next timestep 
        }

        f = f_new; // Update solution
        results.push_back(f); // Save current state
    }

    // Save results to a file
    save_to_csv(filename, results, dx, dt);

    std::string return_message = "Simulations completed, results saved in: ";
    std::cout << return_message << filename << std::endl;

    return 0;
}