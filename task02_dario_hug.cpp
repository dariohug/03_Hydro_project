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
        std::cerr << "Unknown profile:" << profile;
        return;
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
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << "<Output file name>"  << "<starting condition>" << "<method>" << "<resolution>" << "<diffusion coefficient>" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    std::string starting_profile = argv[2];
    std::string method = argv[3];
    int resolution = std::atoi(argv[4]);
    double diffusion_coeficient = std::atof(argv[5]);

    // Parameters
    const double v0 = 1.0;
    const double x_min = 0.0, x_max = 1.0;
    int N = resolution;
    double dx = (x_max - x_min) / N;
    double dt = dx * dx / (4 * diffusion_coeficient); // CFL Condition
    double t_max = 1.0;
    int timesteps = static_cast<int>(t_max / dt);

    std::vector<double> f(N, 0.0); // Initialize the starting condition grid
    std::vector<double> f_new(N, 0.0); // Initialize the next timestep grid

    // use step starting condition or gaussian starting condition
    initialize(f, dx, starting_profile);

    // store data for output as a vector
    std::vector<std::vector<double>> results;
    results.push_back(f);

    if (method == "firstOrderEuler") {
        // Time integration loop
        for (int t = 0; t < timesteps; ++t) {
            for (int i = 0; i < N; ++i) {
                
                int i_upwind = (i == 0) ? N - 1 : i - 1;    // Periodic boundary for i-1
                int i_downwind = (i == N - 1) ? 0 : i + 1;  // Periodic boundary for i+1

                double d2fdx2 = (f[i_downwind] - 2 * f[i] + f[i_upwind]) / dx; // Second spatial derivative 

                f_new[i] = f[i] + dt * d2fdx2; // Firts-order-Euler 
            }
            f = f_new; // f_new has the updated values for the next timestep
            results.push_back(f); // Save current state
        }

    } else if (method == "RK2") {
        // Time integration loop
        for (int t = 0; t < timesteps; ++t) {
            std::vector<double> k1(N, 0.0);
            std::vector<double> k2(N, 0.0);

            // Compute k1 (Euler step)
            for (int i = 0; i < N; ++i) {
                int i_upwind = (i == 0) ? N - 1 : i - 1;    // Periodic boundary for i-1
                int i_downwind = (i == N - 1) ? 0 : i + 1;  // Periodic boundary for i+1

                double d2fdx2 = (f[i_downwind] - 2 * f[i] + f[i_upwind]) / dx; // Second spatial derivative
                k1[i] = dt * d2fdx2; // First step
            }

            // Compute k2 (use intermediate state)
            for (int i = 0; i < N; ++i) {
                int i_upwind = (i == 0) ? N - 1 : i - 1;    // Periodic boundary for i-1
                int i_downwind = (i == N - 1) ? 0 : i + 1;  // Periodic boundary for i+1

                double d2fdx2 = ((f[i] + k1[i_downwind]) - 2 * (f[i] + k1[i]) + (f[i_upwind] + k1[i_upwind])) / dx;
                k2[i] = dt * d2fdx2; // Second step
            }

            // Update solution with RK2 formula
            for (int i = 0; i < N; ++i) {
                f_new[i] = f[i] + 0.5 * (k1[i] + k2[i]);
            }
        f = f_new; // Update f for the next timestep
        results.push_back(f); // Save current state
        }
    } else {
        std::cerr << "Undefined Method. Methods include: [forwardEuler, RK2]" << std::endl;
        return 1;
            }

    // Save results to a file
    save_to_csv(filename, results, dx, dt);

    std::cout << "Simulations completed, results saved in: " << filename << std::endl;

    return 0;
}