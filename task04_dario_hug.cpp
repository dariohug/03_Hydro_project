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

void solve_advection_equation(std::vector<std::vector<double>> &f, int N, double dx, double dt, double v0, const std::string &slope_lim_method)
{
    std::vector<std::vector<double>> f_new(N, std::vector<double>(N, 0.0));
    const double epsilon = 1e-6;

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            int i_upwind = (i - 1 + N) % N;
            int i_downwind = (i + 1) % N;
            int j_upwind = (j - 1 + N) % N;
            int j_downwind = (j + 1) % N;

            double r_right = (f[i][j] - f[i_upwind][j]) / (f[i_downwind][j] - f[i][j] + epsilon);
            double r_left = (f[i_upwind][j] - f[i][j]) / (f[i][j] - f[i_upwind][j] + epsilon);

            double phi = (slope_lim_method == "vanleer")  ? (r_right + std::abs(r_right)) / (1.0 + std::abs(r_right))
                         : (slope_lim_method == "minmod") ? std::max(0.0, std::min(1.0, r_right))
                                                          : 1.0;

            double F_right = v0 * ((v0 > 0) ? (f[i][j] + 0.5 * phi * (f[i_downwind][j] - f[i][j])) : f[i_downwind][j]);
            double F_left = v0 * ((v0 > 0) ? (f[i_upwind][j] + 0.5 * phi * (f[i][j] - f[i_upwind][j])) : f[i][j]);

            f_new[i][j] = f[i][j] - (dt / dx) * (F_right - F_left);
        }
    }
    f = f_new;
}

int main(int argc, char *argv[])
{
    if (argc < 6)
    {
        std::cerr << "Usage: " << argv[0] << " <timesteps> <filename> <Resolution> <starting profile> <slope_limiter_method>" << std::endl;
        return 1;
    }

    int timesteps = std::atoi(argv[1]);
    std::string filename = argv[2];
    int resolution = std::atoi(argv[3]);
    std::string starting_profile = argv[4];
    std::string slope_lim_method = argv[5];

    const double vx0 = 1.0 / std::sqrt(2), vy0 = 1.0 / std::sqrt(2);
    const double x_min = 0.0, x_max = 1.0, y_min = 0.0, y_max = 1.0;
    int N = resolution;
    double dx = (x_max - x_min) / N, dy = (y_max - y_min) / N;
    double dt = 0.5 * (dx * dy) / std::sqrt(((dx * vx0) * (dx * vx0)) + ((dy * vy0) * (dy * vy0)));

    std::vector<std::vector<double>> f(N, std::vector<double>(N, 0.0));
    initialize(f, dx, starting_profile);

    std::vector<std::vector<std::vector<double>>> results;
    results.push_back(f);

    for (int t = 0; t < timesteps; ++t)
    {
        solve_advection_equation(f, N, dy, dt / 2, vy0, slope_lim_method);
        solve_advection_equation(f, N, dx, dt, vx0, slope_lim_method);
        solve_advection_equation(f, N, dy, dt / 2, vy0, slope_lim_method);
        results.push_back(f);
    }

    save_to_csv(filename, results, dx, dt);
    std::cout << "Simulations completed, results saved in: " << filename << std::endl;

    return 0;
}