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

void solve_advection_equation(std::vector<std::vector<double>> &f, int N, double dx, double dy, double dt, double vx0, double vy0, const std::string &slope_lim_method, bool solve_x)
{
    std::vector<std::vector<double>> f_new(N, std::vector<double>(N, 0.0));

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            int i_upwind = (i - 1 + N) % N;
            int i_downwind = (i + 1) % N;
            int j_upwind = (j - 1 + N) % N;
            int j_downwind = (j + 1) % N;

            if (solve_x)
            {
                // Solve MUSCL for x-direction
                double denominator_r_right_x = (f[i_downwind][j] - f[i][j] == 0.0) ? 1e-9 : f[i_downwind][j] - f[i][j];
                double r_right_x = (f[i][j] - f[i_upwind][j]) / denominator_r_right_x;

                double denominator_r_left_x = (f[i][j] - f[i_upwind][j] == 0.0) ? 1e-9 : f[i][j] - f[i_upwind][j];
                double r_left_x = (f[i_upwind][j] - f[i][j]) / denominator_r_left_x;

                double phi_right_x = (slope_lim_method == "vanleer") ? (r_right_x + std::abs(r_right_x)) / (1.0 + std::abs(r_right_x)) : std::max(0.0, std::min(r_right_x, 1.0));
                double phi_left_x = (slope_lim_method == "vanleer") ? (r_left_x + std::abs(r_left_x)) / (1.0 + std::abs(r_left_x)) : std::max(0.0, std::min(r_left_x, 1.0));

                double F_right_x = vx0 * ((vx0 > 0) ? (f[i][j] + 0.5 * phi_right_x * (f[i_downwind][j] - f[i][j])) : f[i_downwind][j]);
                double F_left_x = vx0 * ((vx0 > 0) ? (f[i_upwind][j] + 0.5 * phi_left_x * (f[i][j] - f[i_upwind][j])) : f[i][j]);

                f_new[i][j] = f[i][j] - (dt / dx) * (F_right_x - F_left_x);
            }
            else
            {
                // Solve MSUCL for y-direction
                double denominator_r_right_y = (f[i][j_downwind] - f[i][j] == 0.0) ? 1e-9 : f[i][j_downwind] - f[i][j];
                double r_right_y = (f[i][j] - f[i][j_upwind]) / denominator_r_right_y;

                double denominator_r_left_y = (f[i][j] - f[i][j_upwind] == 0.0) ? 1e-9 : f[i][j] - f[i][j_upwind];
                double r_left_y = (f[i][j_upwind] - f[i][j]) / denominator_r_left_y;

                double phi_right_y = (slope_lim_method == "vanleer") ? (r_right_y + std::abs(r_right_y)) / (1.0 + std::abs(r_right_y)) : std::max(0.0, std::min(r_right_y, 1.0));
                double phi_left_y = (slope_lim_method == "vanleer") ? (r_left_y + std::abs(r_left_y)) / (1.0 + std::abs(r_left_y)) : std::max(0.0, std::min(r_left_y, 1.0));

                double F_right_y = vy0 * ((vy0 > 0) ? (f[i][j] + 0.5 * phi_right_y * (f[i][j_downwind] - f[i][j])) : f[i][j_downwind]);
                double F_left_y = vy0 * ((vy0 > 0) ? (f[i][j_upwind] + 0.5 * phi_left_y * (f[i][j] - f[i][j_upwind])) : f[i][j]);

                f_new[i][j] = f[i][j] - (dt / dy) * (F_right_y - F_left_y);
            }
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

    const double vx0 = 1.0, vy0 = 1.0;
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
        solve_advection_equation(f, N, dx, dy, dt / 2, vx0, vy0, slope_lim_method, true); // Half a step in x direction
        solve_advection_equation(f, N, dx, dy, dt, vx0, vy0, slope_lim_method, false);    // Full step in y Direction
        solve_advection_equation(f, N, dx, dy, dt / 2, vx0, vy0, slope_lim_method, true); // Half a step in x direction
        results.push_back(f);
    }

    save_to_csv(filename, results, dx, dt);
    std::cout << "Simulations completed, results saved in: " << filename << std::endl;

    return 0;
}