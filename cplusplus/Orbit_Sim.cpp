// We want to convert the simulation code, i.e. the n_body_ODE_system function and the runge_kutta_system functions in here

#include <iostream>
#include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
#include <pybind11/stl.h>
// #include <pybind11/complex.h>
#include <vector>
#include <functional>

#include <map>
#include <cmath>

namespace py = pybind11;
using namespace std;

// Constants
const double DEFAULT_GRAVITY_CONSTANT = 6.67430e-11; // Gravitational constant
//
// G = 6.67259*10**-20

// Masses of the bodies (this will be global or passed to the function)


// Function to calculate the ODE system for the N-body problem
std::vector<double> n_body_ODE_system(double time, const std::vector<double>& y,const std::vector<double>& masses, const double G) {
    int num = masses.size(); // Number of bodies
    std::vector<std::pair<std::vector<double>, std::vector<double>>> bodies;
    
    // Extract positions and velocities from the input vector y
    for (int nth = 0; nth < num; ++nth) {
        std::vector<double> pos(3), vel(3);
        for (int v = 0; v < 3; ++v) {
            vel[v] = y[nth * 3 + v];
        }
        for (int x = 0; x < 3; ++x) {
            pos[x] = y[num * 3 + nth * 3 + x];
        }
        bodies.push_back({pos, vel});
    }

    // Map to store distances between bodies
    std::map<std::pair<int, int>, double> distances;
    
    // Vectors to store accelerations and new positions
    std::vector<double> accelerations, new_positions;

    // Calculate accelerations for each body
    for (int first = 0; first < num; ++first) {
        std::vector<double> acceleration(3, 0.0); // Initialize to zero

        for (int second = 0; second < num; ++second) {
            if (first != second) {
                // Check if the distance has already been computed
                auto dist_key = std::make_pair(first, second);
                if (distances.find(dist_key) == distances.end()) {
                    double distance_squared = 0.0;

                    // Calculate distance between first and second body
                    for (int i = 0; i < 3; ++i) {
                        distance_squared += std::pow(bodies[second].first[i] - bodies[first].first[i], 2);
                    }
                    double distance = std::sqrt(distance_squared);
                    distances[dist_key] = distance;
                    distances[std::make_pair(second, first)] = distance; // Symmetric distance
                }

                // Calculate acceleration for the first body due to the second body
                double distance = distances[dist_key];
                for (int i = 0; i < 3; ++i) {
                    acceleration[i] += G * masses[second] * (bodies[second].first[i] - bodies[first].first[i]) / std::pow(distance, 3);
                }
            }
        }
        accelerations.insert(accelerations.end(), acceleration.begin(), acceleration.end());
    }

    // Append new velocities (velocities are equivalent to the time derivative of position)
    for (int nth = 0; nth < num; ++nth) {
        new_positions.insert(new_positions.end(), bodies[nth].second.begin(), bodies[nth].second.end());
    }

    // Return accelerations followed by new positions (as derivatives)
    accelerations.insert(accelerations.end(), new_positions.begin(), new_positions.end());
    return accelerations;
}
namespace py = pybind11;

// Typedef for a function that takes double t and a vector y, and returns a vector
// using SystemFunction = function<vector<double>(double, const vector<double>&)>;

std::tuple<std::vector<double>, std::vector<std::vector<double>>> runge_kutta_system(
    std::vector<double> y0,
    std::vector<double> masses,
    double t0,
    double t_end,
    double h,
    double G = DEFAULT_GRAVITY_CONSTANT
) {
    // Variables to store time and y values
    std::vector<double> t_values;
    std::vector<std::vector<double>> y_values;

    double t = t0;
    std::vector<double> y = y0;

    // Runge-Kutta Loop
    while (t <= t_end) {
        t_values.push_back(t);
        y_values.push_back(y);
        std::vector<double> k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size());
        std::vector<double> temp(y.size());

        // k1
        std::vector<double> f_t_y = n_body_ODE_system(t, y,masses,G);
        for (size_t i = 0; i < y.size(); ++i) {
            k1[i] = h * f_t_y[i];
        }

        // k2
        for (size_t i = 0; i < y.size(); ++i) {
            temp[i] = y[i] + 0.5 * k1[i];
        }
        f_t_y = n_body_ODE_system(t + 0.5 * h, temp,masses,G);
        for (size_t i = 0; i < y.size(); ++i) {
            k2[i] = h * f_t_y[i];
        }

        // k3
        for (size_t i = 0; i < y.size(); ++i) {
            temp[i] = y[i] + 0.5 * k2[i];
        }
        f_t_y = n_body_ODE_system(t + 0.5 * h, temp,masses,G);
        for (size_t i = 0; i < y.size(); ++i) {
            k3[i] = h * f_t_y[i];
        }

        // k4
        for (size_t i = 0; i < y.size(); ++i) {
            temp[i] = y[i] + k3[i];
        }
        f_t_y = n_body_ODE_system(t + h, temp,masses,G);
        for (size_t i = 0; i < y.size(); ++i) {
            k4[i] = h * f_t_y[i];
        }

        // Update y and t
        for (size_t i = 0; i < y.size(); ++i) {
            y[i] += (1.0 / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        }

        t += h;
    }

    return std::make_tuple(t_values, y_values);
}

// pybind11 wrapper
PYBIND11_MODULE(Orbit_Sim, m) {
    // m.def("n_body_ODE_system", &n_body_ODE_system, "ODE system for N-body problem");
    m.def("runge_kutta_system", &runge_kutta_system, "Fourth-order Runge-Kutta for system of ODEs");
}

