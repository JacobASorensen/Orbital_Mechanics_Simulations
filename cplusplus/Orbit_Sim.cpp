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
const double DEFAULT_ERROR_MARGIN = 10e-10;

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

std::tuple<std::vector<double>, std::vector<std::vector<double>>> RKF_4_5(
    std::vector<double> y0,
    std::vector<double> masses,
    double t0,
    double t_end,
    double h,
    double G = DEFAULT_GRAVITY_CONSTANT,
    double acceptable_error = DEFAULT_ERROR_MARGIN,
    double upper_bound_fraction = 1.0,
    double lower_bound_fraction = 1.0

) {
    // Variables to store time and y values
    std::vector<double> t_values;
    std::vector<std::vector<double>> y_values;
    double original_step_size = h;
    double t = t0;
    std::vector<double> y = y0;
    bool record_previous_data = true;
    double old_h = h;
    const double A[7] = {
        0.0,0.0,2.0/9.0,1.0/3.0,3.0/4.0,1.0,5.0/6.0
    };

    const double B[7][6] = {
        {0.0,0.0,0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0,0.0,0.0},
        {0.0,2.0/9.0,0.0,0.0,0.0,0.0},
        {0.0,1.0/12.0,1.0/4.0,0.0,0.0,0.0},
        {0.0,69.0/128.0,-243.0/128.0,135.0/64.0,0.0,0.0},
        {0.0,-17.0/12.0,27.0/4.0,-27.0/5.0,16.0/15.0,0.0},
        {0.0,65.0/432.0,-5.0/16.0,13.0/16.0,4.0/27.0,5.0/144.0}
    };
    const double CH[7] = {
        0.0,
        47.0/450.0,
        0.0,
        12.0/25.0,
        32.0/225.0,
        1.0/30.0,
        6.0/25.0
    };
    const double CT[7] = {
        1.0/150.0,
        0.0,
        -3.0/100.0,
        16.0/75.0,
        1.0/20.0,
        -6.0/25.0
    };
    // for (size_t i = 0; i < 7; ++i) {
    //     cout << "A[" << i << "] is " << A[i] << " ";
    // }
    // cout  << endl;
    // for (size_t i = 0; i < 7; ++i) {
    //     for(size_t x = 0; x < 6; x++) {
    //         if(B[i][x] != 0.0) {
    //             cout << "B[" << i << "][" << x << "] is " << B[i][x] << " ";
    //         }
    //     }
    //     cout << endl;
    // }
    // cout << endl;
    // cout  << endl;
    // for (size_t i = 0; i < 7; ++i) {
    //     cout << "CH[" << i << "] is " << CH[i] << " ";
    // }
    // cout  << endl;
    // for (size_t i = 0; i < 7; ++i) {
    //     cout << "CT[" << i << "] is " << CT[i] << " ";
    // }
    // cout  << endl;
    // Runge-Kutta Loop
    while (t <= t_end) {
        if(record_previous_data) {
            t_values.push_back(t);
            y_values.push_back(y);
        }
        record_previous_data = true;
        std::vector<double> k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size()), k5(y.size()), k6(y.size());
        std::vector<double> temp(y.size());
        

        // k1
        std::vector<double> f_t_y = n_body_ODE_system(t, y,masses,G);
        for (size_t i = 0; i < y.size(); ++i) {
            k1[i] = h * f_t_y[i];
        }
        
        // k2
        for (size_t i = 0; i < y.size(); ++i) {
            temp[i] = y[i] + B[2][1] * k1[i];
        }
        f_t_y = n_body_ODE_system(t + A[2] * h, temp,masses,G);
        for (size_t i = 0; i < y.size(); ++i) {
            k2[i] = h * f_t_y[i];
        }

        // k3
        for (size_t i = 0; i < y.size(); ++i) {
            temp[i] = y[i] + B[3][1] * k1[i] + B[3][2] * k2[i];
        }
        f_t_y = n_body_ODE_system(t + A[3] * h, temp,masses,G);
        for (size_t i = 0; i < y.size(); ++i) {
            k3[i] = h * f_t_y[i];
        }

        // k4
        for (size_t i = 0; i < y.size(); ++i) {
            temp[i] = y[i] + B[4][1] * k1[i] + B[4][2] * k2[i] + B[4][3] * k3[i];
        }
        f_t_y = n_body_ODE_system(t + A[4] * h, temp,masses,G);
        for (size_t i = 0; i < y.size(); ++i) {
            k4[i] = h * f_t_y[i];
        }

        // k5
        for (size_t i = 0; i < y.size(); ++i) {
            temp[i] = y[i] + B[5][1] * k1[i] + B[5][2] * k2[i] + B[5][3] * k3[i] + B[5][4] * k4[i];
        }
        f_t_y = n_body_ODE_system(t + A[5] * h, temp,masses,G);
        for (size_t i = 0; i < y.size(); ++i) {
            k5[i] = h * f_t_y[i];
        }

        // k6
        for (size_t i = 0; i < y.size(); ++i) {
            temp[i] = y[i] + B[6][1] * k1[i] + B[6][2] * k2[i] + B[6][3] * k3[i] + B[6][4] * k4[i] + B[6][5] * k5[i];
        }
        f_t_y = n_body_ODE_system(t + A[6] * h, temp,masses,G);
        for (size_t i = 0; i < y.size(); ++i) {
            k6[i] = h * f_t_y[i];
        }
        // cout << "k1 :";
        // for (size_t i = 0; i < k1.size(); ++i) {
        //     cout << k1[i] << " ";
        // }
        // cout << endl;
        // cout << "k2 :";
        // for (size_t i = 0; i < k1.size(); ++i) {
        //     cout << k2[i] << " ";
        // }
        // cout << endl;
        // cout << "k3 :";
        // for (size_t i = 0; i < k1.size(); ++i) {
        //     cout << k3[i] << " ";
        // }
        // cout << endl;
        // cout << "k4 :";
        // for (size_t i = 0; i < k1.size(); ++i) {
        //     cout << k4[i] << " ";
        // }
        // cout << endl;
        // cout << "k5 :";
        // for (size_t i = 0; i < k1.size(); ++i) {
        //     cout << k5[i] << " ";
        // }
        // cout << endl;
        // cout << "k6 :";
        // for (size_t i = 0; i < k1.size(); ++i) {
        //     cout << k6[i] << " ";
        // }
        // cout << endl;

        
        
        // Update y and t
        for (size_t i = 0; i < y.size(); ++i) {
            y[i] += CH[1] * k1[i] + CH[2] * k2[i] + CH[3] * k3[i] + CH[4] * k4[i] + CH[5] * k5[i] + CH[6] * k6[i];
            // cout << "we are adding: " << CH[1] * k1[i] + CH[2] * k2[i] + CH[3] * k3[i] + CH[4] * k4[i] + CH[5] * k5[i] + CH[6] * k6[i];
            // cout << " for time " << t << " to y " << i << " which is: " << y[i] << endl;
        }

        std::vector<double> t_vec(y.size());
        double truncation_error = 0.0;
        double trunc_temp = 0.0;
        for (size_t i = 0; i < y.size(); ++i) {
            trunc_temp = CT[1] * k1[i] + CT[2] * k2[i] + CT[3] * k3[i] + CT[4] * k4[i] + CT[5] * k5[i] + CT[6] * k6[i];
            truncation_error += trunc_temp*trunc_temp;
        }
        truncation_error = pow(truncation_error,0.5);

        old_h = h;
        if(truncation_error > acceptable_error) {
            record_previous_data = false;
        }
        h = 0.9*h*pow(acceptable_error/truncation_error,1/5);

        if(h < lower_bound_fraction*original_step_size) {
            h = lower_bound_fraction*original_step_size;
            record_previous_data = true;
        } else if(h > upper_bound_fraction*original_step_size) {
            h = upper_bound_fraction*original_step_size;
            record_previous_data = true;
        }
        if(record_previous_data) {
            t += old_h;
        }
        
    }

    return std::make_tuple(t_values, y_values);
}

// pybind11 wrapper
PYBIND11_MODULE(Orbit_Sim, m) {
    // m.def("n_body_ODE_system", &n_body_ODE_system, "ODE system for N-body problem");
    m.def("runge_kutta_system", &runge_kutta_system, "Fourth-order Runge-Kutta for system of ODEs");
    m.def("runge_kutta_fehlberg_system", &RKF_4_5, "Fourth-order Runge-Kutta with fifth order error estimation and variable step-size for system of ODEs");
}


