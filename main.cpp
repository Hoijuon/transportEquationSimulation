#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <utility>
#include <limits>

// path is file path.
// const std::string path = "./txt/1/1-60/";
// const std::string path = "./txt/1/2-40/";
// const std::string path = "./txt/1/3-30/";
// const std::string path = "./txt/3/2x3+3x5+1x22/";
// const std::string path = "./txt/3/2x4+3x5+1x21/";
// const std::string path = "./txt/3/2x5+3x4+1x21/";
// const std::string path = "./txt/3/2x6+3x4+1x20/";
// const std::string path = "./txt/3/2x7+3x3+1x20/";
const std::string path = "./txt/3/2x8+3x3+1x19/";
// const std::string path = "./txt/3/2x9+3x2+1x19/";
// const std::string path = "./txt/3/2x10+3x2+1x18/";
// const std::string path = "./txt/air/";

const double L = 130.0;  // Model size, in units of μm.
const double dy = 1.30;  // Spatial step size, in units of μm.
const int N = static_cast<int>(L / dy);  // Number of spatial discrete points.
const double dt = 0.01;  // Time step size, in units of s.
const double alphaC = 1.5;  // Constant alphaC.
const double binderRadius = 0.1e-6;  // Binder radius, in units of m.
const double Porosity = 0.52;  // Porosity.
const double kB = 1.380649e-23;  // Boltzmann constant, in units of J/K.
const double thetaL_ini = 0.7368989339236663;
const double thetaL_end = 0.01;
const double thetaL_star = 0.02;  // Residual Saturation.
const double kc = 0.8;  // Processing coefficient for velocity.

// Function: Calculate dynamic viscosity.
double dynamicViscosity(double T) {
    return 0.0318e-3 * exp(484.3726 / (T - 273.15 + 120.2202));  // Unit: Pa·s.
}

// Function: Calculate the step function step(thetaL).
double step(double thetaL) {
    // return std::pow(1 + std::exp(-((thetaL -thetaL_star / 2) * 1e4)), -1);  // Sigmoid function.
    if (thetaL > thetaL_star) { return 1; }
    else if (thetaL < thetaL_star && thetaL > 0.005) { return std::pow(thetaL / thetaL_star, 1); }
    else { return 0; }
}

// Function: Calculate the effective diffusion coefficient Deff. The formula is Deff(t) = ϵ(t)/τC(t) * D(t), where τC(t) = ϵ(t)^(1 - αC).
double calculateDeff(double temperature, double thetaL) {
    double D = kB * temperature / (6 * M_PI * dynamicViscosity(temperature) * binderRadius);  //  Stokes - Einstein function.
    double tauC = std::pow(Porosity, 1 - alphaC);
    double Deff = Porosity / tauC * D * step(thetaL); 
    return Deff;
}

// Function: Advection - Diffusion Equation (Law of Mass Conservation), using the Finite Difference Method.
void updateConcentration(int m, std::vector<double>& c, double& liquidLevel, double& prevLiquidLevel, double v, double temperature, double thetaL) {
    int liquidIndex = static_cast<int>(std::round(liquidLevel / dy));
    int prevLiquidIndex = static_cast<int>(std::round(prevLiquidLevel / dy));

    // Convert concentration to amount of substance and calculate the total amount of substances.
    double total_before = 0.0;
    std::vector<double> cm(c.size());
    for (int i = 0; i < c.size(); ++i) {
        cm[i] = c[i] * dy * 1e-6;  // Consider one - dimensional case, where both x and z are 1m.
        total_before += cm[i];
    }

    // Handle the substances between the liquid level at the previous time step and that at the current time step.
    double sumAbove = 0.0;
    for (int i = liquidIndex; i < prevLiquidIndex; ++i) {
        sumAbove += cm[i];
    }
    cm[liquidIndex - 1] += sumAbove;

    std::vector<double> cmNew(c.size());
    double Deff = calculateDeff(temperature, thetaL);
    double Reff = Deff * dt / (dy * dy * 1e-12);
    if (liquidIndex >= 2) {
        for (int i = 1; i < liquidIndex - 1; ++i) {
            // The central difference method: Fick's law states that the diffusion flux is proportional to the concentration gradient, and its expression is J = -D * dc/dx.
            cmNew[i] = cm[i] + Reff * (cm[i + 1] - 2 * cm[i] + cm[i - 1]);
        }
        // Bottom boundary, reflective boundary condition. Assume no substances flow in or out.
        cmNew[0] = cm[0] + Reff * (cm[1] - cm[0]);
        // Top boundary, reflective boundary condition. Assume no substances flow in or out.
        cmNew[liquidIndex - 1] = cm[liquidIndex - 1] + Reff * (cm[liquidIndex - 1] - cm[liquidIndex - 2]);
    } else {
        cmNew[0] = cm[0];
    }

    // Source term: Adjust the last element to ensure the conservation of the total amount of solute.
    double total_after = 0.0;
    for (int i = 0; i < liquidIndex; ++i) {
        total_after += cmNew[i];
    }
    cmNew[liquidIndex - 1] += total_before - total_after;

    // Convert back to concentration.
    for (int i = 0; i < (c.size()); ++i) {
        c[i] = cmNew[i] / dy * 1e6;
    }

    // Update the liquid level position.
    prevLiquidLevel = liquidLevel;
    liquidLevel -= v * dt * 1e6; 
    if (liquidLevel < 0) {
        liquidLevel = 0;
    }
}

// Function: Read data from file.
std::vector<std::pair<double, double>> readFileData(const std::string& filename) {
    std::vector<std::pair<double, double>> data;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double time, value;
        if (iss >> time >> value) {
            data.emplace_back(time, value);
        }
    }
    return data;
}

// Function: Transpose Matrix.
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    std::vector<std::vector<double>> result(cols, std::vector<double>(rows));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}

int main() {
    std::vector<double> u(N, 1);  // Initial Concentration.
    std::vector<std::pair<double, double>> temperatureData = readFileData(path + "temperature.txt");  // Read temperature from temperature.txt, in units of K.
    std::vector<std::pair<double, double>> velocityData_ini = readFileData(path + "v.txt");  // Read velocity from v_cf.txt, in units of m/s.
    std::vector<std::pair<double, double>> thetaLData = readFileData(path + "thetaL.txt");  // Read MC from thetaL.txt, in units of 1.

    // Processing velocity.
    size_t row = velocityData_ini.size();
    std::vector<std::pair<double, double>> velocityData_new(row);
    for (size_t i = 0; i < row; ++i) {
        velocityData_new[i].first = velocityData_ini[i].first;
        velocityData_new[i].second = 0;
    }
    for (size_t i = 0; i < row; ++i) {
        velocityData_new[static_cast<int>(velocityData_ini[i].first * kc / dt) + 1].second = velocityData_ini[i].second / kc;
    }

    // main.
    double t_out = 0;
    for (size_t i = 0; i < thetaLData.size(); ++i) {
        if (thetaLData[i].second <= thetaL_end) {
            t_out = thetaLData[i].first;
            break;
        }
    }
    std::cout << "t_out = " << t_out << " s." << std::endl;
   
    std::vector<std::vector<double>> concentrationMatrix;  // Store the concentration distribution.
    concentrationMatrix.push_back(u);

    int min_size = std::min(std::min(temperatureData.size(), velocityData_new.size()), thetaLData.size());
    double y = L;  double prevY = y;  // Store the liquid level position.
    for (int m = 1; m < min_size; ++m) {
        double v = velocityData_new[m].second;
        double temperature = temperatureData[m].second;
        double thetaL = thetaLData[m].second;

        updateConcentration(m, u, y, prevY, v, temperature, thetaL);
        concentrationMatrix.push_back(u);

        // Output the concentration distribution at the current moment.
        if (m * dt == t_out) {
            double sum = 0;
            std::cout << "Time step " << m << ": ";
            for (int j = 0; j < N; ++j) {
                std::cout << u[j] << "  ";
                sum += u[j] * dy;
            }
            std::cout << sum << "  ";
            std::cout << std::endl;
        }
    }

    std::vector<std::vector<double>> transposedMatrix = transpose(concentrationMatrix);

    std::string input;
    std::cout << "Save concentration? (y/n)" << std::endl;
    std::getline(std::cin, input);
    if (input.empty()) {
    } else if (input =="y") {  // Output the results to a file.
        std::ofstream outputFile(path + "concentration_data.txt");  // Open file.
        if (!outputFile.is_open()) {
            std::cerr << "Failed to open the output file." << std::endl;
            return 1;
        }
        for (int i = 0; i < transposedMatrix.size(); ++i) {
            outputFile << dy / 2.0 + (i * dy) << "\t";
            for (int j = 0; j < transposedMatrix[i].size(); ++j) {
                if (j * dt == t_out) {
                    outputFile << transposedMatrix[i][j] << "\t";
                }
            }
            outputFile << std::endl;
        }
        outputFile.close();  // Close file.
    } else {
    }

    return 0;
}
