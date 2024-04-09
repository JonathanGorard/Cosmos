#include <iostream>
#include <cmath>

const double g = 9.81; // Acceleration due to gravity (m/s^2)

// Function to calculate the time of flight
double timeOfFlight(double v0, double angle) {
    double angleRad = (angle * M_PI) / 180; // Convert angle to radians
    return (2 * v0 * sin(angleRad)) / g;
}

// Function to calculate maximum height
double maxHeight(double v0, double angle) {
    double angleRad = (angle * M_PI) / 180; // Convert angle to radians
    return (pow(v0, 2) * pow(sin(angleRad), 2)) / (2 * g);
}

// Function to calculate range
double range(double v0, double angle) {
    double angleRad = (angle * M_PI) / 180; // Convert angle to radians
    return (pow(v0, 2) * sin(2 * angleRad)) / g;
}

int main() {
    double v0, angle;
    std::cout << "Enter initial velocity (m/s): ";
    std::cin >> v0;
    std::cout << "Enter launch angle (degrees): ";
    std::cin >> angle;

    std::cout << "Time of Flight: " << timeOfFlight(v0, angle) << " seconds" << std::endl;
    std::cout << "Maximum Height: " << maxHeight(v0, angle) << " meters" << std::endl;
    std::cout << "Range: " << range(v0, angle) << " meters" << std::endl;

    return 0;
}
