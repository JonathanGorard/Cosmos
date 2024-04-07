#include "RelativisticEulerEquationOfState.hpp"

double RelativisticEulerEquationOfState::computeSpecificEnthalpy(double density, double pressure, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    
    return 1.0 + ((pressure / density) * (adiabaticIndex / (adiabaticIndex - 1.0)));
}

double RelativisticEulerEquationOfState::computeSoundSpeed(double density, double pressure, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    
    double numerator = (adiabaticIndex * pressure) / density;
    double denominator = 1.0 + ((pressure / density) * (adiabaticIndex / (adiabaticIndex - 1.0)));
    
    return sqrt(numerator / denominator);
}
