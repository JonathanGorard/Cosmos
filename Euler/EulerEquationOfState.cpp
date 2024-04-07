#include "EulerEquationOfState.hpp"

double EulerEquationOfState::computeSpecificInternalEnergy(double density, double pressure, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    
    return pressure / ((adiabaticIndex - 1.0) * density);
}

double EulerEquationOfState::computePressure(double density, double specificInternalEnergy, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    
    return specificInternalEnergy * (adiabaticIndex - 1.0) * density;
}

double EulerEquationOfState::computeSoundSpeed(double density, double pressure, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    
    return sqrt((adiabaticIndex * pressure) / density);
}
