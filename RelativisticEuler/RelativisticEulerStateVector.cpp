#include "RelativisticEulerStateVector.hpp"

RelativisticEulerStateVector::RelativisticEulerStateVector()
{
    density = 1.0;
    xVelocity = 0.0;
    yVelocity = 0.0;
    zVelocity = 0.0;
    pressure = 1.0;
}

RelativisticEulerStateVector::RelativisticEulerStateVector(double newDensity, double newXVelocity, double newYVelocity, double newZVelocity, double newPressure)
{
    density = newDensity;
    xVelocity = newXVelocity;
    yVelocity = newYVelocity;
    zVelocity = newZVelocity;
    pressure = newPressure;
}

void RelativisticEulerStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    density = newPrimitiveVariableVector[0];
    xVelocity = newPrimitiveVariableVector[1];
    yVelocity = newPrimitiveVariableVector[2];
    zVelocity = newPrimitiveVariableVector[3];
    pressure = newPrimitiveVariableVector[4];
}

void RelativisticEulerStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    
    double relativisticDensity = newConservedVariableVector[0];
    double xMomentum = newConservedVariableVector[1];
    double yMomentum = newConservedVariableVector[2];
    double zMomentum = newConservedVariableVector[3];
    double relativisticEnergy = newConservedVariableVector[4];
    
    double densityConstant = relativisticDensity / sqrt(((relativisticEnergy + relativisticDensity) * (relativisticEnergy + relativisticDensity)) -
                                                        ((xMomentum * xMomentum) + (yMomentum * yMomentum) + (zMomentum * zMomentum)));
    double energyConstant = (relativisticDensity + relativisticEnergy) / sqrt(((relativisticEnergy + relativisticDensity) * (relativisticEnergy + relativisticDensity)) -
                                                                              ((xMomentum * xMomentum) + (yMomentum * yMomentum) + (zMomentum * zMomentum)));
    
    if (((relativisticEnergy + relativisticDensity) * (relativisticEnergy + relativisticDensity)) -
        ((xMomentum * xMomentum) + (yMomentum * yMomentum) + (zMomentum * zMomentum)) < pow(10.0, -8.0))
    {
        densityConstant = relativisticDensity / sqrt(pow(10.0, -8.0));
        energyConstant = (relativisticDensity + relativisticEnergy) / sqrt(pow(10.0, -8.0));
    }
    
    double constantTerm = -1.0 / (adiabaticIndex * adiabaticIndex);
    double linearTerm = -2.0 * densityConstant * ((adiabaticIndex - 1.0) / (adiabaticIndex * adiabaticIndex));
    double quadraticTerm = ((adiabaticIndex - 2.0) / adiabaticIndex) * ((energyConstant * energyConstant) - 1.0) + 1.0 -
    (densityConstant * densityConstant) * ((adiabaticIndex - 1.0) / adiabaticIndex) * ((adiabaticIndex - 1.0) / adiabaticIndex);
    double quarticTerm = (energyConstant * energyConstant) - 1.0;
    double etaConstant = 2.0 * densityConstant * ((adiabaticIndex - 1.0) / adiabaticIndex);
    
    double currentGuess = 1.0;
    int iterationCount = 0;
    
    while (iterationCount < 1000)
    {
        double function = (quarticTerm * (currentGuess * currentGuess * currentGuess) * (currentGuess - etaConstant)) +
        (quadraticTerm * (currentGuess * currentGuess)) + (linearTerm * currentGuess) + constantTerm;
        double functionDerivative = linearTerm + (2.0 * quadraticTerm * currentGuess) + (4.0 * quarticTerm * (currentGuess * currentGuess * currentGuess)) -
        (3.0 * etaConstant * quarticTerm * (currentGuess * currentGuess));
        
        double newGuess = currentGuess - (function / functionDerivative);
        
        if (abs(currentGuess - newGuess) < pow(10.0, -8.0))
        {
            iterationCount = 1000;
        }
        else
        {
            iterationCount += 1;
            currentGuess = newGuess;
        }
    }
    
    double lorentzFactor = 0.5 * energyConstant * currentGuess * (1.0 + sqrt(1.0 + (4.0 * ((adiabaticIndex - 1.0) / adiabaticIndex) *
                                                                                    ((1.0 - (densityConstant * currentGuess)) /
                                                                                     ((energyConstant * energyConstant) * (currentGuess * currentGuess))))));
    double specificEnthalpy = 1.0 / (densityConstant * currentGuess);
    
    density = relativisticDensity / lorentzFactor;
    xVelocity = xMomentum / (density * specificEnthalpy * (lorentzFactor * lorentzFactor));
    yVelocity = yMomentum / (density * specificEnthalpy * (lorentzFactor * lorentzFactor));
    zVelocity = zMomentum / (density * specificEnthalpy * (lorentzFactor * lorentzFactor));
    pressure = (density * specificEnthalpy * (lorentzFactor * lorentzFactor)) - relativisticDensity - relativisticEnergy;
    
    if (density < pow(10.0, -8.0))
    {
        density = pow(10.0, -8.0);
    }
    if (pressure < pow(10.0, -8.0))
    {
        pressure = pow(10.0, -8.0);
    }
}

vector<double> RelativisticEulerStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(5);
    
    primitiveVariableVector[0] = density;
    primitiveVariableVector[1] = xVelocity;
    primitiveVariableVector[2] = yVelocity;
    primitiveVariableVector[3] = zVelocity;
    primitiveVariableVector[4] = pressure;
    
    return primitiveVariableVector;
}

vector<double> RelativisticEulerStateVector::computeConservedVariableVector(EulerMaterialParameters materialParameters)
{
    double lorentzFactor = 1.0 / sqrt(1.0 - ((xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity)));
    
    if ((xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity) > 1.0 - pow(10.0, -8.0))
    {
        lorentzFactor = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
    }
    
    double specificEnthalpy = RelativisticEulerEquationOfState::computeSpecificEnthalpy(density, pressure, materialParameters);
    
    vector<double> conservedVariableVector(5);
    
    conservedVariableVector[0] = density * lorentzFactor;
    conservedVariableVector[1] = density * specificEnthalpy * (lorentzFactor * lorentzFactor) * xVelocity;
    conservedVariableVector[2] = density * specificEnthalpy * (lorentzFactor * lorentzFactor) * yVelocity;
    conservedVariableVector[3] = density * specificEnthalpy * (lorentzFactor * lorentzFactor) * zVelocity;
    conservedVariableVector[4] = (density * specificEnthalpy * (lorentzFactor * lorentzFactor)) - pressure - (density * lorentzFactor);
    
    return conservedVariableVector;
}

vector<double> RelativisticEulerStateVector::computeXFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    
    double relativisticDensity = conservedVariableVector[0];
    double xMomentum = conservedVariableVector[1];
    double yMomentum = conservedVariableVector[2];
    double zMomentum = conservedVariableVector[3];
    double relativisticEnergy = conservedVariableVector[4];
    
    double densityConstant = relativisticDensity / sqrt(((relativisticEnergy + relativisticDensity) * (relativisticEnergy + relativisticDensity)) -
                                                        ((xMomentum * xMomentum) + (yMomentum * yMomentum) + (zMomentum * zMomentum)));
    double energyConstant = (relativisticDensity + relativisticEnergy) / sqrt(((relativisticEnergy + relativisticDensity) * (relativisticEnergy + relativisticDensity)) -
                                                                              ((xMomentum * xMomentum) + (yMomentum * yMomentum) + (zMomentum * zMomentum)));
    
    if (((relativisticEnergy + relativisticDensity) * (relativisticEnergy + relativisticDensity)) -
        ((xMomentum * xMomentum) + (yMomentum * yMomentum) + (zMomentum * zMomentum)) < pow(10.0, -8.0))
    {
        densityConstant = relativisticDensity / sqrt(pow(10.0, -8.0));
        energyConstant = (relativisticDensity + relativisticEnergy) / sqrt(pow(10.0, -8.0));
    }
    
    double constantTerm = -1.0 / (adiabaticIndex * adiabaticIndex);
    double linearTerm = -2.0 * densityConstant * ((adiabaticIndex - 1.0) / (adiabaticIndex * adiabaticIndex));
    double quadraticTerm = ((adiabaticIndex - 2.0) / adiabaticIndex) * ((energyConstant * energyConstant) - 1.0) + 1.0 -
    (densityConstant * densityConstant) * ((adiabaticIndex - 1.0) / adiabaticIndex) * ((adiabaticIndex - 1.0) / adiabaticIndex);
    double quarticTerm = (energyConstant * energyConstant) - 1.0;
    double etaConstant = 2.0 * densityConstant * ((adiabaticIndex - 1.0) / adiabaticIndex);
    
    double currentGuess = 1.0;
    int iterationCount = 0;
    
    while (iterationCount < 1000)
    {
        double function = (quarticTerm * (currentGuess * currentGuess * currentGuess) * (currentGuess - etaConstant)) +
        (quadraticTerm * (currentGuess * currentGuess)) + (linearTerm * currentGuess) + constantTerm;
        double functionDerivative = linearTerm + (2.0 * quadraticTerm * currentGuess) + (4.0 * quarticTerm * (currentGuess * currentGuess * currentGuess)) -
        (3.0 * etaConstant * quarticTerm * (currentGuess * currentGuess));
        
        double newGuess = currentGuess - (function / functionDerivative);
        
        if (abs(currentGuess - newGuess) < pow(10.0, -8.0))
        {
            iterationCount = 1000;
        }
        else
        {
            iterationCount += 1;
            currentGuess = newGuess;
        }
    }
    
    double lorentzFactor = 0.5 * energyConstant * currentGuess * (1.0 + sqrt(1.0 + (4.0 * ((adiabaticIndex - 1.0) / adiabaticIndex) *
                                                                                    ((1.0 - (densityConstant * currentGuess)) /
                                                                                     ((energyConstant * energyConstant) * (currentGuess * currentGuess))))));
    double specificEnthalpy = 1.0 / (densityConstant * currentGuess);
    
    double computedDensity = relativisticDensity / lorentzFactor;
    double computedXVelocity = xMomentum / (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor));
    double computedYVelocity = yMomentum / (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor));
    double computedZVelocity = zMomentum / (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor));
    double computedPressure = (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor)) - relativisticDensity - relativisticEnergy;
    
    if (computedDensity < pow(10.0, -8.0))
    {
        computedDensity = pow(10.0, -8.0);
    }
    if (computedPressure < pow(10.0, -8.0))
    {
        computedPressure = pow(10.0, -8.0);
    }
    
    vector<double> fluxVector(5);
    
    fluxVector[0] = computedDensity * lorentzFactor * computedXVelocity;
    fluxVector[1] = computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) * (computedXVelocity * computedXVelocity) + computedPressure;
    fluxVector[2] = computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) * (computedYVelocity * computedXVelocity);
    fluxVector[3] = computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) * (computedZVelocity * computedXVelocity);
    fluxVector[4] = ((computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor)) - computedPressure - (computedDensity * lorentzFactor)) *
    computedXVelocity + (computedPressure * computedXVelocity);
    
    return fluxVector;
}

vector<double> RelativisticEulerStateVector::computeXFluxVector(EulerMaterialParameters materialParameters)
{
    return computeXFluxVector(computeConservedVariableVector(materialParameters), materialParameters);
}

vector<double> RelativisticEulerStateVector::computeYFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    
    double relativisticDensity = conservedVariableVector[0];
    double xMomentum = conservedVariableVector[1];
    double yMomentum = conservedVariableVector[2];
    double zMomentum = conservedVariableVector[3];
    double relativisticEnergy = conservedVariableVector[4];
    
    double densityConstant = relativisticDensity / sqrt(((relativisticEnergy + relativisticDensity) * (relativisticEnergy + relativisticDensity)) -
                                                        ((xMomentum * xMomentum) + (yMomentum * yMomentum) + (zMomentum * zMomentum)));
    double energyConstant = (relativisticDensity + relativisticEnergy) / sqrt(((relativisticEnergy + relativisticDensity) * (relativisticEnergy + relativisticDensity)) -
                                                                              ((xMomentum * xMomentum) + (yMomentum * yMomentum) + (zMomentum * zMomentum)));
    
    if (((relativisticEnergy + relativisticDensity) * (relativisticEnergy + relativisticDensity)) -
        ((xMomentum * xMomentum) + (yMomentum * yMomentum) + (zMomentum * zMomentum)) < pow(10.0, -8.0))
    {
        densityConstant = relativisticDensity / sqrt(pow(10.0, -8.0));
        energyConstant = (relativisticDensity + relativisticEnergy) / sqrt(pow(10.0, -8.0));
    }
    
    double constantTerm = -1.0 / (adiabaticIndex * adiabaticIndex);
    double linearTerm = -2.0 * densityConstant * ((adiabaticIndex - 1.0) / (adiabaticIndex * adiabaticIndex));
    double quadraticTerm = ((adiabaticIndex - 2.0) / adiabaticIndex) * ((energyConstant * energyConstant) - 1.0) + 1.0 -
    (densityConstant * densityConstant) * ((adiabaticIndex - 1.0) / adiabaticIndex) * ((adiabaticIndex - 1.0) / adiabaticIndex);
    double quarticTerm = (energyConstant * energyConstant) - 1.0;
    double etaConstant = 2.0 * densityConstant * ((adiabaticIndex - 1.0) / adiabaticIndex);
    
    double currentGuess = 1.0;
    int iterationCount = 0;
    
    while (iterationCount < 1000)
    {
        double function = (quarticTerm * (currentGuess * currentGuess * currentGuess) * (currentGuess - etaConstant)) +
        (quadraticTerm * (currentGuess * currentGuess)) + (linearTerm * currentGuess) + constantTerm;
        double functionDerivative = linearTerm + (2.0 * quadraticTerm  * currentGuess) + (4.0 * quarticTerm * (currentGuess * currentGuess * currentGuess)) -
        (3.0 * etaConstant * quarticTerm * (currentGuess * currentGuess));
        
        double newGuess = currentGuess - (function / functionDerivative);
        
        if (abs(currentGuess - newGuess) < pow(10.0, -8.0))
        {
            iterationCount = 1000;
        }
        else
        {
            iterationCount += 1;
            currentGuess = newGuess;
        }
    }
    
    double lorentzFactor = 0.5 * energyConstant * currentGuess * (1.0 + sqrt(1.0 + (4.0 * ((adiabaticIndex - 1.0) / adiabaticIndex) *
                                                                                    ((1.0 - (densityConstant * currentGuess)) /
                                                                                     ((energyConstant * energyConstant) * (currentGuess * currentGuess))))));
    double specificEnthalpy = 1.0 / (densityConstant * currentGuess);
    
    double computedDensity = relativisticDensity / lorentzFactor;
    double computedXVelocity = xMomentum / (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor));
    double computedYVelocity = yMomentum / (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor));
    double computedZVelocity = zMomentum / (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor));
    double computedPressure = (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor)) - relativisticDensity - relativisticEnergy;
    
    if (computedDensity < pow(10.0, -8.0))
    {
        computedDensity = pow(10.0, -8.0);
    }
    if (computedPressure < pow(10.0, -8.0))
    {
        computedPressure = pow(10.0, -8.0);
    }
    
    vector<double> fluxVector(5);
    
    fluxVector[0] = computedDensity * lorentzFactor * computedYVelocity;
    fluxVector[1] = computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) * (computedXVelocity * computedYVelocity);
    fluxVector[2] = computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) * (computedYVelocity * computedYVelocity) + computedPressure;
    fluxVector[3] = computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) * (computedZVelocity * computedYVelocity);
    fluxVector[4] = ((computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor)) - computedPressure - (computedDensity * lorentzFactor)) *
    computedYVelocity + (computedPressure * computedYVelocity);
    
    return fluxVector;
}

vector<double> RelativisticEulerStateVector::computeYFluxVector(EulerMaterialParameters materialParameters)
{
    return computeYFluxVector(computeConservedVariableVector(materialParameters), materialParameters);
}

double RelativisticEulerStateVector::computeSoundSpeed(EulerMaterialParameters materialParameters)
{
    return RelativisticEulerEquationOfState::computeSoundSpeed(density, pressure, materialParameters);
}

void RelativisticEulerStateVector::setDensity(double newDensity)
{
    density = newDensity;
}

void RelativisticEulerStateVector::setXVelocity(double newXVelocity)
{
    xVelocity = newXVelocity;
}

void RelativisticEulerStateVector::setYVelocity(double newYVelocity)
{
    yVelocity = newYVelocity;
}

void RelativisticEulerStateVector::setZVelocity(double newZVelocity)
{
    zVelocity = newZVelocity;
}

void RelativisticEulerStateVector::setPressure(double newPressure)
{
    pressure = newPressure;
}

double RelativisticEulerStateVector::getDensity()
{
    return density;
}

double RelativisticEulerStateVector::getXVelocity()
{
    return xVelocity;
}

double RelativisticEulerStateVector::getYVelocity()
{
    return yVelocity;
}

double RelativisticEulerStateVector::getZVelocity()
{
    return zVelocity;
}

double RelativisticEulerStateVector::getPressure()
{
    return pressure;
}
