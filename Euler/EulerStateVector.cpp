#include "EulerStateVector.hpp"

EulerStateVector::EulerStateVector()
{
    density = 1.0;
    xVelocity = 0.0;
    yVelocity = 0.0;
    zVelocity = 0.0;
    pressure = 1.0;
}

EulerStateVector::EulerStateVector(double newDensity, double newXVelocity, double newYVelocity, double newZVelocity, double newPressure)
{
    density = newDensity;
    xVelocity = newXVelocity;
    yVelocity = newYVelocity;
    zVelocity = newZVelocity;
    pressure = newPressure;
}

void EulerStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    density = newPrimitiveVariableVector[0];
    xVelocity = newPrimitiveVariableVector[1];
    yVelocity = newPrimitiveVariableVector[2];
    zVelocity = newPrimitiveVariableVector[3];
    pressure = newPrimitiveVariableVector[4];
}

void EulerStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, EulerMaterialParameters materialParameters)
{
    density = newConservedVariableVector[0];
    xVelocity = newConservedVariableVector[1] / density;
    yVelocity = newConservedVariableVector[2] / density;
    zVelocity = newConservedVariableVector[3] / density;
    
    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);
    double specificInternalEnergy = (newConservedVariableVector[4] / density) - (0.5 * velocitySquared);
    pressure = EulerEquationOfState::computePressure(density, specificInternalEnergy, materialParameters);
}

vector<double> EulerStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(5);
    
    primitiveVariableVector[0] = density;
    primitiveVariableVector[1] = xVelocity;
    primitiveVariableVector[2] = yVelocity;
    primitiveVariableVector[3] = zVelocity;
    primitiveVariableVector[4] = pressure;
    
    return primitiveVariableVector;
}

vector<double> EulerStateVector::computeConservedVariableVector(EulerMaterialParameters materialParameters)
{
    vector<double> conservedVariableVector(5);
    
    conservedVariableVector[0] = density;
    conservedVariableVector[1] = density * xVelocity;
    conservedVariableVector[2] = density * yVelocity;
    conservedVariableVector[3] = density * zVelocity;
    conservedVariableVector[4] = computeTotalEnergy(materialParameters);
    
    return conservedVariableVector;
}

vector<double> EulerStateVector::computeXFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters materialParameters)
{
    vector<double> fluxVector(5);
    
    double computedDensity = conservedVariableVector[0];
    double computedXVelocity = conservedVariableVector[1] / computedDensity;
    double computedYVelocity = conservedVariableVector[2] / computedDensity;
    double computedZVelocity = conservedVariableVector[3] / computedDensity;
    
    double velocitySquared = (computedXVelocity * computedXVelocity) + (computedYVelocity * computedYVelocity) + (computedZVelocity * computedZVelocity);
    double totalEnergy = conservedVariableVector[4] / computedDensity;
    double specificInternalEnergy = totalEnergy - (0.5 * velocitySquared);
    double computedPressure = EulerEquationOfState::computePressure(computedDensity, specificInternalEnergy, materialParameters);
    
    fluxVector[0] = computedDensity * computedXVelocity;
    fluxVector[1] = (computedDensity * (computedXVelocity * computedXVelocity)) + computedPressure;
    fluxVector[2] = computedDensity * (computedXVelocity * computedYVelocity);
    fluxVector[3] = computedDensity * (computedXVelocity * computedZVelocity);
    fluxVector[4] = (computedDensity * (computedXVelocity * totalEnergy)) + (computedXVelocity * computedPressure);
    
    return fluxVector;
}

vector<double> EulerStateVector::computeXFluxVector(EulerMaterialParameters materialParameters)
{
    return computeXFluxVector(computeConservedVariableVector(materialParameters), materialParameters);
}

vector<double> EulerStateVector::computeYFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters materialParameters)
{
    vector<double> fluxVector(5);
    
    double computedDensity = conservedVariableVector[0];
    double computedXVelocity = conservedVariableVector[1] / computedDensity;
    double computedYVelocity = conservedVariableVector[2] / computedDensity;
    double computedZVelocity = conservedVariableVector[3] / computedDensity;
    
    double velocitySquared = (computedXVelocity * computedXVelocity) + (computedYVelocity * computedYVelocity) + (computedZVelocity * computedZVelocity);
    double totalEnergy = conservedVariableVector[4] / computedDensity;
    double specificInternalEnergy = totalEnergy - (0.5 * velocitySquared);
    double computedPressure = EulerEquationOfState::computePressure(computedDensity, specificInternalEnergy, materialParameters);
    
    fluxVector[0] = computedDensity * computedYVelocity;
    fluxVector[1] = computedDensity * (computedYVelocity * computedXVelocity);
    fluxVector[2] = (computedDensity * (computedYVelocity * computedYVelocity)) + computedPressure;
    fluxVector[3] = computedDensity * (computedYVelocity * computedZVelocity);
    fluxVector[4] = (computedDensity * (computedYVelocity * totalEnergy)) + (computedYVelocity * computedPressure);
    
    return fluxVector;
}

vector<double> EulerStateVector::computeYFluxVector(EulerMaterialParameters materialParameters)
{
    return computeYFluxVector(computeConservedVariableVector(materialParameters), materialParameters);
}

double EulerStateVector::computeSpecificInternalEnergy(EulerMaterialParameters materialParameters)
{
    return EulerEquationOfState::computeSpecificInternalEnergy(density, pressure, materialParameters);
}

double EulerStateVector::computeTotalEnergy(EulerMaterialParameters materialParameters)
{
    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);
    
    return density * ((0.5 * velocitySquared) + computeSpecificInternalEnergy(materialParameters));
}

double EulerStateVector::computeSoundSpeed(EulerMaterialParameters materialParameters)
{
    return EulerEquationOfState::computeSoundSpeed(density, pressure, materialParameters);
}

void EulerStateVector::setDensity(double newDensity)
{
    density = newDensity;
}

void EulerStateVector::setXVelocity(double newXVelocity)
{
    xVelocity = newXVelocity;
}

void EulerStateVector::setYVelocity(double newYVelocity)
{
    yVelocity = newYVelocity;
}

void EulerStateVector::setZVelocity(double newZVelocity)
{
    zVelocity = newZVelocity;
}

void EulerStateVector::setPressure(double newPressure)
{
    pressure = newPressure;
}

double EulerStateVector::getDensity()
{
    return density;
}

double EulerStateVector::getXVelocity()
{
    return xVelocity;
}

double EulerStateVector::getYVelocity()
{
    return yVelocity;
}

double EulerStateVector::getZVelocity()
{
    return zVelocity;
}

double EulerStateVector::getPressure()
{
    return pressure;
}
