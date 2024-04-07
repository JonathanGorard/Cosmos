#include "BlackHoleEulerStateVector.hpp"

BlackHoleEulerStateVector::BlackHoleEulerStateVector()
{
    density = 1.0;
    xVelocity = 0.0;
    yVelocity = 0.0;
    zVelocity = 0.0;
    pressure = 1.0;
}

BlackHoleEulerStateVector::BlackHoleEulerStateVector(double newDensity, double newXVelocity, double newYVelocity, double newZVelocity, double newPressure)
{
    density = newDensity;
    xVelocity = newXVelocity;
    yVelocity = newYVelocity;
    zVelocity = newZVelocity;
    pressure = newPressure;
}

void BlackHoleEulerStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    density = newPrimitiveVariableVector[0];
    xVelocity = newPrimitiveVariableVector[1];
    yVelocity = newPrimitiveVariableVector[2];
    zVelocity = newPrimitiveVariableVector[3];
    pressure = newPrimitiveVariableVector[4];
}

void BlackHoleEulerStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, double xCoordinate, double yCoordinate, double zCoordinate,
                                                           EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double spatialMetricDeterminant = blackHole.computeSpatialMetricDeterminant(xCoordinate, yCoordinate, zCoordinate);
    
    double relativisticDensity = newConservedVariableVector[0] / sqrt(spatialMetricDeterminant);
    double xMomentum = newConservedVariableVector[1] / sqrt(spatialMetricDeterminant);
    double yMomentum = newConservedVariableVector[2] / sqrt(spatialMetricDeterminant);
    double zMomentum = newConservedVariableVector[3] / sqrt(spatialMetricDeterminant);
    double relativisticEnergy = newConservedVariableVector[4] / sqrt(spatialMetricDeterminant);
    
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

vector<double> BlackHoleEulerStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(5);
    
    primitiveVariableVector[0] = density;
    primitiveVariableVector[1] = xVelocity;
    primitiveVariableVector[2] = yVelocity;
    primitiveVariableVector[3] = zVelocity;
    primitiveVariableVector[4] = pressure;
    
    return primitiveVariableVector;
}

vector<double> BlackHoleEulerStateVector::computeConservedVariableVector(double xCoordinate, double yCoordinate, double zCoordinate,
                                                                         EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<vector<double> > spatialMetricTensor = blackHole.computeSpatialMetricTensor(xCoordinate, yCoordinate, zCoordinate);
    double spatialMetricDeterminant = blackHole.computeSpatialMetricDeterminant(xCoordinate, yCoordinate, zCoordinate);
    
    vector<double> velocityVector(3);
    
    velocityVector[0] = xVelocity;
    velocityVector[1] = yVelocity;
    velocityVector[2] = zVelocity;
    
    double velocitySquared = 0.0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            velocitySquared += spatialMetricTensor[i][j] * velocityVector[i] * velocityVector[j];
        }
    }
    
    double lorentzFactor = 1.0 / sqrt(1.0 - velocitySquared);
    
    if (velocitySquared > 1.0 - pow(10.0, -8.0))
    {
        lorentzFactor = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
    }
    
    double specificEnthalpy = RelativisticEulerEquationOfState::computeSpecificEnthalpy(density, pressure, materialParameters);
    
    vector<double> conservedVariableVector(5);
    
    conservedVariableVector[0] = sqrt(spatialMetricDeterminant) * density * lorentzFactor;
    conservedVariableVector[1] = sqrt(spatialMetricDeterminant) * density * specificEnthalpy * (lorentzFactor * lorentzFactor) * xVelocity;
    conservedVariableVector[2] = sqrt(spatialMetricDeterminant) * density * specificEnthalpy * (lorentzFactor * lorentzFactor) * yVelocity;
    conservedVariableVector[3] = sqrt(spatialMetricDeterminant) * density * specificEnthalpy * (lorentzFactor * lorentzFactor) * zVelocity;
    conservedVariableVector[4] = sqrt(spatialMetricDeterminant) * ((density * specificEnthalpy * (lorentzFactor)) - pressure - (density * lorentzFactor));
    
    return conservedVariableVector;
}

vector<double> BlackHoleEulerStateVector::computeXFluxVector(vector<double> conservedVariableVector, double xCoordinate, double yCoordinate, double zCoordinate,
                                                             EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double spatialMetricDeterminant = blackHole.computeSpatialMetricDeterminant(xCoordinate, yCoordinate, zCoordinate);
    
    double lapseFunction = blackHole.computeLapseFunction(xCoordinate, yCoordinate, zCoordinate);
    vector<double> shiftVector = blackHole.computeShiftVector(xCoordinate, yCoordinate, zCoordinate);
    
    double relativisticDensity = conservedVariableVector[0] / sqrt(spatialMetricDeterminant);
    double xMomentum = conservedVariableVector[1] / sqrt(spatialMetricDeterminant);
    double yMomentum = conservedVariableVector[2] / sqrt(spatialMetricDeterminant);
    double zMomentum = conservedVariableVector[3] / sqrt(spatialMetricDeterminant);
    double relativisticEnergy = conservedVariableVector[4] / sqrt(spatialMetricDeterminant);
    
    double densityConstant = relativisticDensity / sqrt(((relativisticEnergy + relativisticDensity) * (relativisticEnergy + relativisticDensity)) -
                                                        ((xMomentum * xMomentum ) + (yMomentum * yMomentum) + (zMomentum * zMomentum)));
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
        
        if (abs(currentGuess - newGuess) < pow(10.0, -10.0))
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
    
    fluxVector[0] = (lapseFunction * sqrt(spatialMetricDeterminant)) * (computedDensity * lorentzFactor * (computedXVelocity - (shiftVector[0] / lapseFunction)));
    fluxVector[1] = (lapseFunction * sqrt(spatialMetricDeterminant)) * (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) *
                                                                        (computedXVelocity * (computedXVelocity - (shiftVector[0] / lapseFunction))) + computedPressure);
    fluxVector[2] = (lapseFunction * sqrt(spatialMetricDeterminant)) * (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) *
                                                                        (computedYVelocity * (computedXVelocity - (shiftVector[0] / lapseFunction))));
    fluxVector[3] = (lapseFunction * sqrt(spatialMetricDeterminant)) * (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) *
                                                                        (computedZVelocity * (computedXVelocity - (shiftVector[0] / lapseFunction))));
    fluxVector[4] = (lapseFunction * sqrt(spatialMetricDeterminant)) * (((computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor)) - computedPressure -
                                                                         (computedDensity * lorentzFactor)) * (computedXVelocity - (shiftVector[0] / lapseFunction)) +
                                                                        (computedPressure * computedXVelocity));
    
    return fluxVector;
}

vector<double> BlackHoleEulerStateVector::computeXFluxVector(double xCoordinate, double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                             BlackHoleSpacetime blackHole)
{
    return computeXFluxVector(computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole),
                              xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
}

vector<double> BlackHoleEulerStateVector::computeYFluxVector(vector<double> conservedVariableVector, double xCoordinate, double yCoordinate, double zCoordinate,
                                                             EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double spatialMetricDeterminant = blackHole.computeSpatialMetricDeterminant(xCoordinate, yCoordinate, zCoordinate);
    
    double lapseFunction = blackHole.computeLapseFunction(xCoordinate, yCoordinate, zCoordinate);
    vector<double> shiftVector = blackHole.computeShiftVector(xCoordinate, yCoordinate, zCoordinate);
    
    double relativisticDensity = conservedVariableVector[0] / sqrt(spatialMetricDeterminant);
    double xMomentum = conservedVariableVector[1] / sqrt(spatialMetricDeterminant);
    double yMomentum = conservedVariableVector[2] / sqrt(spatialMetricDeterminant);
    double zMomentum = conservedVariableVector[3] / sqrt(spatialMetricDeterminant);
    double relativisticEnergy = conservedVariableVector[4] / sqrt(spatialMetricDeterminant);
    
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
    
    fluxVector[0] = (lapseFunction * sqrt(spatialMetricDeterminant)) * (computedDensity * lorentzFactor * (computedYVelocity - (shiftVector[1] / lapseFunction)));
    fluxVector[1] = (lapseFunction * sqrt(spatialMetricDeterminant)) * (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) *
                                                                        (computedXVelocity * (computedYVelocity - (shiftVector[1] / lapseFunction))));
    fluxVector[2] = (lapseFunction * sqrt(spatialMetricDeterminant)) * (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) *
                                                                        (computedYVelocity * (computedYVelocity - (shiftVector[1] / lapseFunction))) + computedPressure);
    fluxVector[3] = (lapseFunction * sqrt(spatialMetricDeterminant)) * (computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) *
                                                                        (computedZVelocity * (computedYVelocity - (shiftVector[1] / lapseFunction))));
    fluxVector[4] = (lapseFunction * sqrt(spatialMetricDeterminant)) * (((computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor)) - computedPressure -
                                                                         (computedDensity * lorentzFactor)) * (computedYVelocity - (shiftVector[1] / lapseFunction)) +
                                                                        (computedPressure * computedYVelocity));
    
    return fluxVector;
}

vector<double> BlackHoleEulerStateVector::computeYFluxVector(double xCoordinate, double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                             BlackHoleSpacetime blackHole)
{
    return computeYFluxVector(computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole),
                              xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
}

vector<double> BlackHoleEulerStateVector::computeSourceTermVector(vector<double> conservedVariableVector, double xCoordinate, double yCoordinate, double zCoordinate,
                                                                  double cellSpacing, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double spatialMetricDeterminant = blackHole.computeSpatialMetricDeterminant(xCoordinate, yCoordinate, zCoordinate);
    
    double lapseFunction = blackHole.computeLapseFunction(xCoordinate, yCoordinate, zCoordinate);
    vector<double> shiftVector = blackHole.computeShiftVector(xCoordinate, yCoordinate, zCoordinate);
    
    double relativisticDensity = conservedVariableVector[0] / sqrt(spatialMetricDeterminant);
    double xMomentum = conservedVariableVector[1] / sqrt(spatialMetricDeterminant);
    double yMomentum = conservedVariableVector[2] / sqrt(spatialMetricDeterminant);
    double zMomentum = conservedVariableVector[3] / sqrt(spatialMetricDeterminant);
    double relativisticEnergy = conservedVariableVector[4] / sqrt(spatialMetricDeterminant);
    
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
    
    vector<vector<double> > inverseSpacetimeMetricTensor = blackHole.computeInverseSpacetimeMetricTensor(xCoordinate, yCoordinate, zCoordinate);
    
    vector<double> spacetimeVelocityVector(4);
    spacetimeVelocityVector[0] = lorentzFactor / lapseFunction;
    spacetimeVelocityVector[1] = (lapseFunction * spacetimeVelocityVector[0] * computedXVelocity) - (shiftVector[0] * spacetimeVelocityVector[0]);
    spacetimeVelocityVector[2] = (lapseFunction * spacetimeVelocityVector[0] * computedYVelocity) - (shiftVector[1] * spacetimeVelocityVector[0]);
    spacetimeVelocityVector[3] = (lapseFunction * spacetimeVelocityVector[0] * computedZVelocity) - (shiftVector[2] * spacetimeVelocityVector[0]);
    
    vector<vector<double> > stressEnergyTensor(4, vector<double>(4));
    
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            stressEnergyTensor[i][j] = (computedDensity * specificEnthalpy * spacetimeVelocityVector[i] * spacetimeVelocityVector[j]) +
            (computedPressure * inverseSpacetimeMetricTensor[i][j]);
        }
    }
    
    vector<double> lapseFunctionDerivative = blackHole.computeLapseFunctionDerivative(xCoordinate, yCoordinate, zCoordinate, cellSpacing);
    vector<vector<double> > shiftVectorDerivative = blackHole.computeShiftVectorDerivative(xCoordinate, yCoordinate, zCoordinate, cellSpacing);
    
    vector<vector<double> > spatialMetricTensor = blackHole.computeSpatialMetricTensor(xCoordinate, yCoordinate, zCoordinate);
    vector<vector<vector<double> > > spatialMetricTensorDerivative = blackHole.computeSpatialMetricTensorDerivative(xCoordinate, yCoordinate, zCoordinate, cellSpacing);
    vector<vector<double> > extrinsicCurvatureTensor = blackHole.computeExtrinsicCurvatureTensor(xCoordinate, yCoordinate, zCoordinate, cellSpacing);
    
    
    vector<double> sourceTermVector(5);
    
    for (int i = 0; i < 5; i++)
    {
        sourceTermVector[i] = 0.0;
    }
    
    vector<double> velocityCovector(3);
    
    for (int i = 0; i < 3; i++)
    {
        velocityCovector[i] = (spatialMetricTensor[i][0] * computedXVelocity) + (spatialMetricTensor[i][1] * computedYVelocity) +
        (spatialMetricTensor[i][2] * computedZVelocity);
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                sourceTermVector[i + 1] += stressEnergyTensor[0][0] * (0.5 * shiftVector[j] * shiftVector[k] * spatialMetricTensorDerivative[i][j][k]);
                sourceTermVector[i + 1] += stressEnergyTensor[0][j + 1] * shiftVector[k] * spatialMetricTensorDerivative[i][j][k];
                sourceTermVector[i + 1] += 0.5 * (stressEnergyTensor[j + 1][k + 1] * spatialMetricTensorDerivative[i][j][k]);
            }
            
            sourceTermVector[i + 1] += ((computedDensity * specificEnthalpy * (lorentzFactor * lorentzFactor) * velocityCovector[j]) / lapseFunction) *
            shiftVectorDerivative[i][j];
        }
        
        sourceTermVector[i + 1] -= stressEnergyTensor[0][0] * lapseFunction * lapseFunctionDerivative[i];
    }
    
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sourceTermVector[4] += stressEnergyTensor[0][0] * (shiftVector[i] * shiftVector[j] * extrinsicCurvatureTensor[i][j]);
            sourceTermVector[4] += stressEnergyTensor[0][i + 1] * (2.0 * shiftVector[j] * extrinsicCurvatureTensor[i][j]);
            sourceTermVector[4] += stressEnergyTensor[i + 1][j + 1] * extrinsicCurvatureTensor[i][j];
        }
        
        sourceTermVector[4] -= stressEnergyTensor[0][0] * (shiftVector[i] * lapseFunctionDerivative[i]);
        sourceTermVector[4] -= stressEnergyTensor[0][i + 1] * lapseFunctionDerivative[i];
     }
    
    return VectorAlgebra::multiplyVector((lapseFunction * sqrt(spatialMetricDeterminant)), sourceTermVector);
}

vector<double> BlackHoleEulerStateVector::computeSourceTermVector(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing,
                                                                  EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    return computeSourceTermVector(computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole),
                                   xCoordinate, yCoordinate, zCoordinate, cellSpacing, materialParameters, blackHole);
}

double BlackHoleEulerStateVector::computeSoundSpeed(double xCoordinate, double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                    BlackHoleSpacetime blackHole)
{
    vector<double> materialWaveEigenvalues = BlackHoleEulerEigenvalues::computeMaterialWaveEigenvalues(xCoordinate, yCoordinate, zCoordinate, xVelocity, yVelocity,
                                                                                                       zVelocity, materialParameters, blackHole);
    vector<double> fastAcousticWaveEigenvalues = BlackHoleEulerEigenvalues::computeFastAcousticWaveEigenvalues(xCoordinate, yCoordinate, zCoordinate, density,
                                                                                                               pressure, xVelocity, yVelocity, zVelocity,
                                                                                                               materialParameters, blackHole);
    vector<double> slowAcousticWaveEigenvalues = BlackHoleEulerEigenvalues::computeSlowAcousticWaveEigenvalues(xCoordinate, yCoordinate, zCoordinate, density,
                                                                                                               pressure, xVelocity, yVelocity, zVelocity,
                                                                                                               materialParameters, blackHole);
    
    double soundSpeed = 0.0;
    for (int i = 0; i < 3; i++)
    {
        if (materialWaveEigenvalues[i] > soundSpeed)
        {
            soundSpeed = materialWaveEigenvalues[i];
        }
        if (fastAcousticWaveEigenvalues[i] > soundSpeed)
        {
            soundSpeed = fastAcousticWaveEigenvalues[i];
        }
        if (slowAcousticWaveEigenvalues[i] > soundSpeed)
        {
            soundSpeed = slowAcousticWaveEigenvalues[i];
        }
    }
    
    return soundSpeed;
}

void BlackHoleEulerStateVector::setDensity(double newDensity)
{
    density = newDensity;
}

void BlackHoleEulerStateVector::setXVelocity(double newXVelocity)
{
    xVelocity = newXVelocity;
}

void BlackHoleEulerStateVector::setYVelocity(double newYVelocity)
{
    yVelocity = newYVelocity;
}

void BlackHoleEulerStateVector::setZVelocity(double newZVelocity)
{
    zVelocity = newZVelocity;
}

void BlackHoleEulerStateVector::setPressure(double newPressure)
{
    pressure = newPressure;
}

double BlackHoleEulerStateVector::getDensity()
{
    return density;
}

double BlackHoleEulerStateVector::getXVelocity()
{
    return xVelocity;
}

double BlackHoleEulerStateVector::getYVelocity()
{
    return yVelocity;
}

double BlackHoleEulerStateVector::getZVelocity()
{
    return zVelocity;
}

double BlackHoleEulerStateVector::getPressure()
{
    return pressure;
}
