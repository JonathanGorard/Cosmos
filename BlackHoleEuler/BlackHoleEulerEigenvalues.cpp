#include "BlackHoleEulerEigenvalues.hpp"

vector<double> BlackHoleEulerEigenvalues::computeMaterialWaveEigenvalues(double xCoordinate, double yCoordinate, double zCoordinate, double xVelocity, double yVelocity,
                                                                         double zVelocity, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    double lapseFunction = blackHole.computeLapseFunction(xCoordinate, yCoordinate, zCoordinate);
    vector<double> shiftVector = blackHole.computeShiftVector(xCoordinate, yCoordinate, zCoordinate);
    
    vector<double> materialWaveEigenvalues(3);
    
    materialWaveEigenvalues[0] = (lapseFunction * xVelocity) - shiftVector[0];
    materialWaveEigenvalues[1] = (lapseFunction * yVelocity) - shiftVector[1];
    materialWaveEigenvalues[2] = (lapseFunction * zVelocity) - shiftVector[2];
    
    return materialWaveEigenvalues;
}

vector<double> BlackHoleEulerEigenvalues::computeFastAcousticWaveEigenvalues(double xCoordinate, double yCoordinate, double zCoordinate, double density,
                                                                             double pressure, double xVelocity, double yVelocity, double zVelocity,
                                                                             EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    double lapseFunction = blackHole.computeLapseFunction(xCoordinate, yCoordinate, zCoordinate);
    vector<double> shiftVector = blackHole.computeShiftVector(xCoordinate, yCoordinate, zCoordinate);
    
    vector<vector<double> > spatialMetricTensor = blackHole.computeSpatialMetricTensor(xCoordinate, yCoordinate, zCoordinate);
    vector<vector<double> > inverseSpatialMetricTensor = blackHole.computeInverseSpatialMetricTensor(xCoordinate, yCoordinate, zCoordinate);
    
    double soundSpeed = RelativisticEulerEquationOfState::computeSoundSpeed(density, pressure, materialParameters);
    
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
    
    vector<double> fastAcousticWaveEigenvalues(3);
    
    for (int i = 0; i < 3; i++)
    {
        fastAcousticWaveEigenvalues[i] = (lapseFunction / (1.0 - (velocitySquared * (soundSpeed * soundSpeed)))) *
        ((velocityVector[i] * (1.0 - (soundSpeed * soundSpeed))) +
         (soundSpeed * sqrt((1.0 - velocitySquared) * (inverseSpatialMetricTensor[i][i] * (1.0 - (velocitySquared * (soundSpeed * soundSpeed))) -
                                                       (velocityVector[i] * velocityVector[i]) * (1.0 - (soundSpeed * soundSpeed)))))) - shiftVector[i];
    }
    
    return fastAcousticWaveEigenvalues;
}

vector<double> BlackHoleEulerEigenvalues::computeSlowAcousticWaveEigenvalues(double xCoordinate, double yCoordinate, double zCoordinate, double density,
                                                                             double pressure, double xVelocity, double yVelocity, double zVelocity,
                                                                             EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    double lapseFunction = blackHole.computeLapseFunction(xCoordinate, yCoordinate, zCoordinate);
    vector<double> shiftVector = blackHole.computeShiftVector(xCoordinate, yCoordinate, zCoordinate);
    
    vector<vector<double> > spatialMetricTensor = blackHole.computeSpatialMetricTensor(xCoordinate, yCoordinate, zCoordinate);
    vector<vector<double> > inverseSpatialMetricTensor = blackHole.computeInverseSpatialMetricTensor(xCoordinate, yCoordinate, zCoordinate);
    
    double soundSpeed = RelativisticEulerEquationOfState::computeSoundSpeed(density, pressure, materialParameters);
    
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
    
    vector<double> slowAcousticWaveEigenvalues(3);
    
    for (int i = 0; i < 3; i++)
    {
        slowAcousticWaveEigenvalues[i] = (lapseFunction / (1.0 - (velocitySquared * (soundSpeed * soundSpeed)))) *
        ((velocityVector[i] * (1.0 - (soundSpeed * soundSpeed))) -
         (soundSpeed * sqrt((1.0 - velocitySquared) * (inverseSpatialMetricTensor[i][i] * (1.0 - (velocitySquared * (soundSpeed * soundSpeed))) -
                                                       (velocityVector[i] * velocityVector[i]) * (1.0 - (soundSpeed * soundSpeed)))))) - shiftVector[i];
    }
    
    return slowAcousticWaveEigenvalues;
}
