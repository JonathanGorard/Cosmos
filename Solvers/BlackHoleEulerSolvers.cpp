#include "BlackHoleEulerSolvers.hpp"

vector<BlackHoleEulerStateVector> BlackHoleEulerSolvers::insertBoundaryCells(vector<BlackHoleEulerStateVector> & currentCells, int boundarySize)
{
    long cellCount = currentCells.size();
    vector<BlackHoleEulerStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));
    
    if (boundarySize == 1)
    {
        currentCellsWithBoundary[0] = currentCells[0];
        currentCellsWithBoundary[cellCount + 1] = currentCells[cellCount - 1];
    }
    else if (boundarySize == 2)
    {
        currentCellsWithBoundary[0] = currentCells[1];
        currentCellsWithBoundary[1] = currentCells[0];
        currentCellsWithBoundary[cellCount + 2] = currentCells[cellCount - 1];
        currentCellsWithBoundary[cellCount + 3] = currentCells[cellCount - 2];
    }
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        currentCellsWithBoundary[i + boundarySize] = currentCells[i];
    }
    
    return currentCellsWithBoundary;
}

vector<vector<BlackHoleEulerStateVector> > BlackHoleEulerSolvers::insertBoundaryCells2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, int boundarySize)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<BlackHoleEulerStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize),
                                                                        vector<BlackHoleEulerStateVector>(columnCount + (2 * boundarySize)));
    
    if (boundarySize == 1)
    {
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            currentCellsWithBoundary[i + 1][0] = currentCells[i][0];
            currentCellsWithBoundary[i + 1][columnCount + 1] = currentCells[i][columnCount - 1];
        }
        
#pragma omp parallel for
        for (int i = 0; i < columnCount; i++)
        {
            currentCellsWithBoundary[0][i + 1] = currentCells[0][i];
            currentCellsWithBoundary[rowCount + 1][i + 1] = currentCells[rowCount - 1][i];
        }
    }
    else if (boundarySize == 2)
    {
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            currentCellsWithBoundary[i + 2][0] = currentCells[i][1];
            currentCellsWithBoundary[i + 2][1] = currentCells[i][0];
            
            currentCellsWithBoundary[i + 2][columnCount + 2] = currentCells[i][columnCount - 1];
            currentCellsWithBoundary[i + 2][columnCount + 3] = currentCells[i][columnCount - 2];
        }
        
#pragma omp parallel for
        for (int i = 0; i < columnCount; i++)
        {
            currentCellsWithBoundary[0][i + 2] = currentCells[1][i];
            currentCellsWithBoundary[1][i + 2] = currentCells[0][i];
            
            currentCellsWithBoundary[rowCount + 2][i + 2] = currentCells[rowCount - 1][i];
            currentCellsWithBoundary[rowCount + 3][i + 2] = currentCells[rowCount - 2][i];
        }
    }
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            currentCellsWithBoundary[i + boundarySize][j + boundarySize] = currentCells[i][j];
        }
    }
    
    return currentCellsWithBoundary;
}

double BlackHoleEulerSolvers::computeMaximumWaveSpeed(vector<BlackHoleEulerStateVector> & currentCells, EulerMaterialParameters materialParameters,
                                                      BlackHoleSpacetime blackHole)
{
    double maximumWaveSpeed = 0.0;
    long cellCount = currentCells.size();
    double cellSpacing = 1.0 / cellCount;
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        if (!blackHole.inExcisionRegion(i * cellSpacing, 0.0, 0.0))
        {
            double waveSpeed = abs(currentCells[i].computeSoundSpeed(i * cellSpacing, 0.0, 0.0, materialParameters, blackHole));
            
            if(waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }
    
    return maximumWaveSpeed;
}

double BlackHoleEulerSolvers::computeMaximumWaveSpeed2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, EulerMaterialParameters materialParameters,
                                                        BlackHoleSpacetime blackHole)
{
    double maximumWaveSpeed = 0.0;
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    
    double rowCellSpacing = 1.0 / rowCount;
    double columnCellSpacing = 1.0 / columnCount;
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (!blackHole.inExcisionRegion(j * columnCellSpacing, i * rowCellSpacing, 0.0))
            {
                double waveSpeed = abs(currentCells[i][j].computeSoundSpeed(i * rowCellSpacing, j * columnCellSpacing, 0.0, materialParameters, blackHole));
                
                if (waveSpeed > maximumWaveSpeed)
                {
                    maximumWaveSpeed = waveSpeed;
                }
            }
        }
    }
    
    return maximumWaveSpeed;
}

double BlackHoleEulerSolvers::computeStableTimeStep(vector<BlackHoleEulerStateVector> currentCells, double cellSpacing, double CFLCoefficient, double currentTime,
                                                    double finalTime, int currentIteration, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, materialParameters, blackHole));
    
    return EulerSolvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double BlackHoleEulerSolvers::computeStableTimeStep2D(vector<vector<BlackHoleEulerStateVector> > currentCells, double cellSpacing, double CFLCoefficient,
                                                      double currentTime, double finalTime, int currentIteration, EulerMaterialParameters materialParameters,
                                                      BlackHoleSpacetime blackHole)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, materialParameters, blackHole));
    
    return EulerSolvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

BlackHoleEulerStateVector BlackHoleEulerSolvers::evolveStateByFractionalTimeStep(vector<double> middleConservedVariableVector,
                                                                                 vector<double> conservedVariableVectorEvolution, double xCoordinate,
                                                                                 double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                                                 BlackHoleSpacetime blackHole)
{
    BlackHoleEulerStateVector evolvedStateVector;
    evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(middleConservedVariableVector, conservedVariableVectorEvolution), xCoordinate,
                                                  yCoordinate, zCoordinate, materialParameters, blackHole);
    
    return evolvedStateVector;
}

BlackHoleEulerStateVector BlackHoleEulerSolvers::evolveStateByFractionalXTimeStep(double stepFraction, BlackHoleEulerStateVector leftStateVector,
                                                                                  BlackHoleEulerStateVector middleStateVector,
                                                                                  BlackHoleEulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                                                  double bias, int slopeLimiter, double xCoordinate, double yCoordinate,
                                                                                  double zCoordinate, EulerMaterialParameters materialParameters,
                                                                                  BlackHoleSpacetime blackHole)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVectorX(leftStateVector, middleStateVector, rightStateVector, cellSpacing, bias, slopeLimiter,
                                                                    xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate, materialParameters,
                                                                                                    blackHole);
    vector<double> leftConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    
    vector<double> leftFluxVector = BlackHoleEulerStateVector::computeXFluxVector(leftConservedVariableVector, xCoordinate - cellSpacing, yCoordinate,
                                                                                  zCoordinate, materialParameters, blackHole);
    vector<double> rightFluxVector = BlackHoleEulerStateVector::computeXFluxVector(rightConservedVariableVector, xCoordinate + cellSpacing, yCoordinate,
                                                                                   zCoordinate, materialParameters, blackHole);
    vector<double> evolutionVector = EulerSolvers::computeFractionalEvolutionVector(stepFraction, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
        
    return evolveStateByFractionalTimeStep(middleConservedVariableVector, evolutionVector, xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
}

BlackHoleEulerStateVector BlackHoleEulerSolvers::evolveStateByFractionalYTimeStep(double stepFraction, BlackHoleEulerStateVector topStateVector,
                                                                                  BlackHoleEulerStateVector middleStateVector,
                                                                                  BlackHoleEulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                                                  double bias, int slopeLimiter, double xCoordinate, double yCoordinate,
                                                                                  double zCoordinate, EulerMaterialParameters materialParameters,
                                                                                  BlackHoleSpacetime blackHole)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVectorY(topStateVector, middleStateVector, bottomStateVector, cellSpacing, bias, slopeLimiter,
                                                                    xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate, materialParameters,
                                                                                                    blackHole);
    vector<double> topConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    
    vector<double> topFluxVector = BlackHoleEulerStateVector::computeYFluxVector(topConservedVariableVector, xCoordinate, yCoordinate - cellSpacing,
                                                                                 zCoordinate, materialParameters, blackHole);
    vector<double> bottomFluxVector = BlackHoleEulerStateVector::computeYFluxVector(bottomConservedVariableVector, xCoordinate, yCoordinate + cellSpacing,
                                                                                    zCoordinate, materialParameters, blackHole);
    vector<double> evolutionVector = EulerSolvers::computeFractionalEvolutionVector(stepFraction, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
    
    return evolveStateByFractionalTimeStep(middleConservedVariableVector, evolutionVector, xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
}

BlackHoleEulerStateVector BlackHoleEulerSolvers::evolveStateByHalfXTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue,
                                                                            vector<double> evolutionVector, int side, double cellSpacing, double xCoordinate,
                                                                            double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                                            BlackHoleSpacetime blackHole)
{
    BlackHoleEulerStateVector evolvedStateVector;
    
    if (side == 0)
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(leftExtrapolatedValue, evolutionVector), xCoordinate - (cellSpacing * 0.5),
                                                      yCoordinate, zCoordinate, materialParameters, blackHole);
    }
    else
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(rightExtrapolatedValue, evolutionVector), xCoordinate + (cellSpacing * 0.5),
                                                      yCoordinate, zCoordinate, materialParameters, blackHole);
    }
    
    return evolvedStateVector;
}

BlackHoleEulerStateVector BlackHoleEulerSolvers::evolveStateByHalfXTimeStep(BlackHoleEulerStateVector leftStateVector, BlackHoleEulerStateVector middleStateVector,
                                                                            BlackHoleEulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                                            double bias, int slopeLimiter, int side, double xCoordinate, double yCoordinate,
                                                                            double zCoordinate, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVectorX(leftStateVector, middleStateVector, rightStateVector, cellSpacing, bias, slopeLimiter, xCoordinate,
                                                                    yCoordinate, zCoordinate, materialParameters, blackHole);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate,
                                                                                                                           materialParameters, blackHole),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate,
                                                                                                                       materialParameters, blackHole),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));
    
    vector<double> leftFluxVector = BlackHoleEulerStateVector::computeXFluxVector(leftExtrapolatedValue, xCoordinate - (cellSpacing * 0.5), yCoordinate,
                                                                                  zCoordinate, materialParameters, blackHole);
    vector<double> rightFluxVector = BlackHoleEulerStateVector::computeXFluxVector(rightExtrapolatedValue, xCoordinate + (cellSpacing * 0.5), yCoordinate,
                                                                                   zCoordinate, materialParameters, blackHole);
    vector<double> evolutionVector = EulerSolvers::computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);
    
    return evolveStateByHalfXTimeStep(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, cellSpacing, xCoordinate, yCoordinate, zCoordinate,
                                      materialParameters, blackHole);
}

BlackHoleEulerStateVector BlackHoleEulerSolvers::evolveStateByHalfYTimeStep(vector<double> topExtrapolatedValue, vector<double> bottomExtrapolatedValue,
                                                                            vector<double> evolutionVector, int side, double cellSpacing, double xCoordinate,
                                                                            double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                                            BlackHoleSpacetime blackHole)
{
    BlackHoleEulerStateVector evolvedStateVector;
    
    if (side == 0)
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(topExtrapolatedValue, evolutionVector), xCoordinate,
                                                      yCoordinate - (cellSpacing * 0.5), zCoordinate, materialParameters, blackHole);
    }
    else
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(bottomExtrapolatedValue, evolutionVector), xCoordinate,
                                                      yCoordinate + (cellSpacing * 0.5), zCoordinate, materialParameters, blackHole);
    }
    
    return evolvedStateVector;
}

BlackHoleEulerStateVector BlackHoleEulerSolvers::evolveStateByHalfYTimeStep(BlackHoleEulerStateVector topStateVector, BlackHoleEulerStateVector middleStateVector,
                                                                            BlackHoleEulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                                            double bias, int slopeLimiter, int side, double xCoordinate, double yCoordinate,
                                                                            double zCoordinate, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVectorY(topStateVector, middleStateVector, bottomStateVector, cellSpacing, bias, slopeLimiter,
                                                                    xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate,
                                                                                                                          materialParameters, blackHole),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate,
                                                                                                                        materialParameters, blackHole),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));
    
    vector<double> topFluxVector = BlackHoleEulerStateVector::computeYFluxVector(topExtrapolatedValue, xCoordinate, yCoordinate - (cellSpacing * 0.5),
                                                                                 zCoordinate, materialParameters, blackHole);
    vector<double> bottomFluxVector = BlackHoleEulerStateVector::computeYFluxVector(bottomExtrapolatedValue, xCoordinate, yCoordinate + (cellSpacing * 0.5),
                                                                                    zCoordinate, materialParameters, blackHole);
    vector<double> evolutionVector = EulerSolvers::computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);
    
    return evolveStateByHalfYTimeStep(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, cellSpacing, xCoordinate, yCoordinate, zCoordinate,
                                      materialParameters, blackHole);
}
