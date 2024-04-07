#include "RelativisticEulerSolvers.hpp"

vector<RelativisticEulerStateVector> RelativisticEulerSolvers::insertBoundaryCells(vector<RelativisticEulerStateVector> & currentCells, int boundarySize)
{
    long cellCount = currentCells.size();
    vector<RelativisticEulerStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));
    
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

vector<vector<RelativisticEulerStateVector> > RelativisticEulerSolvers::insertBoundaryCells2D(vector<vector<RelativisticEulerStateVector> > & currentCells,
                                                                                              int boundarySize)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<RelativisticEulerStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize),
                                                                           vector<RelativisticEulerStateVector>(columnCount + (2 * boundarySize)));
    
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

double RelativisticEulerSolvers::computeMaximumWaveSpeed(vector<RelativisticEulerStateVector> & currentCells, EulerMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    long cellCount = currentCells.size();
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getXVelocity()) + abs(currentCells[i].computeSoundSpeed(materialParameters));
        
        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }
    
    return maximumWaveSpeed;
}

double RelativisticEulerSolvers::computeMaximumWaveSpeed2D(vector<vector<RelativisticEulerStateVector> > & currentCells, EulerMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getXVelocity()), abs(currentCells[i][j].getYVelocity())) +
            abs(currentCells[i][j].computeSoundSpeed(materialParameters));
            
            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }
    
    return maximumWaveSpeed;
}

double RelativisticEulerSolvers::computeStableTimeStep(vector<RelativisticEulerStateVector> currentCells, double cellSpacing, double CFLCoefficient,
                                                       double currentTime, double finalTime, int currentIteration, EulerMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, materialParameters));
    
    return EulerSolvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double RelativisticEulerSolvers::computeStableTimeStep2D(vector<vector<RelativisticEulerStateVector> > currentCells, double cellSpacing, double CFLCoefficient,
                                                         double currentTime, double finalTime, int currentIteration, EulerMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, materialParameters));
    
    return EulerSolvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

RelativisticEulerStateVector RelativisticEulerSolvers::evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue,
                                                                                 vector<double> evolutionVector, int side, EulerMaterialParameters materialParameters)
{
    RelativisticEulerStateVector evolvedStateVector;
    
    if (side == 0)
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(leftExtrapolatedValue, evolutionVector), materialParameters);
    }
    else
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(rightExtrapolatedValue, evolutionVector), materialParameters);
    }
    
    return evolvedStateVector;
}

RelativisticEulerStateVector RelativisticEulerSolvers::evolveStateByHalfXTimeStep(RelativisticEulerStateVector leftStateVector,
                                                                                  RelativisticEulerStateVector middleStateVector,
                                                                                  RelativisticEulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                                                  double bias, int slopeLimiter, int side, EulerMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, materialParameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));
    
    vector<double> leftFluxVector = RelativisticEulerStateVector::computeXFluxVector(leftExtrapolatedValue, materialParameters);
    vector<double> rightFluxVector = RelativisticEulerStateVector::computeXFluxVector(rightExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = EulerSolvers::computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);
    
    return evolveStateByHalfTimeStep(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, materialParameters);
}

RelativisticEulerStateVector RelativisticEulerSolvers::evolveStateByHalfYTimeStep(RelativisticEulerStateVector topStateVector,
                                                                                  RelativisticEulerStateVector middleStateVector,
                                                                                  RelativisticEulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                                                  double bias, int slopeLimiter, int side, EulerMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, materialParameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));
    
    vector<double> topFluxVector = RelativisticEulerStateVector::computeYFluxVector(topExtrapolatedValue, materialParameters);
    vector<double> bottomFluxVector = RelativisticEulerStateVector::computeYFluxVector(bottomExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = EulerSolvers::computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);
    
    return evolveStateByHalfTimeStep(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, materialParameters);
}
