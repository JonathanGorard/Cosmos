#include "EulerSolvers.hpp"

vector<EulerStateVector> EulerSolvers::insertBoundaryCells(vector<EulerStateVector> & currentCells, int boundarySize)
{
    long cellCount = currentCells.size();
    vector<EulerStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));
    
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

vector<vector<EulerStateVector> > EulerSolvers::insertBoundaryCells2D(vector<vector<EulerStateVector> > & currentCells, int boundarySize)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<EulerStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<EulerStateVector>(columnCount + (2 * boundarySize)));
    
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

double EulerSolvers::computeMaximumWaveSpeed(vector<EulerStateVector> & currentCells, EulerMaterialParameters materialParameters)
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

double EulerSolvers::computeMaximumWaveSpeed2D(vector<vector<EulerStateVector> > & currentCells, EulerMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getXVelocity()), abs(currentCells[i][j].getYVelocity()) +
                                   abs(currentCells[i][j].computeSoundSpeed(materialParameters)));
            
            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }
    
    return maximumWaveSpeed;
}

double EulerSolvers::computeStableTimeStep(double timeStep, double currentTime, double finalTime, int currentIteration)
{
    double newTimeStep = timeStep;
    
    if (currentIteration <= 5)
    {
        newTimeStep *= 0.2;
    }
    
    if ((currentTime + newTimeStep) > finalTime)
    {
        newTimeStep = finalTime - currentTime;
    }
    
    return newTimeStep;
}

double EulerSolvers::computeStableTimeStep(vector<EulerStateVector> currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                           int currentIteration, EulerMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, materialParameters));
    
    return computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double EulerSolvers::computeStableTimeStep2D(vector<vector<EulerStateVector> > currentCells, double cellSpacing, double CFLCoefficient, double currentTime,
                                             double finalTime, int currentIteration, EulerMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, materialParameters));
    
    return computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

vector<double> EulerSolvers::computeFractionalEvolutionVector(double stepFraction, vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing,
                                                              double timeStep)
{
    return VectorAlgebra::multiplyVector(stepFraction * (timeStep / cellSpacing), VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector));
}

vector<double> EulerSolvers::computeEvolutionVector(vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing, double timeStep)
{
    return computeFractionalEvolutionVector(0.5, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

EulerStateVector EulerSolvers::evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector,
                                                         int side, EulerMaterialParameters materialParameters)
{
    EulerStateVector evolvedStateVector;
    
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

EulerStateVector EulerSolvers::evolveStateByHalfXTimeStep(EulerStateVector leftStateVector, EulerStateVector middleStateVector, EulerStateVector rightStateVector,
                                                          double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                          EulerMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, materialParameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));
    
    vector<double> leftFluxVector = EulerStateVector::computeXFluxVector(leftExtrapolatedValue, materialParameters);
    vector<double> rightFluxVector = EulerStateVector::computeXFluxVector(rightExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);
    
    return evolveStateByHalfTimeStep(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, materialParameters);
}

EulerStateVector EulerSolvers::evolveStateByHalfYTimeStep(EulerStateVector topStateVector, EulerStateVector middleStateVector, EulerStateVector bottomStateVector,
                                                          double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                          EulerMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, materialParameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));
    
    vector<double> topFluxVector = EulerStateVector::computeYFluxVector(topExtrapolatedValue, materialParameters);
    vector<double> bottomFluxVector = EulerStateVector::computeYFluxVector(bottomExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);
    
    return evolveStateByHalfTimeStep(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, materialParameters);
}

void EulerSolvers::outputStatus(int currentIteration, double currentTime, double timeStep)
{
    cout << "Iteration = " << currentIteration << "; Time = " << currentTime << "; Timestep = " << timeStep << endl;
}
