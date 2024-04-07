#include "BlackHoleEulerSecondOrderSolver.hpp"

vector<double> BlackHoleEulerSecondOrderSolver::computeXSLICFlux(BlackHoleEulerStateVector leftLeftStateVector, BlackHoleEulerStateVector leftStateVector,
                                                                 BlackHoleEulerStateVector rightStateVector, BlackHoleEulerStateVector rightRightStateVector,
                                                                 double cellSpacing, double timeStep, double bias, int slopeLimiter, double xCoordinate,
                                                                 double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                                 BlackHoleSpacetime blackHole)
{
    BlackHoleEulerStateVector evolvedRightStateVector = BlackHoleEulerSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector,
                                                                                                          cellSpacing, timeStep, bias, slopeLimiter, 1,
                                                                                                          xCoordinate - (cellSpacing * 0.5), yCoordinate, zCoordinate,
                                                                                                          materialParameters, blackHole);
    BlackHoleEulerStateVector evolvedLeftStateVector = BlackHoleEulerSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector,
                                                                                                         cellSpacing, timeStep, bias, slopeLimiter, 0,
                                                                                                         xCoordinate + (cellSpacing * 0.5), yCoordinate, zCoordinate,
                                                                                                         materialParameters, blackHole);
    
    return BlackHoleEulerFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, xCoordinate, yCoordinate, zCoordinate, cellSpacing,
                                                             timeStep, materialParameters, blackHole);
}

vector<double>  BlackHoleEulerSecondOrderSolver::computeYSLICFlux(BlackHoleEulerStateVector topTopStateVector, BlackHoleEulerStateVector topStateVector,
                                                                  BlackHoleEulerStateVector bottomStateVector, BlackHoleEulerStateVector bottomBottomStateVector,
                                                                  double cellSpacing, double timeStep, double bias, int slopeLimiter, double xCoordinate,
                                                                  double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                                  BlackHoleSpacetime blackHole)
{
    BlackHoleEulerStateVector evolvedBottomStateVector = BlackHoleEulerSolvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector,
                                                                                                           cellSpacing, timeStep, bias, slopeLimiter, 1,
                                                                                                           xCoordinate, yCoordinate - (cellSpacing * 0.5), zCoordinate,
                                                                                                           materialParameters, blackHole);
    BlackHoleEulerStateVector evolvedTopStateVector = BlackHoleEulerSolvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector,
                                                                                                        cellSpacing, timeStep, bias, slopeLimiter, 0,
                                                                                                        xCoordinate, yCoordinate + (cellSpacing * 0.5), zCoordinate,
                                                                                                        materialParameters, blackHole);
    
    return BlackHoleEulerFirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, xCoordinate, yCoordinate, zCoordinate, cellSpacing,
                                                             timeStep, materialParameters, blackHole);
}

void BlackHoleEulerSecondOrderSolver::computeSLICTimeStep(vector<BlackHoleEulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias,
                                                          int slopeLimiter, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long cellCount = currentCells.size();
    vector<BlackHoleEulerStateVector> currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        if (!blackHole.inExcisionRegion(i * cellSpacing, 0.0, 0.0) && !blackHole.inExcisionRegion((i - 1) * cellSpacing, 0.0, 0.0) &&
            !blackHole.inExcisionRegion((i - 2) * cellSpacing, 0.0, 0.0) && !blackHole.inExcisionRegion((i + 1) * cellSpacing, 0.0, 0.0) &&
            !blackHole.inExcisionRegion((i + 2) * cellSpacing, 0.0, 0.0))
        {
            vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(i * cellSpacing, 0.0, 0.0, materialParameters, blackHole);
            
            vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2],
                                                             currentCellsWithBoundary[i + 3], cellSpacing, timeStep, bias, slopeLimiter, (i - 0.5) * cellSpacing, 0.0,
                                                             0.0, materialParameters, blackHole);
            vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                                                              currentCellsWithBoundary[i + 4], cellSpacing, timeStep, bias, slopeLimiter, (i + 0.5) * cellSpacing, 0.0,
                                                              0.0, materialParameters, blackHole);
            
            currentCells[i].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing,
                                                                                                 timeStep), i * cellSpacing, 0.0, 0.0, materialParameters, blackHole);
        }
    }
}

void BlackHoleEulerSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                             double bias, int slopeLimiter, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<BlackHoleEulerStateVector> > currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (!blackHole.inExcisionRegion(j * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j - 1) * cellSpacing, i * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion((j - 2) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j + 1) * cellSpacing, i * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion((j + 2) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i - 1) * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion(j * cellSpacing, (i - 2) * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i + 1) * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion(j * cellSpacing, (i + 2) * cellSpacing, 0.0))
            {
                vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(j * cellSpacing, i * cellSpacing, 0.0, materialParameters,
                                                                                                           blackHole);
                
                vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 2][j], currentCellsWithBoundary[i + 2][j + 1],
                                                                 currentCellsWithBoundary[i + 2][j + 2], currentCellsWithBoundary[i + 2][j + 3],
                                                                 cellSpacing, timeStep, bias, slopeLimiter, (j - 0.5) * cellSpacing, i * cellSpacing, 0.0,
                                                                 materialParameters, blackHole);
                vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 2][j + 1], currentCellsWithBoundary[i + 2][j + 2],
                                                                  currentCellsWithBoundary[i + 2][j + 3], currentCellsWithBoundary[i + 2][j + 4],
                                                                  cellSpacing, timeStep, bias, slopeLimiter, (j + 0.5) * cellSpacing, i * cellSpacing, 0.0,
                                                                  materialParameters, blackHole);
                
                currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                        cellSpacing, timeStep), j * cellSpacing, i * cellSpacing,
                                                              0.0, materialParameters, blackHole);
            }
        }
    }
}

void BlackHoleEulerSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                             double bias, int slopeLimiter, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<BlackHoleEulerStateVector> > currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (!blackHole.inExcisionRegion(j * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j - 1) * cellSpacing, i * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion((j - 2) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j + 1) * cellSpacing, i * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion((j + 2) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i - 1) * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion(j * cellSpacing, (i - 2) * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i + 1) * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion(j * cellSpacing, (i + 2) * cellSpacing, 0.0))
            {
                vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(j * cellSpacing, i * cellSpacing, 0.0, materialParameters,
                                                                                                           blackHole);
                
                vector<double> topFluxVector = computeYSLICFlux(currentCellsWithBoundary[i][j + 2], currentCellsWithBoundary[i + 1][j + 2],
                                                                currentCellsWithBoundary[i + 2][j + 2], currentCellsWithBoundary[i + 3][j + 2],
                                                                cellSpacing, timeStep, bias, slopeLimiter, j * cellSpacing, (i - 0.5) * cellSpacing, 0.0,
                                                                materialParameters, blackHole);
                vector<double> bottomFluxVector = computeYSLICFlux(currentCellsWithBoundary[i + 1][j + 2], currentCellsWithBoundary[i + 2][j + 2],
                                                                   currentCellsWithBoundary[i + 3][j + 2], currentCellsWithBoundary[i + 4][j + 2],
                                                                   cellSpacing, timeStep, bias, slopeLimiter, j * cellSpacing, (i + 0.5) * cellSpacing, 0.0,
                                                                   materialParameters, blackHole);
                
                currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector,
                                                                                                        cellSpacing, timeStep), j * cellSpacing, i * cellSpacing,
                                                              0.0, materialParameters, blackHole);
            }
        }
    }
}

vector<BlackHoleEulerStateVector> BlackHoleEulerSecondOrderSolver::solve(vector<BlackHoleEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient,
                                                                         double finalTime, double bias, int slopeLimiter, int subcyclingIterations,
                                                                         EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<BlackHoleEulerStateVector> currentCells = initialCells;
    
    while (currentTime < finalTime)
    {
        double timeStep = BlackHoleEulerSolvers::computeStableTimeStep(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration,
                                                                       materialParameters, blackHole);
        
        for (int i = 0; i < subcyclingIterations; i++)
        {
            BlackHoleEulerForcingSolver::computeRungeKuttaTimeStep(currentCells, cellSpacing, 0.5 * (timeStep / subcyclingIterations), 0.0, 1,
                                                                   materialParameters, blackHole);
        }
        
        computeSLICTimeStep(currentCells, cellSpacing, timeStep, bias, slopeLimiter, materialParameters, blackHole);
        
        for (int i = 0; i < subcyclingIterations; i++)
        {
            BlackHoleEulerForcingSolver::computeRungeKuttaTimeStep(currentCells, cellSpacing, 0.5 * (timeStep / subcyclingIterations), 0.0, 1,
                                                                   materialParameters, blackHole);
        }
        
        currentTime += timeStep;
        currentIteration += 1;
        
        EulerSolvers::outputStatus(currentIteration, currentTime, timeStep);
    }
    
    return currentCells;
}

vector<vector<BlackHoleEulerStateVector> > BlackHoleEulerSecondOrderSolver::solve2D(vector<vector<BlackHoleEulerStateVector> > initialCells, double cellSpacing,
                                                                                    double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                                                    int subcyclingIterations, EulerMaterialParameters materialParameters,
                                                                                    BlackHoleSpacetime blackHole)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<BlackHoleEulerStateVector> > currentCells = initialCells;
    
    while (currentTime < finalTime)
    {
        double timeStep = BlackHoleEulerSolvers::computeStableTimeStep2D(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration,
                                                                         materialParameters, blackHole);
        
        for (int i = 0; i < subcyclingIterations; i++)
        {
            BlackHoleEulerForcingSolver::computeRungeKuttaTimeStep2D(currentCells, cellSpacing, 0.5 * (timeStep / subcyclingIterations), 0.0, 1,
                                                                     materialParameters, blackHole);
        }
        
        computeXSLICTimeStep2D(currentCells, cellSpacing, timeStep * 0.5, bias, slopeLimiter, materialParameters, blackHole);
        computeYSLICTimeStep2D(currentCells, cellSpacing, timeStep, bias, slopeLimiter, materialParameters, blackHole);
        computeXSLICTimeStep2D(currentCells, cellSpacing, timeStep * 0.5, bias, slopeLimiter, materialParameters, blackHole);
        
        for (int i = 0; i < subcyclingIterations; i++)
        {
            BlackHoleEulerForcingSolver::computeRungeKuttaTimeStep2D(currentCells, cellSpacing, 0.5 * (timeStep / subcyclingIterations), 0.0, 1,
                                                                     materialParameters, blackHole);
        }
        
        currentTime += timeStep;
        currentIteration += 1;
        
        EulerSolvers::outputStatus(currentIteration, currentTime, timeStep);
    }
    
    return currentCells;
}
