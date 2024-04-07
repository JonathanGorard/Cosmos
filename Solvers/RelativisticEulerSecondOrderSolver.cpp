#include "RelativisticEulerSecondOrderSolver.hpp"

vector<double> RelativisticEulerSecondOrderSolver::computeXSLICFlux(RelativisticEulerStateVector leftLeftStateVector, RelativisticEulerStateVector leftStateVector,
                                                                    RelativisticEulerStateVector rightStateVector, RelativisticEulerStateVector rightRightStateVector,
                                                                    double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                    EulerMaterialParameters materialParameters)
{
    RelativisticEulerStateVector evolvedRightStateVector = RelativisticEulerSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector,
                                                                                                                rightStateVector, cellSpacing, timeStep, bias,
                                                                                                                slopeLimiter, 1, materialParameters);
    RelativisticEulerStateVector evolvedLeftStateVector = RelativisticEulerSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector,
                                                                                                               rightRightStateVector, cellSpacing, timeStep, bias,
                                                                                                               slopeLimiter, 0, materialParameters);
    
    return RelativisticEulerFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> RelativisticEulerSecondOrderSolver::computeYSLICFlux(RelativisticEulerStateVector topTopStateVector, RelativisticEulerStateVector topStateVector,
                                                                    RelativisticEulerStateVector bottomStateVector, RelativisticEulerStateVector bottomBottomStateVector,
                                                                    double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                    EulerMaterialParameters materialParameters)
{
    RelativisticEulerStateVector evolvedBottomStateVector = RelativisticEulerSolvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector,
                                                                                                                 cellSpacing, timeStep, bias, slopeLimiter, 1,
                                                                                                                 materialParameters);
    RelativisticEulerStateVector evolvedTopStateVector = RelativisticEulerSolvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector,
                                                                                                              cellSpacing, timeStep, bias, slopeLimiter, 0,
                                                                                                              materialParameters);
    
    return RelativisticEulerFirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, materialParameters);
}

void RelativisticEulerSecondOrderSolver::computeSLICTimeStep(vector<RelativisticEulerStateVector> & currentCells, double cellSpacing, double timeStep,
                                                             double bias, int slopeLimiter, EulerMaterialParameters materialParameters)
{
    long cellCount = currentCells.size();
    vector<RelativisticEulerStateVector> currentCellsWithBoundary = RelativisticEulerSolvers::insertBoundaryCells(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);
        
        vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2],
                                                         currentCellsWithBoundary[i + 3], cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
        vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                                                          currentCellsWithBoundary[i + 4], cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
        
        currentCells[i].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing,
                                                                                             timeStep), materialParameters);
    }
}

void RelativisticEulerSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                                double bias, int slopeLimiter, EulerMaterialParameters materialParameters)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<RelativisticEulerStateVector> > currentCellsWithBoundary = RelativisticEulerSolvers::insertBoundaryCells2D(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);
            
            vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 2][j], currentCellsWithBoundary[i + 2][j + 1],
                                                             currentCellsWithBoundary[i + 2][j + 2], currentCellsWithBoundary[i + 2][j + 3], cellSpacing, timeStep,
                                                             bias, slopeLimiter, materialParameters);
            vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 2][j + 1], currentCellsWithBoundary[i + 2][j + 2],
                                                              currentCellsWithBoundary[i + 2][j + 3], currentCellsWithBoundary[i + 2][j + 4], cellSpacing, timeStep,
                                                              bias, slopeLimiter, materialParameters);
            
            currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                    cellSpacing, timeStep), materialParameters);
        }
    }
}

void RelativisticEulerSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                                double bias, int slopeLimiter, EulerMaterialParameters materialParameters)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<RelativisticEulerStateVector> > currentCellsWithBoundary = RelativisticEulerSolvers::insertBoundaryCells2D(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);
            
            vector<double> topFluxVector = computeYSLICFlux(currentCellsWithBoundary[i][j + 2], currentCellsWithBoundary[i + 1][j + 2],
                                                            currentCellsWithBoundary[i + 2][j + 2], currentCellsWithBoundary[i + 3][j + 2], cellSpacing, timeStep,
                                                            bias, slopeLimiter, materialParameters);
            vector<double> bottomFluxVector = computeYSLICFlux(currentCellsWithBoundary[i + 1][j + 2], currentCellsWithBoundary[i + 2][j + 2],
                                                               currentCellsWithBoundary[i + 3][j + 2], currentCellsWithBoundary[i + 4][j + 2], cellSpacing, timeStep,
                                                               bias, slopeLimiter, materialParameters);
            
            currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector,
                                                                                                    cellSpacing, timeStep), materialParameters);
        }
    }
}

vector<RelativisticEulerStateVector> RelativisticEulerSecondOrderSolver::solve(vector<RelativisticEulerStateVector> initialCells, double cellSpacing,
                                                                               double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                                               EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<RelativisticEulerStateVector> currentCells = initialCells;
    
    while (currentTime < finalTime)
    {
        double timeStep = RelativisticEulerSolvers::computeStableTimeStep(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration,
                                                                          materialParameters);
        
        computeSLICTimeStep(currentCells, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
        
        currentTime += timeStep;
        currentIteration += 1;
        
        EulerSolvers::outputStatus(currentIteration, currentTime, timeStep);
    }
    
    return currentCells;
}

vector<vector<RelativisticEulerStateVector> > RelativisticEulerSecondOrderSolver::solve2D(vector<vector<RelativisticEulerStateVector> > initialCells, double cellSpacing,
                                                                                          double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                                                          EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<RelativisticEulerStateVector> > currentCells = initialCells;
    
    while (currentTime < finalTime)
    {
        double timeStep = RelativisticEulerSolvers::computeStableTimeStep2D(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration,
                                                                            materialParameters);
        
        computeXSLICTimeStep2D(currentCells, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, materialParameters);
        computeYSLICTimeStep2D(currentCells, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
        computeXSLICTimeStep2D(currentCells, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, materialParameters);
        
        currentTime += timeStep;
        currentIteration += 1;
        
        EulerSolvers::outputStatus(currentIteration, currentTime, timeStep);
    }
    
    return currentCells;
}
