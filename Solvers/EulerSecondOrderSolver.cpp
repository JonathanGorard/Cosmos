#include "EulerSecondOrderSolver.hpp"

vector<double> EulerSecondOrderSolver::computeXSLICFlux(EulerStateVector leftLeftStateVector, EulerStateVector leftStateVector, EulerStateVector rightStateVector,
                                                        EulerStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                        EulerMaterialParameters materialParameters)
{
    EulerStateVector evolvedRightStateVector = EulerSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing,
                                                                                        timeStep, bias, slopeLimiter, 1, materialParameters);
    EulerStateVector evolvedLeftStateVector = EulerSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing,
                                                                                       timeStep, bias, slopeLimiter, 0, materialParameters);
    
    return EulerFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> EulerSecondOrderSolver::computeYSLICFlux(EulerStateVector topTopStateVector, EulerStateVector topStateVector, EulerStateVector bottomStateVector,
                                                        EulerStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                        EulerMaterialParameters materialParameters)
{
    EulerStateVector evolvedBottomStateVector = EulerSolvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing,
                                                                                         timeStep, bias, slopeLimiter, 1, materialParameters);
    EulerStateVector evolvedTopStateVector = EulerSolvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing,
                                                                                      timeStep, bias, slopeLimiter, 0, materialParameters);
    
    return EulerFirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, materialParameters);
}

void EulerSecondOrderSolver::computeSLICTimeStep(vector<EulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                 EulerMaterialParameters materialParameters)
{
    long cellCount = currentCells.size();
    vector<EulerStateVector> currentCellsWithBoundary = EulerSolvers::insertBoundaryCells(currentCells, 2);
    
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

 void EulerSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                                     int slopeLimiter, EulerMaterialParameters materialParameters)
{
     long rowCount = currentCells.size();
     long columnCount = currentCells[0].size();
     vector<vector<EulerStateVector> > currentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(currentCells, 2);
     
#pragma omp parallel for
     for (int i= 0; i < rowCount; i++)
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

void EulerSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                                    int slopeLimiter, EulerMaterialParameters materialParameters)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<EulerStateVector> > currentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(currentCells, 2);
    
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

vector<EulerStateVector> EulerSecondOrderSolver::solve(vector<EulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                       double bias, int slopeLimiter, EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<EulerStateVector> currentCells = initialCells;
    
    while (currentTime < finalTime)
    {
        double timeStep = EulerSolvers::computeStableTimeStep(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);
        
        computeSLICTimeStep(currentCells, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
        
        currentTime += timeStep;
        currentIteration += 1;
        
        EulerSolvers::outputStatus(currentIteration, currentTime, timeStep);
    }
    
    return currentCells;
}

vector<vector<EulerStateVector> > EulerSecondOrderSolver::solve2D(vector<vector<EulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient,
                                                                  double finalTime, double bias, int slopeLimiter, EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<EulerStateVector> > currentCells = initialCells;
    
    while (currentTime < finalTime)
    {
        double timeStep = EulerSolvers::computeStableTimeStep2D(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration,
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
