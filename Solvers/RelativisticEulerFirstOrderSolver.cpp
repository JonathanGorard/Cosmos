#include "RelativisticEulerFirstOrderSolver.hpp"

vector<double> RelativisticEulerFirstOrderSolver::computeXLaxFriedrichsFlux(RelativisticEulerStateVector leftStateVector, RelativisticEulerStateVector rightStateVector,
                                                                            double cellSpacing, double timeStep, EulerMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);
    
    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);
    
    return EulerFirstOrderSolver::computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector,
                                                           cellSpacing, timeStep);
}

vector<double> RelativisticEulerFirstOrderSolver::computeXRichtmyerFlux(RelativisticEulerStateVector leftStateVector, RelativisticEulerStateVector rightStateVector,
                                                                        double cellSpacing, double timeStep, EulerMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);
    
    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);
    
    vector<double> intermediateStateVector = EulerFirstOrderSolver::computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector,
                                                                                         rightFluxVector, cellSpacing, timeStep);
    return RelativisticEulerStateVector::computeXFluxVector(intermediateStateVector, materialParameters);
}

vector<double> RelativisticEulerFirstOrderSolver::computeXFORCEFlux(RelativisticEulerStateVector leftStateVector, RelativisticEulerStateVector rightStateVector,
                                                                    double cellSpacing, double timeStep, EulerMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);
    
    return EulerFirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> RelativisticEulerFirstOrderSolver::computeYLaxFriedrichsFlux(RelativisticEulerStateVector topStateVector, RelativisticEulerStateVector bottomStateVector,
                                                                            double cellSpacing, double timeStep, EulerMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);
    
    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);
    
    return EulerFirstOrderSolver::computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector,
                                                           cellSpacing, timeStep);
}

vector<double> RelativisticEulerFirstOrderSolver::computeYRichtmyerFlux(RelativisticEulerStateVector topStateVector, RelativisticEulerStateVector bottomStateVector,
                                                                        double cellSpacing, double timeStep, EulerMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);
    
    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);
    
    vector<double> intermediateStateVector = EulerFirstOrderSolver::computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector,
                                                                                         bottomFluxVector, cellSpacing, timeStep);
    return RelativisticEulerStateVector::computeYFluxVector(intermediateStateVector, materialParameters);
}

vector<double> RelativisticEulerFirstOrderSolver::computeYFORCEFlux(RelativisticEulerStateVector topStateVector, RelativisticEulerStateVector bottomStateVector,
                                                                    double cellSpacing, double timeStep, EulerMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);
    
    return EulerFirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

void RelativisticEulerFirstOrderSolver::computeFORCETimeStep(vector<RelativisticEulerStateVector> & currentCells, double cellSpacing, double timeStep,
                                                             EulerMaterialParameters materialParameters)
{
    long cellCount = currentCells.size();
    vector<RelativisticEulerStateVector> currentCellsWithBoundary = RelativisticEulerSolvers::insertBoundaryCells(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);
        
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, materialParameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, materialParameters);
        
        currentCells[i].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                             cellSpacing, timeStep), materialParameters);
    }
}

void RelativisticEulerFirstOrderSolver::computeXFORCETimeStep2D(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                                EulerMaterialParameters materialParameters)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<RelativisticEulerStateVector> > currentCellsWithBoundary = RelativisticEulerSolvers::insertBoundaryCells2D(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);
            
            vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j], currentCellsWithBoundary[i + 1][j + 1], cellSpacing,
                                                              timeStep, materialParameters);
            vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 1][j + 2], cellSpacing,
                                                               timeStep, materialParameters);
            
            currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                    cellSpacing, timeStep), materialParameters);
        }
    }
}

void RelativisticEulerFirstOrderSolver::computeYFORCETimeStep2D(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                                EulerMaterialParameters materialParameters)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<RelativisticEulerStateVector> > currentCellsWithBoundary = RelativisticEulerSolvers::insertBoundaryCells2D(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);
            
            vector<double> topFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i][j + 1], currentCellsWithBoundary[i + 1][j + 1], cellSpacing,
                                                             timeStep, materialParameters);
            vector<double> bottomFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 2][j + 1], cellSpacing,
                                                                timeStep, materialParameters);
            
            currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector,
                                                                                                    cellSpacing, timeStep), materialParameters);
        }
    }
}

vector<RelativisticEulerStateVector> RelativisticEulerFirstOrderSolver::solve(vector<RelativisticEulerStateVector> initialCells, double cellSpacing,
                                                                              double CFLCoefficient, double finalTime, EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<RelativisticEulerStateVector> currentCells = initialCells;
    
    while (currentTime < finalTime)
    {
        double timeStep = RelativisticEulerSolvers::computeStableTimeStep(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration,
                                                                          materialParameters);
        
        computeFORCETimeStep(currentCells, cellSpacing, timeStep, materialParameters);
        
        currentTime += timeStep;
        currentIteration += 1;
        
        EulerSolvers::outputStatus(currentIteration, currentTime, timeStep);
    }
    
    return currentCells;
}

vector<vector<RelativisticEulerStateVector> > RelativisticEulerFirstOrderSolver::solve2D(vector<vector<RelativisticEulerStateVector> > initialCells,
                                                                                         double cellSpacing, double CFLCoefficient, double finalTime,
                                                                                         EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<RelativisticEulerStateVector> > currentCells = initialCells;
    
    while (currentTime < finalTime)
    {
        double timeStep = RelativisticEulerSolvers::computeStableTimeStep2D(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration,
                                                                            materialParameters);
        
        computeXFORCETimeStep2D(currentCells, cellSpacing, timeStep * 0.5, materialParameters);
        computeYFORCETimeStep2D(currentCells, cellSpacing, timeStep, materialParameters);
        computeXFORCETimeStep2D(currentCells, cellSpacing, timeStep * 0.5, materialParameters);
        
        currentTime += timeStep;
        currentIteration += 1;
        
        EulerSolvers::outputStatus(currentIteration, currentTime, timeStep);
    }
    
    return currentCells;
}
