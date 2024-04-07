#include "BlackHoleEulerFirstOrderSolver.hpp"

vector<double> BlackHoleEulerFirstOrderSolver::computeXLaxFriedrichsFlux(BlackHoleEulerStateVector leftStateVector, BlackHoleEulerStateVector rightStateVector,
                                                                         double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing,
                                                                         double timeStep, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(xCoordinate - (cellSpacing * 0.5), yCoordinate, zCoordinate,
                                                                                                materialParameters, blackHole);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(xCoordinate + (cellSpacing * 0.5), yCoordinate, zCoordinate,
                                                                                                  materialParameters, blackHole);
    
    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(xCoordinate - (cellSpacing * 0.5), yCoordinate, zCoordinate, materialParameters, blackHole);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(xCoordinate + (cellSpacing * 0.5), yCoordinate, zCoordinate, materialParameters, blackHole);
    
    return EulerFirstOrderSolver::computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector,
                                                           cellSpacing, timeStep);
}

vector<double> BlackHoleEulerFirstOrderSolver::computeXRichtmyerFlux(BlackHoleEulerStateVector leftStateVector, BlackHoleEulerStateVector rightStateVector,
                                                                     double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing,
                                                                     double timeStep, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(xCoordinate - (cellSpacing * 0.5), yCoordinate, zCoordinate,
                                                                                                materialParameters, blackHole);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(xCoordinate + (cellSpacing * 0.5), yCoordinate, zCoordinate,
                                                                                                  materialParameters, blackHole);
    
    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(xCoordinate - (cellSpacing * 0.5), yCoordinate, zCoordinate, materialParameters, blackHole);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(xCoordinate + (cellSpacing * 0.5), yCoordinate, zCoordinate, materialParameters, blackHole);
    
    vector<double> intermediateStateVector = EulerFirstOrderSolver::computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector,
                                                                                         rightFluxVector, cellSpacing, timeStep);
    return BlackHoleEulerStateVector::computeXFluxVector(intermediateStateVector, xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
}

vector<double> BlackHoleEulerFirstOrderSolver::computeXFORCEFlux(BlackHoleEulerStateVector leftStateVector, BlackHoleEulerStateVector rightStateVector,
                                                                 double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing, double timeStep,
                                                                 EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, xCoordinate, yCoordinate, zCoordinate, cellSpacing, timeStep,
                                                                 materialParameters, blackHole);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, xCoordinate, yCoordinate, zCoordinate, cellSpacing, timeStep,
                                                         materialParameters, blackHole);
    
    return EulerFirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> BlackHoleEulerFirstOrderSolver::computeYLaxFriedrichsFlux(BlackHoleEulerStateVector topStateVector, BlackHoleEulerStateVector bottomStateVector,
                                                                         double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing,
                                                                         double timeStep, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(xCoordinate, yCoordinate - (cellSpacing * 0.5), zCoordinate,
                                                                                              materialParameters, blackHole);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(xCoordinate, yCoordinate + (cellSpacing * 0.5), zCoordinate,
                                                                                                    materialParameters, blackHole);
    
    vector<double> topFluxVector = topStateVector.computeYFluxVector(xCoordinate, yCoordinate - (cellSpacing * 0.5), zCoordinate, materialParameters, blackHole);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(xCoordinate, yCoordinate + (cellSpacing * 0.5), zCoordinate, materialParameters, blackHole);
    
    return EulerFirstOrderSolver::computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector,
                                                           cellSpacing, timeStep);
}

vector<double> BlackHoleEulerFirstOrderSolver::computeYRichtmyerFlux(BlackHoleEulerStateVector topStateVector, BlackHoleEulerStateVector bottomStateVector,
                                                                     double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing,
                                                                     double timeStep, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(xCoordinate, yCoordinate - (cellSpacing * 0.5), zCoordinate,
                                                                                              materialParameters, blackHole);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(xCoordinate, yCoordinate + (cellSpacing * 0.5), zCoordinate,
                                                                                                    materialParameters, blackHole);
    
    vector<double> topFluxVector = topStateVector.computeYFluxVector(xCoordinate, yCoordinate - (cellSpacing * 0.5), zCoordinate, materialParameters, blackHole);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(xCoordinate, yCoordinate + (cellSpacing * 0.5), zCoordinate, materialParameters, blackHole);
    
    vector<double> intermediateStateVector = EulerFirstOrderSolver::computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector,
                                                                                         bottomFluxVector, cellSpacing, timeStep);
    return BlackHoleEulerStateVector::computeYFluxVector(intermediateStateVector, xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
}

vector<double> BlackHoleEulerFirstOrderSolver::computeYFORCEFlux(BlackHoleEulerStateVector topStateVector, BlackHoleEulerStateVector bottomStateVector,
                                                                 double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing, double timeStep,
                                                                 EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, xCoordinate, yCoordinate, zCoordinate, cellSpacing, timeStep,
                                                                 materialParameters, blackHole);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, xCoordinate, yCoordinate, zCoordinate, cellSpacing, timeStep,
                                                         materialParameters, blackHole);
    
    return EulerFirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

void BlackHoleEulerFirstOrderSolver::computeFORCETimeStep(vector<BlackHoleEulerStateVector> & currentCells, double cellSpacing, double timeStep,
                                                          EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long cellCount = currentCells.size();
    vector<BlackHoleEulerStateVector> currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        if (!blackHole.inExcisionRegion(i * cellSpacing, 0.0, 0.0) && !blackHole.inExcisionRegion((i - 1) * cellSpacing, 0.0, 0.0) &&
            !blackHole.inExcisionRegion((i + 1) * cellSpacing, 0.0, 0.0))
        {
            vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(i * cellSpacing, 0.0, 0.0, materialParameters, blackHole);
            
            vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], (i - 0.5) * cellSpacing, 0.0, 0.0,
                                                              cellSpacing, timeStep, materialParameters, blackHole);
            vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], (i + 0.5) * cellSpacing, 0.0, 0.0,
                                                               cellSpacing, timeStep, materialParameters, blackHole);
            
            currentCells[i].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                 cellSpacing, timeStep), i * cellSpacing, 0.0, 0.0,
                                                       materialParameters, blackHole);
        }
    }
}

void BlackHoleEulerFirstOrderSolver::computeXFORCETimeStep2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                             EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<BlackHoleEulerStateVector> > currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (!blackHole.inExcisionRegion(j * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j - 1) * cellSpacing, i * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion((j + 1) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i - 1) * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion(j * cellSpacing, (i + 1) * cellSpacing, 0.0))
            {
                vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(j * cellSpacing, i * cellSpacing, 0.0, materialParameters,
                                                                                                           blackHole);
                
                vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j], currentCellsWithBoundary[i + 1][j + 1], (j - 0.5) * cellSpacing,
                                                                  i * cellSpacing, 0.0, cellSpacing, timeStep, materialParameters, blackHole);
                vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 1][j + 2], (j + 0.5) * cellSpacing,
                                                                   i * cellSpacing, 0.0, cellSpacing, timeStep, materialParameters, blackHole);
                
                currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                        cellSpacing, timeStep), j * cellSpacing, i * cellSpacing,
                                                              0.0, materialParameters, blackHole);
            }
        }
    }
}

void BlackHoleEulerFirstOrderSolver::computeYFORCETimeStep2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                             EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<BlackHoleEulerStateVector> > currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (!blackHole.inExcisionRegion(j * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j - 1) * cellSpacing, i * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion((j + 1) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i - 1) * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion(j * cellSpacing, (i + 1) * cellSpacing, 0.0))
            {
                vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(j * cellSpacing, i * cellSpacing, 0.0, materialParameters,
                                                                                                           blackHole);
                
                vector<double> topFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i][j + 1], currentCellsWithBoundary[i + 1][j + 1], j * cellSpacing,
                                                                 (i - 0.5) * cellSpacing, 0.0, cellSpacing, timeStep, materialParameters, blackHole);
                vector<double> bottomFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 2][j + 1], j * cellSpacing,
                                                                    (i + 0.5) * cellSpacing, 0.0, cellSpacing, timeStep, materialParameters, blackHole);
                
                currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector,
                                                                                                        cellSpacing, timeStep), j * cellSpacing, i * cellSpacing,
                                                              0.0, materialParameters, blackHole);
            }
        }
    }
}

vector<BlackHoleEulerStateVector> BlackHoleEulerFirstOrderSolver::solve(vector<BlackHoleEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient,
                                                                        double finalTime, int subcyclingIterations, EulerMaterialParameters materialParameters,
                                                                        BlackHoleSpacetime blackHole)
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
            BlackHoleEulerForcingSolver::computeRungeKuttaTimeStep(currentCells, cellSpacing, 0.5 * (timeStep / subcyclingIterations), 0.0, 1, materialParameters, blackHole);
        }
        
        computeFORCETimeStep(currentCells, cellSpacing, timeStep, materialParameters, blackHole);
        
        for (int i = 0; i < subcyclingIterations; i++)
        {
            BlackHoleEulerForcingSolver::computeRungeKuttaTimeStep(currentCells, cellSpacing, 0.5 * (timeStep / subcyclingIterations), 0.0, 1, materialParameters, blackHole);
        }
        
        currentTime += timeStep;
        currentIteration += 1;
        
        EulerSolvers::outputStatus(currentIteration, currentTime, timeStep);
    }
    
    return currentCells;
}

vector<vector<BlackHoleEulerStateVector> > BlackHoleEulerFirstOrderSolver::solve2D(vector<vector<BlackHoleEulerStateVector> > initialCells, double cellSpacing,
                                                                                   double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                                                   EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
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
        
        computeXFORCETimeStep2D(currentCells, cellSpacing, timeStep * 0.5, materialParameters, blackHole);
        computeYFORCETimeStep2D(currentCells, cellSpacing, timeStep, materialParameters, blackHole);
        computeXFORCETimeStep2D(currentCells, cellSpacing, timeStep * 0.5, materialParameters, blackHole);
        
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
