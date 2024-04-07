#include "EulerFirstOrderSolver.hpp"

vector<double> EulerFirstOrderSolver::computeLaxFriedrichsFlux(vector<double> leftConservedVariableVector, vector<double> rightConservedVariableVector,
                                                               vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing, double timeStep)
{
    return VectorAlgebra::multiplyVector(0.5,
                                         VectorAlgebra::addVectors(VectorAlgebra::addVectors(leftFluxVector, rightFluxVector),
                                                                   VectorAlgebra::multiplyVector((cellSpacing / timeStep),
                                                                                                 VectorAlgebra::subtractVectors(leftConservedVariableVector,
                                                                                                                                rightConservedVariableVector))));
}

vector<double> EulerFirstOrderSolver::computeRichtmyerFlux(vector<double> leftConservedVariableVector, vector<double> rightConservedVariableVector,
                                                           vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing, double timeStep)
{
    return VectorAlgebra::multiplyVector(0.5, VectorAlgebra::addVectors(VectorAlgebra::addVectors(leftConservedVariableVector, rightConservedVariableVector),
                                                                        VectorAlgebra::multiplyVector((timeStep / cellSpacing),
                                                                                                      VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector))));
}

vector<double> EulerFirstOrderSolver::computeFORCEFlux(vector<double> laxFriedrichsFlux, vector<double> richtmyerFlux)
{
    return VectorAlgebra::multiplyVector(0.5, VectorAlgebra::addVectors(laxFriedrichsFlux, richtmyerFlux));
}

vector<double> EulerFirstOrderSolver::computeXLaxFriedrichsFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                                EulerMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);
    
    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);
    
    return computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> EulerFirstOrderSolver::computeXRichtmyerFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                            EulerMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);
    
    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);
    
    vector<double> intermediateStateVector = computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector,
                                                                  cellSpacing, timeStep);
    return EulerStateVector::computeXFluxVector(intermediateStateVector, materialParameters);
}

vector<double> EulerFirstOrderSolver::computeXFORCEFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                        EulerMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);
    
    return computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> EulerFirstOrderSolver::computeYLaxFriedrichsFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                                EulerMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);
    
    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);
    
    return computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
}

vector<double> EulerFirstOrderSolver::computeYRichtmyerFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                            EulerMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);
    
    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);
    
    vector<double> intermediateStateVector = computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector,
                                                                  cellSpacing, timeStep);
    return EulerStateVector::computeYFluxVector(intermediateStateVector, materialParameters);
}

vector<double> EulerFirstOrderSolver::computeYFORCEFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                        EulerMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);
    
    return computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> EulerFirstOrderSolver::computeFORCEUpdate(vector<double> conservedVariableVector, vector<double> leftFluxVector, vector<double> rightFluxVector,
                                                         double cellSpacing, double timeStep)
{
    return VectorAlgebra::addVectors(conservedVariableVector, VectorAlgebra::multiplyVector((timeStep / cellSpacing),
                                                                                            VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector)));
}

void EulerFirstOrderSolver::computeFORCETimeStep(vector<EulerStateVector> & currentCells, double cellSpacing, double timeStep, EulerMaterialParameters materialParameters)
{
    long cellCount = currentCells.size();
    vector<EulerStateVector> currentCellsWithBoundary = EulerSolvers::insertBoundaryCells(currentCells, 1);

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);
        
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, materialParameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, materialParameters);
        
        currentCells[i].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                   materialParameters);
    }
}

void EulerFirstOrderSolver::computeXFORCETimeStep2D(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters materialParameters)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<EulerStateVector> > currentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(currentCells, 1);
    
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
            
            currentCells[i][j].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                          materialParameters);
        }
    }
}

void EulerFirstOrderSolver::computeYFORCETimeStep2D(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters materialParameters)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<EulerStateVector> > currentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(currentCells, 1);
    
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
            
            currentCells[i][j].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep),
                                                          materialParameters);
        }
    }
}

vector<EulerStateVector> EulerFirstOrderSolver::solve(vector<EulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                      EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<EulerStateVector> currentCells = initialCells;
    
    while (currentTime < finalTime)
    {
        double timeStep = EulerSolvers::computeStableTimeStep(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);
        
        computeFORCETimeStep(currentCells, cellSpacing, timeStep, materialParameters);
        
        currentTime += timeStep;
        currentIteration += 1;
        
        EulerSolvers::outputStatus(currentIteration, currentTime, timeStep);
    }
    
    return currentCells;
}

vector<vector<EulerStateVector> > EulerFirstOrderSolver::solve2D(vector<vector<EulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient,
                                                                 double finalTime, EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<EulerStateVector> > currentCells = initialCells;
    
    while (currentTime < finalTime)
    {
        double timeStep = EulerSolvers::computeStableTimeStep2D(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration,
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
