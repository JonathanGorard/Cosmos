#include "BlackHoleEulerAMR.hpp"

void BlackHoleEulerAMR::computeFORCETimeStepAMR(vector<BlackHoleEulerStateVector> & currentCells, double cellSpacing, double timeStep, vector<bool> AMRStructure,
                                                EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long cellCount = currentCells.size();
    vector<BlackHoleEulerStateVector> currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        if (AMRStructure[i])
        {
            if (!blackHole.inExcisionRegion(i * cellSpacing, 0.0, 0.0) && !blackHole.inExcisionRegion((i - 1) * cellSpacing, 0.0, 0.0) &&
                !blackHole.inExcisionRegion((i + 1) * cellSpacing, 0.0, 0.0))
            {
                vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(i * cellSpacing, 0.0, 0.0, materialParameters, blackHole);
                
                vector<double> leftFluxVector = BlackHoleEulerFirstOrderSolver::computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1],
                                                                                                  (i - 0.5) * cellSpacing, 0.0, 0.0, cellSpacing, timeStep,
                                                                                                  materialParameters, blackHole);
                vector<double> rightFluxVector = BlackHoleEulerFirstOrderSolver::computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2],
                                                                                                   (i + 0.5) * cellSpacing, 0.0, 0.0, cellSpacing, timeStep,
                                                                                                   materialParameters, blackHole);
                
                currentCells[i].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                     cellSpacing, timeStep), i * cellSpacing, 0.0, 0.0,
                                                           materialParameters, blackHole);
            }
        }
    }
}

void BlackHoleEulerAMR::computeXFORCETimeStep2DAMR(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                   vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<BlackHoleEulerStateVector> > currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (AMRStructure[i][j])
            {
                if (!blackHole.inExcisionRegion(j * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j - 1) * cellSpacing, i * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion((j + 1) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i - 1) * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion(j * cellSpacing, (i + 1) * cellSpacing, 0.0))
                {
                    vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(j * cellSpacing, i * cellSpacing, 0.0,
                                                                                                               materialParameters, blackHole);
                    
                    vector<double> leftFluxVector = BlackHoleEulerFirstOrderSolver::computeXFORCEFlux(currentCellsWithBoundary[i + 1][j],
                                                                                                      currentCellsWithBoundary[i + 1][j + 1], (j - 0.5) * cellSpacing,
                                                                                                      i * cellSpacing, 0.0, cellSpacing, timeStep, materialParameters,
                                                                                                      blackHole);
                    vector<double> rightFluxVector = BlackHoleEulerFirstOrderSolver::computeXFORCEFlux(currentCellsWithBoundary[i + 1][j + 1],
                                                                                                       currentCellsWithBoundary[i + 1][j + 2], (j + 0.5) * cellSpacing,
                                                                                                       i * cellSpacing, 0.0, cellSpacing, timeStep, materialParameters,
                                                                                                       blackHole);
                    
                    currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                            cellSpacing, timeStep), j * cellSpacing, i * cellSpacing, 0.0,
                                                                  materialParameters, blackHole);
                }
            }
        }
    }
}

void BlackHoleEulerAMR::computeYFORCETimeStep2DAMR(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                   vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<BlackHoleEulerStateVector> > currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (AMRStructure[i][j])
            {
                if (!blackHole.inExcisionRegion(j * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j - 1) * cellSpacing, i * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion((j + 1) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i - 1) * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion(j * cellSpacing, (i + 1) * cellSpacing, 0.0))
                {
                    vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(j * cellSpacing, i * cellSpacing, 0.0,
                                                                                                               materialParameters, blackHole);
                    
                    vector<double> topFluxVector = BlackHoleEulerFirstOrderSolver::computeYFORCEFlux(currentCellsWithBoundary[i][j + 1],
                                                                                                     currentCellsWithBoundary[i + 1][j + 1], j * cellSpacing,
                                                                                                     (i - 0.5) * cellSpacing, 0.0, cellSpacing, timeStep,
                                                                                                     materialParameters, blackHole);
                    vector<double> bottomFluxVector = BlackHoleEulerFirstOrderSolver::computeYFORCEFlux(currentCellsWithBoundary[i + 1][j + 1],
                                                                                                        currentCellsWithBoundary[i + 2][j + 1], j * cellSpacing,
                                                                                                        (i + 0.5) * cellSpacing, 0.0, cellSpacing, timeStep,
                                                                                                        materialParameters, blackHole);
                    
                    currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector,
                                                                                                            cellSpacing, timeStep), j * cellSpacing, i * cellSpacing, 0.0,
                                                                  materialParameters, blackHole);
                }
            }
        }
    }
}

void BlackHoleEulerAMR::computeSLICTimeStepAMR(vector<BlackHoleEulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias,
                                               int slopeLimiter, vector<bool> AMRStructure, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long cellCount = currentCells.size();
    vector<BlackHoleEulerStateVector> currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        if (AMRStructure[i])
        {
            if (!blackHole.inExcisionRegion(i * cellSpacing, 0.0, 0.0) && !blackHole.inExcisionRegion((i - 1) * cellSpacing, 0.0, 0.0) &&
                !blackHole.inExcisionRegion((i - 2) * cellSpacing, 0.0, 0.0) && !blackHole.inExcisionRegion((i + 1) * cellSpacing, 0.0, 0.0) &&
                !blackHole.inExcisionRegion((i + 2) * cellSpacing, 0.0, 0.0))
            {
                vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(i * cellSpacing, 0.0, 0.0, materialParameters, blackHole);
                
                vector<double> leftFluxVector = BlackHoleEulerSecondOrderSolver::computeXSLICFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1],
                                                                                                  currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                                                                                                  cellSpacing, timeStep, bias, slopeLimiter, (i - 0.5) * cellSpacing,
                                                                                                  0.0, 0.0, materialParameters, blackHole);
                vector<double> rightFluxVector = BlackHoleEulerSecondOrderSolver::computeXSLICFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2],
                                                                                                   currentCellsWithBoundary[i + 3], currentCellsWithBoundary[i + 4],
                                                                                                   cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                   (i + 0.5) * cellSpacing, 0.0, 0.0, materialParameters, blackHole);
                
                currentCells[i].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                     cellSpacing, timeStep), i * cellSpacing, 0.0, 0.0,
                                                           materialParameters, blackHole);
            }
        }
    }
}

void BlackHoleEulerAMR::computeXSLICTimeStep2DAMR(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                                  int slopeLimiter, vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters,
                                                  BlackHoleSpacetime blackHole)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<BlackHoleEulerStateVector> > currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (AMRStructure[i][j])
            {
                if (!blackHole.inExcisionRegion(j * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j - 1) * cellSpacing, i * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion((j - 2) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j + 1) * cellSpacing, i * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion((j + 2) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i - 1) * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion(j * cellSpacing, (i - 2) * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i + 1) * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion(j * cellSpacing, (i + 2) * cellSpacing, 0.0))
                {
                    vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(j * cellSpacing, i * cellSpacing, 0.0,
                                                                                                               materialParameters, blackHole);
                    
                    vector<double> leftFluxVector = BlackHoleEulerSecondOrderSolver::computeXSLICFlux(currentCellsWithBoundary[i + 2][j],
                                                                                                      currentCellsWithBoundary[i + 2][j + 1],
                                                                                                      currentCellsWithBoundary[i + 2][j + 2],
                                                                                                      currentCellsWithBoundary[i + 2][j + 3], cellSpacing, timeStep,
                                                                                                      bias, slopeLimiter, (j - 0.5) * cellSpacing, i * cellSpacing, 0.0,
                                                                                                      materialParameters, blackHole);
                    vector<double> rightFluxVector = BlackHoleEulerSecondOrderSolver::computeXSLICFlux(currentCellsWithBoundary[i + 2][j + 1],
                                                                                                       currentCellsWithBoundary[i + 2][j + 2],
                                                                                                       currentCellsWithBoundary[i + 2][j + 3],
                                                                                                       currentCellsWithBoundary[i + 2][j + 4], cellSpacing, timeStep,
                                                                                                       bias, slopeLimiter, (j + 0.5) * cellSpacing, i * cellSpacing, 0.0,
                                                                                                       materialParameters, blackHole);
                    
                    currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                            cellSpacing, timeStep), j * cellSpacing, i * cellSpacing,
                                                                  0.0, materialParameters, blackHole);
                }
            }
        }
    }
}

void BlackHoleEulerAMR::computeYSLICTimeStep2DAMR(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                                  int slopeLimiter, vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters,
                                                  BlackHoleSpacetime blackHole)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<BlackHoleEulerStateVector> >  currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (AMRStructure[i][j])
            {
                if (!blackHole.inExcisionRegion(j * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j - 1) * cellSpacing, i * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion((j - 2) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j + 1) * cellSpacing, i * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion((j + 2) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i - 1) * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion(j * cellSpacing, (i - 2) * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i + 1) * cellSpacing, 0.0) &&
                    !blackHole.inExcisionRegion(j * cellSpacing, (i + 2) * cellSpacing, 0.0))
                {
                    vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(j * cellSpacing, i * cellSpacing, 0.0,
                                                                                                               materialParameters, blackHole);
                    
                    vector<double> topFluxVector = BlackHoleEulerSecondOrderSolver::computeYSLICFlux(currentCellsWithBoundary[i][j + 2],
                                                                                                     currentCellsWithBoundary[i + 1][j + 2],
                                                                                                     currentCellsWithBoundary[i + 2][j + 2],
                                                                                                     currentCellsWithBoundary[i + 3][j + 2], cellSpacing, timeStep,
                                                                                                     bias, slopeLimiter, j * cellSpacing, (i - 0.5) * cellSpacing, 0.0,
                                                                                                     materialParameters, blackHole);
                    vector<double> bottomFluxVector = BlackHoleEulerSecondOrderSolver::computeYSLICFlux(currentCellsWithBoundary[i + 1][j + 2],
                                                                                                        currentCellsWithBoundary[i + 2][j + 2],
                                                                                                        currentCellsWithBoundary[i + 3][j + 2],
                                                                                                        currentCellsWithBoundary[i + 4][j + 2], cellSpacing, timeStep,
                                                                                                        bias, slopeLimiter, j * cellSpacing, (i + 0.5) * cellSpacing, 0.0,
                                                                                                        materialParameters, blackHole);
                    
                    currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector,
                                                                                                            cellSpacing, timeStep), j * cellSpacing, i * cellSpacing,
                                                                  0.0, materialParameters, blackHole);
                }
            }
        }
    }
}

tuple<vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<bool>, vector<bool>>
BlackHoleEulerAMR::solveLevel1AMR(vector<BlackHoleEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                  double AMRTolerance, int order, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long cellCount = initialCells.size();
    
    vector<BlackHoleEulerStateVector> coarseCurrentCells(cellCount);
    vector<BlackHoleEulerStateVector> fineCurrentCells(cellCount * 2);
    vector<bool> coarseAMRStructure(cellCount);
    vector<bool> fineAMRStructure(cellCount * 2);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        coarseAMRStructure[i] = true;
        
        fineAMRStructure[(i * 2)] = false;
        fineAMRStructure[(i * 2) + 1] = false;
    }
    
    double coarseCurrentTime = 0.0;
    int coarseCurrentIteration = 0;
    int fineTotalIteration = 0;
    coarseCurrentCells = initialCells;
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        fineCurrentCells[(i * 2)] = coarseCurrentCells[i];
        fineCurrentCells[(i * 2) + 1] = coarseCurrentCells[i];
    }
    
    while (coarseCurrentTime < finalTime)
    {
        double coarseTimeStep = BlackHoleEulerSolvers::computeStableTimeStep(coarseCurrentCells, cellSpacing, CFLCoefficient, coarseCurrentTime, finalTime,
                                                                             fineTotalIteration, materialParameters, blackHole);
        
        vector<BlackHoleEulerStateVector> coarseCurrentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(coarseCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            vector<double> leftConservedVariableVector = coarseCurrentCellsWithBoundary[i].computeConservedVariableVector((i - 1) * cellSpacing, 0.0, 0.0,
                                                                                                                          materialParameters, blackHole);
            vector<double> conservedVariableVector = coarseCurrentCellsWithBoundary[i + 1].computeConservedVariableVector(i * cellSpacing, 0.0, 0.0,
                                                                                                                          materialParameters, blackHole);
            vector<double> rightConservedVariableVector = coarseCurrentCellsWithBoundary[i + 2].computeConservedVariableVector((i + 1) * cellSpacing, 0.0, 0.0,
                                                                                                                               materialParameters, blackHole);
            
            long conservedVariableCount = conservedVariableVector.size();
            
            for (int j = 0; j < conservedVariableCount; j++)
            {
                if (abs((conservedVariableVector[j] - leftConservedVariableVector[j]) / cellSpacing) >= AMRTolerance ||
                    abs((rightConservedVariableVector[j] - conservedVariableVector[j]) / cellSpacing) >= AMRTolerance)
                {
                    coarseAMRStructure[i] = false;
                    
                    fineAMRStructure[(i * 2)] = true;
                    fineAMRStructure[(i * 2) + 1] = true;
                }
            }
        }
        
        if (order == 1)
        {
            computeFORCETimeStepAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, coarseAMRStructure, materialParameters, blackHole);
        }
        else
        {
            computeSLICTimeStepAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, 0.0, 1, coarseAMRStructure, materialParameters, blackHole);
        }
        
        AMRHelper::outputCoarseStatus(coarseCurrentIteration + 1, coarseCurrentTime + coarseTimeStep, coarseTimeStep);
        
        bool refined = false;
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            if (fineAMRStructure[(i * 2)] || fineAMRStructure[(i * 2) + 1])
            {
                refined = true;
            }
        }
        
        if (refined)
        {
            double fineCurrentTime = coarseCurrentTime;
            int fineCurrentIteration = 0;
            
            while (fineCurrentTime < (coarseCurrentTime + coarseTimeStep - pow(10.0, -4.0)))
            {
                double fineTimeStep = BlackHoleEulerSolvers::computeStableTimeStep(fineCurrentCells, cellSpacing * 0.5, CFLCoefficient, fineCurrentTime,
                                                                                   coarseCurrentTime + coarseTimeStep, fineTotalIteration, materialParameters, blackHole);
                
                if (order == 1)
                {
                    computeFORCETimeStepAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep, fineAMRStructure, materialParameters, blackHole);
                }
                else
                {
                    computeSLICTimeStepAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep, 0.0, 1, fineAMRStructure, materialParameters, blackHole);
                }
                
                fineCurrentTime += fineTimeStep;
                fineCurrentIteration += 1;
                fineTotalIteration += 1;
                
                AMRHelper::outputIntermediateStatus(fineCurrentIteration, fineCurrentTime, fineTimeStep);
            }
        }
        else
        {
            fineTotalIteration += 1;
        }
        
        coarseCurrentTime += coarseTimeStep;
        coarseCurrentIteration += 1;
        
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            if (!coarseAMRStructure[i])
            {
                vector<double> leftFineConservedVariableVector = fineCurrentCells[(i * 2)].computeConservedVariableVector((i * 2) * (cellSpacing * 0.5), 0.0, 0.0,
                                                                                                                          materialParameters, blackHole);
                vector<double> rightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1].computeConservedVariableVector(((i * 2) + 1) * (cellSpacing * 0.5),
                                                                                                                               0.0, 0.0, materialParameters, blackHole);
                
                coarseCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5,
                                                                                               VectorAlgebra::addVectors(leftFineConservedVariableVector,
                                                                                                                         rightFineConservedVariableVector)),
                                                                 i * cellSpacing, 0.0, 0.0, materialParameters, blackHole);
            }
            
            if (!fineAMRStructure[(i * 2)])
            {
                fineCurrentCells[(i * 2)] = coarseCurrentCells[i];
            }
            if (!fineAMRStructure[(i * 2) + 1])
            {
                fineCurrentCells[(i * 2) + 1] = coarseCurrentCells[i];
            }
        }
        
        vector<BlackHoleEulerStateVector> fineCurrentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(fineCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            vector<double> leftConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 1].computeConservedVariableVector(((i * 2) - 1) * (cellSpacing * 0.5), 0.0,
                                                                                                                                  0.0, materialParameters, blackHole);
            vector<double> conservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2].computeConservedVariableVector((i * 2) * (cellSpacing * 0.5), 0.0, 0.0,
                                                                                                                              materialParameters, blackHole);
            vector<double> rightConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 3].computeConservedVariableVector(((i * 2) + 1) * (cellSpacing * 0.5),
                                                                                                                                   0.0, 0.0, materialParameters, blackHole);
            
            long conservedVariableCount = conservedVariableVector.size();
            bool markedCell = false;
            
            for (int j = 0; j < conservedVariableCount; j++)
            {
                if (abs((conservedVariableVector[j] - leftConservedVariableVector[j]) / cellSpacing) >= (AMRTolerance * 0.5) ||
                    abs((rightConservedVariableVector[j] - conservedVariableVector[j]) / cellSpacing) >= (AMRTolerance * 0.5))
                {
                    markedCell = true;
                }
            }
            
            if (!markedCell && fineAMRStructure[(i * 2)])
            {
                fineAMRStructure[(i * 2)] = false;
                fineAMRStructure[(i * 2) + 1] = false;
                
                coarseAMRStructure[i] = true;
                
                vector<double> leftFineConservedVariableVector = fineCurrentCells[(i * 2)].computeConservedVariableVector((i * 2) * (cellSpacing * 0.5), 0.0, 0.0,
                                                                                                                          materialParameters, blackHole);
                vector<double> rightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1].computeConservedVariableVector(((i * 2) + 1) * (cellSpacing * 0.5), 0.0,
                                                                                                                               0.0, materialParameters, blackHole);
                
                coarseCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5, VectorAlgebra::addVectors(leftFineConservedVariableVector,
                                                                                                                              rightFineConservedVariableVector)),
                                                                 i * cellSpacing, 0.0, 0.0, materialParameters, blackHole);
            }
        }
    }
    
    return tuple<vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<bool>, vector<bool>>(coarseCurrentCells, fineCurrentCells,
                                                                                                                   coarseAMRStructure, fineAMRStructure);
}

tuple<vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >, vector<vector<bool> >,vector<vector<bool> >>
BlackHoleEulerAMR::solve2DLevel1AMR(vector<vector<BlackHoleEulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                    double AMRTolerance, int order, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long rowCount = initialCells.size();
    long columnCount = initialCells[0].size();
    
    vector<vector<BlackHoleEulerStateVector> > coarseCurrentCells(rowCount, vector<BlackHoleEulerStateVector>(columnCount));
    vector<vector<BlackHoleEulerStateVector> > fineCurrentCells(rowCount * 2, vector<BlackHoleEulerStateVector>(columnCount * 2));
    vector<vector<bool> > coarseAMRStructure(rowCount, vector<bool>(columnCount));
    vector<vector<bool> > fineAMRStructure(rowCount * 2, vector<bool>(columnCount * 2));
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            coarseAMRStructure[i][j] = true;
            
            fineAMRStructure[(i * 2)][(j * 2)] = false;
            fineAMRStructure[(i * 2)][(j * 2) + 1] = false;
            fineAMRStructure[(i * 2) + 1][(j * 2)] = false;
            fineAMRStructure[(i * 2) + 1][(j * 2) + 1] = false;
        }
    }
    
    double coarseCurrentTime = 0.0;
    int coarseCurrentIteration = 0;
    int fineTotalIteration = 0;
    coarseCurrentCells = initialCells;
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            fineCurrentCells[(i * 2)][(j * 2)] = coarseCurrentCells[i][j];
            fineCurrentCells[(i * 2)][(j * 2) + 1] = coarseCurrentCells[i][j];
            fineCurrentCells[(i * 2) + 1][(j * 2)] = coarseCurrentCells[i][j];
            fineCurrentCells[(i * 2) + 1][(j * 2) + 1] = coarseCurrentCells[i][j];
        }
    }
    
    while (coarseCurrentTime < finalTime)
    {
        double coarseTimeStep = BlackHoleEulerSolvers::computeStableTimeStep2D(coarseCurrentCells, cellSpacing, CFLCoefficient, coarseCurrentTime, finalTime,
                                                                               fineTotalIteration, materialParameters, blackHole);
        
        vector<vector<BlackHoleEulerStateVector> > coarseCurrentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(coarseCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                vector<double> topConservedVariableVector = coarseCurrentCellsWithBoundary[i][j + 1].computeConservedVariableVector(j * cellSpacing, (i - 1) * cellSpacing,
                                                                                                                                    0.0, materialParameters, blackHole);
                vector<double> leftConservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j].computeConservedVariableVector((j - 1) * cellSpacing, i * cellSpacing,
                                                                                                                                     0.0, materialParameters, blackHole);
                vector<double> conservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j + 1].computeConservedVariableVector(j * cellSpacing, i * cellSpacing, 0.0,
                                                                                                                                     materialParameters, blackHole);
                vector<double> rightConservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j + 2].computeConservedVariableVector((j + 1) * cellSpacing,
                                                                                                                                          i * cellSpacing, 0.0,
                                                                                                                                          materialParameters, blackHole);
                vector<double> bottomConservedVariableVector = coarseCurrentCellsWithBoundary[i + 2][j + 2].computeConservedVariableVector((j + 1) * cellSpacing,
                                                                                                                                           (i + 1) * cellSpacing, 0.0,
                                                                                                                                           materialParameters, blackHole);
                
                long conservedVariableCount = conservedVariableVector.size();
                
                for (int k = 0; k < conservedVariableCount; k++)
                {
                    if (abs((conservedVariableVector[k] - topConservedVariableVector[k]) / cellSpacing) >= AMRTolerance ||
                        abs((bottomConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= AMRTolerance ||
                        abs((conservedVariableVector[k] - leftConservedVariableVector[k]) / cellSpacing) >= AMRTolerance ||
                        abs((rightConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= AMRTolerance)
                    {
                        coarseAMRStructure[i][j] = false;
                        
                        fineAMRStructure[(i * 2)][(j * 2)] = true;
                        fineAMRStructure[(i * 2)][(j * 2) + 1] = true;
                        fineAMRStructure[(i * 2) + 1][(j * 2)] = true;
                        fineAMRStructure[(i * 2) + 1][(j * 2) + 1] = true;
                    }
                }
            }
        }
        
        if (order == 1)
        {
            computeXFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, coarseAMRStructure, materialParameters, blackHole);
            computeYFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, coarseAMRStructure, materialParameters, blackHole);
            computeXFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, coarseAMRStructure, materialParameters, blackHole);
        }
        else
        {
            computeXSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, 0.0, 1, coarseAMRStructure, materialParameters, blackHole);
            computeYSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, 0.0, 1, coarseAMRStructure, materialParameters, blackHole);
            computeXSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, 0.0, 1, coarseAMRStructure, materialParameters, blackHole);
        }
        
        AMRHelper::outputCoarseStatus(coarseCurrentIteration + 1, coarseCurrentTime + coarseTimeStep, coarseTimeStep);
        
        bool refined = true;
#pragma omp parallel for
        for (int i =0 ; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                if (fineAMRStructure[(i * 2)][(j * 2)] || fineAMRStructure[(i * 2)][(j * 2) + 1] || fineAMRStructure[(i * 2) + 1][(j * 2)] ||
                    fineAMRStructure[(i * 2) + 1][(j * 2) + 1])
                {
                    refined = true;
                }
            }
        }
        
        if (refined)
        {
            double fineCurrentTime = coarseCurrentTime;
            int fineCurrentIteration = 0;
            
            while (fineCurrentTime < (coarseCurrentTime + coarseTimeStep - pow(10.0, -4.0)))
            {
                double fineTimeStep = BlackHoleEulerSolvers::computeStableTimeStep2D(fineCurrentCells, cellSpacing * 0.5, CFLCoefficient, fineCurrentTime,
                                                                                     coarseCurrentTime + coarseTimeStep, fineTotalIteration, materialParameters, blackHole);
                
                if (order == 1)
                {
                    computeXFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep * 0.5, fineAMRStructure, materialParameters, blackHole);
                    computeYFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep, fineAMRStructure, materialParameters, blackHole);
                    computeXFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep * 0.5, fineAMRStructure, materialParameters, blackHole);
                }
                else
                {
                    computeXSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep * 0.5, 0.0, 1, fineAMRStructure, materialParameters, blackHole);
                    computeYSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep, 0.0, 1, fineAMRStructure, materialParameters, blackHole);
                    computeXSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep * 0.5, 0.0, 1, fineAMRStructure, materialParameters, blackHole);
                }
                
                fineCurrentTime += fineTimeStep;
                fineCurrentIteration += 1;
                fineTotalIteration += 1;
                
                AMRHelper::outputIntermediateStatus(fineCurrentIteration, fineCurrentTime, fineTimeStep);
            }
        }
        else
        {
            fineTotalIteration += 1;
        }
        
        coarseCurrentTime += coarseTimeStep;
        coarseCurrentIteration += 1;
        
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                if (!coarseAMRStructure[i][j])
                {
                    vector<double> topLeftFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2)].computeConservedVariableVector((j * 2) * (cellSpacing * 0.5),
                                                                                                                                          (i * 2) * (cellSpacing * 0.5),
                                                                                                                                          0.0, materialParameters,
                                                                                                                                          blackHole);
                    vector<double> topRightFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2) + 1].computeConservedVariableVector(
                        ((j * 2) + 1) * (cellSpacing * 0.5), (i * 2) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                    vector<double> bottomLeftFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2)].computeConservedVariableVector(
                        (j * 2) * (cellSpacing * 0.5), ((i * 2) + 1) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                    vector<double> bottomRightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2) + 1].computeConservedVariableVector(
                        ((j * 2) + 1) * (cellSpacing * 0.5), ((i * 2) + 1) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                    
                    coarseCurrentCells[i][j].setConservedVariableVector(VectorAlgebra::multiplyVector(0.25,
                        VectorAlgebra::addVectors(VectorAlgebra::addVectors(VectorAlgebra::addVectors(topLeftFineConservedVariableVector,
                                                                                                      topRightFineConservedVariableVector),
                                                                            bottomLeftFineConservedVariableVector), bottomRightFineConservedVariableVector)),
                                                                        j * cellSpacing, i * cellSpacing, 0.0, materialParameters, blackHole);
                }
                
                if (!fineAMRStructure[(i * 2)][(j * 2)])
                {
                    fineCurrentCells[(i * 2)][(j * 2)] = coarseCurrentCells[i][j];
                }
                if (!fineAMRStructure[(i * 2)][(j * 2) + 1])
                {
                    fineCurrentCells[(i * 2)][(j * 2) + 1] = coarseCurrentCells[i][j];
                }
                if (!fineAMRStructure[(i * 2) + 1][(j * 2)])
                {
                    fineCurrentCells[(i * 2) + 1][(j * 2)] = coarseCurrentCells[i][j];
                }
                if (!fineAMRStructure[(i * 2) + 1][(j * 2) + 1])
                {
                    fineCurrentCells[(i * 2) + 1][(j * 2) + 1] = coarseCurrentCells[i][j];
                }
            }
        }
        
        vector<vector<BlackHoleEulerStateVector> > fineCurrentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(fineCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                vector<double> topConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 1][(j * 2) + 2].computeConservedVariableVector(
                    ((j * 2) + 1) * (cellSpacing * 0.5), (i * 2) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                vector<double> leftConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 1].computeConservedVariableVector(
                    (j * 2) * (cellSpacing * 0.5), ((i * 2) + 1) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                vector<double> conservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 2].computeConservedVariableVector(
                    ((j * 2) + 1) * (cellSpacing * 0.5), ((i * 2) + 1) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                vector<double> rightConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 3].computeConservedVariableVector(
                    ((j * 2) + 2) * (cellSpacing * 0.5), ((i * 2) + 1) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                vector<double> bottomConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 3][(j * 2) + 2].computeConservedVariableVector(
                    ((j * 2) + 1) * (cellSpacing * 0.5), ((i * 2) + 2) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                
                long conservedVariableCount = conservedVariableVector.size();
                bool markedCell = false;
                
                for (int k = 0; k < conservedVariableCount; k++)
                {
                    if (abs((conservedVariableVector[k] - topConservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 0.5) ||
                        abs((bottomConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 0.5) ||
                        abs((conservedVariableVector[k] - leftConservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 0.5) ||
                        abs((rightConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 0.5))
                    {
                        markedCell = true;
                    }
                }
                
                if (!markedCell && fineAMRStructure[(i * 2)][(j * 2)])
                {
                    fineAMRStructure[(i * 2)][(j * 2)] = false;
                    fineAMRStructure[(i * 2)][(j * 2) + 1] = false;
                    fineAMRStructure[(i * 2) + 1][(j * 2)] = false;
                    fineAMRStructure[(i * 2) + 1][(j * 2) + 1] = false;
                    
                    coarseAMRStructure[i][j] = true;
                    
                    vector<double> topLeftFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2)].computeConservedVariableVector((j * 2) * (cellSpacing * 0.5),
                                                                                                                                          (i * 2) * (cellSpacing * 0.5),
                                                                                                                                          0.0, materialParameters,
                                                                                                                                          blackHole);
                    vector<double> topRightFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2) + 1].computeConservedVariableVector(
                        ((j * 2) + 1) * (cellSpacing * 0.5), (i * 2) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                    vector<double> bottomLeftFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2)].computeConservedVariableVector(
                        (j * 2) * (cellSpacing * 0.5), ((i * 2) + 1) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                    vector<double> bottomRightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2) + 1].computeConservedVariableVector(
                        ((j * 2) + 1) * (cellSpacing * 0.5), ((i * 2) + 1) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                    
                    coarseCurrentCells[i][j].setConservedVariableVector(VectorAlgebra::multiplyVector(0.25,
                        VectorAlgebra::addVectors(VectorAlgebra::addVectors(VectorAlgebra::addVectors(topLeftFineConservedVariableVector,
                                                                                                      topRightFineConservedVariableVector),
                                                                            bottomLeftFineConservedVariableVector), bottomRightFineConservedVariableVector)),
                                                                        j * cellSpacing, i * cellSpacing, 0.0, materialParameters, blackHole);
                }
            }
        }
    }
    
    return tuple<vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >, vector<vector<bool> >,
    vector<vector<bool> >>(coarseCurrentCells, fineCurrentCells, coarseAMRStructure, fineAMRStructure);
}

tuple<vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<bool>, vector<bool>, vector<bool>>
BlackHoleEulerAMR::solveLevel2AMR(vector<BlackHoleEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance,
                                  int order, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long cellCount = initialCells.size();
    
    vector<BlackHoleEulerStateVector> coarseCurrentCells(cellCount);
    vector<BlackHoleEulerStateVector> intermediateCurrentCells(cellCount * 2);
    vector<BlackHoleEulerStateVector> fineCurrentCells(cellCount * 4);
    
    vector<bool> coarseAMRStructure(cellCount);
    vector<bool> intermediateAMRStructure(cellCount * 2);
    vector<bool> fineAMRStructure(cellCount * 4);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        coarseAMRStructure[i] = true;
        
        intermediateAMRStructure[(i * 2)] = false;
        intermediateAMRStructure[(i * 2) + 1] = false;
        
        for (int j = 0; j < 4; j++)
        {
            fineAMRStructure[(i * 4) + j] = false;
        }
    }
    
    double coarseCurrentTime = 0.0;
    int coarseCurrentIteration = 0;
    int fineTotalIteration = 0;
    coarseCurrentCells = initialCells;
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        intermediateCurrentCells[(i * 2)] = coarseCurrentCells[i];
        intermediateCurrentCells[(i * 2) + 1] = coarseCurrentCells[i];
        
        for (int j = 0; j < 4; j++)
        {
            fineCurrentCells[(i * 4) + j] = coarseCurrentCells[i];
        }
    }
    
    while (coarseCurrentTime < finalTime)
    {
        double coarseTimeStep = BlackHoleEulerSolvers::computeStableTimeStep(coarseCurrentCells, cellSpacing, CFLCoefficient, coarseCurrentTime, finalTime,
                                                                             fineTotalIteration, materialParameters, blackHole);
        
        vector<BlackHoleEulerStateVector> coarseCurrentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(coarseCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            vector<double> leftConservedVariableVector = coarseCurrentCellsWithBoundary[i].computeConservedVariableVector((i - 1) * cellSpacing, 0.0, 0.0,
                                                                                                                          materialParameters, blackHole);
            vector<double> conservedVariableVector = coarseCurrentCellsWithBoundary[i + 1].computeConservedVariableVector(i * cellSpacing, 0.0, 0.0,
                                                                                                                      materialParameters, blackHole);
            vector<double> rightConservedVariableVector = coarseCurrentCellsWithBoundary[i + 2].computeConservedVariableVector((i + 1) * cellSpacing, 0.0, 0.0,
                                                                                                                           materialParameters, blackHole);
            
            long conservedVariableCount = conservedVariableVector.size();
            
            for (int j = 0; j < conservedVariableCount; j++)
            {
                if (abs((conservedVariableVector[j] - leftConservedVariableVector[j]) / cellSpacing) >= AMRTolerance ||
                    abs((rightConservedVariableVector[j] - conservedVariableVector[j]) / cellSpacing) >= AMRTolerance)
                {
                    coarseAMRStructure[i] = false;
                    
                    intermediateAMRStructure[(i * 2)] = true;
                    intermediateAMRStructure[(i * 2) + 1] = true;
                }
            }
        }
        
        if (order == 1)
        {
            computeFORCETimeStepAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, coarseAMRStructure, materialParameters, blackHole);
        }
        else
        {
            computeSLICTimeStepAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, 0.0, 1, coarseAMRStructure, materialParameters, blackHole);
        }
        
        AMRHelper::outputCoarseStatus(coarseCurrentIteration + 1, coarseCurrentTime + coarseTimeStep, coarseTimeStep);
        
        bool refined = false;
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            if (intermediateAMRStructure[(i * 2)] || intermediateAMRStructure[(i * 2) + 1])
            {
                refined = true;
            }
        }
        
        if (refined)
        {
            double intermediateCurrentTime = coarseCurrentTime;
            int intermediateCurrentIteration = 0;
            
            while (intermediateCurrentTime < (coarseCurrentTime + coarseTimeStep - pow(10.0, -4.0)))
            {
                double intermediateTimeStep = BlackHoleEulerSolvers::computeStableTimeStep(intermediateCurrentCells, cellSpacing * 0.5, CFLCoefficient,
                                                                                           intermediateCurrentTime, coarseCurrentTime + coarseTimeStep,
                                                                                           fineTotalIteration, materialParameters, blackHole);
                
                vector<BlackHoleEulerStateVector> intermediateCurrentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(intermediateCurrentCells, 1);
                
#pragma omp parallel for
                for (int i = 0; i < (cellCount * 2); i++)
                {
                    vector<double> leftConservedVariableVector = intermediateCurrentCellsWithBoundary[i].computeConservedVariableVector((i - 1) * (0.5 * cellSpacing), 0.0,
                                                                                                                                        0.0, materialParameters, blackHole);
                    vector<double> conservedVariableVector = intermediateCurrentCellsWithBoundary[i + 1].computeConservedVariableVector(i * (0.5 * cellSpacing), 0.0, 0.0,
                                                                                                                                        materialParameters, blackHole);
                    vector<double> rightConservedVariableVector = intermediateCurrentCellsWithBoundary[i + 2].computeConservedVariableVector((i + 1) * (0.5 * cellSpacing),
                                                                                                                                             0.0, 0.0, materialParameters,
                                                                                                                                             blackHole);
                    
                    long conservedVariableCount = conservedVariableVector.size();
                    
                    for (int j = 0; j < conservedVariableCount; j++)
                    {
                        if ((abs((conservedVariableVector[j] - leftConservedVariableVector[j]) / cellSpacing) >= (AMRTolerance * 2.0 * order) ||
                             abs((rightConservedVariableVector[j] - conservedVariableVector[j]) / cellSpacing) >= (AMRTolerance * 2.0 * order)) &&
                            intermediateAMRStructure[i])
                        {
                            intermediateAMRStructure[i] = false;
                            
                            fineAMRStructure[(i * 2)] = true;
                            fineAMRStructure[(i * 2) + 1] = true;
                        }
                    }
                }
                
                if (order == 1)
                {
                    computeFORCETimeStepAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep, intermediateAMRStructure, materialParameters, blackHole);
                }
                else
                {
                    computeSLICTimeStepAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep, 0.0, 1, intermediateAMRStructure,
                                           materialParameters, blackHole);
                }
                
                AMRHelper::outputIntermediateStatus(intermediateCurrentIteration + 1, intermediateCurrentTime + intermediateTimeStep, intermediateTimeStep);
                
                bool refinedLevel2 = false;
#pragma omp parallel for
                for (int i = 0; i < (cellCount * 2); i++)
                {
                    if (fineAMRStructure[(i * 2)] || fineAMRStructure[(i * 2) + 1])
                    {
                        refinedLevel2 = true;
                    }
                }
                
                if (refinedLevel2)
                {
                    double fineCurrentTime = intermediateCurrentTime;
                    int fineCurrentIteration = 0;
                    
                    while (fineCurrentTime < (intermediateCurrentTime + intermediateTimeStep - pow(10.0, -4.0)))
                    {
                        double fineTimeStep = BlackHoleEulerSolvers::computeStableTimeStep(fineCurrentCells, cellSpacing * 0.25, CFLCoefficient, fineCurrentTime,
                                                                                           intermediateCurrentTime + intermediateTimeStep, fineTotalIteration,
                                                                                           materialParameters, blackHole);
                        
                        if (order == 1)
                        {
                            computeFORCETimeStepAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep, fineAMRStructure, materialParameters, blackHole);
                        }
                        else
                        {
                            computeSLICTimeStepAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep, 0.0, 1, fineAMRStructure, materialParameters, blackHole);
                        }
                        
                        fineCurrentTime += fineTimeStep;
                        fineCurrentIteration += 1;
                        fineTotalIteration += 1;
                        
                        AMRHelper::outputFineStatus(fineCurrentIteration, fineCurrentTime, fineTimeStep);
                    }
                }
                else
                {
                    fineTotalIteration += 1;
                }
                
                intermediateCurrentTime += intermediateTimeStep;
                intermediateCurrentIteration += 1;
                fineTotalIteration += 1;
                
#pragma omp parallel for
                for (int i = 0; i < (cellCount * 2); i++)
                {
                    if (!intermediateAMRStructure[i] && (fineAMRStructure[(i * 2)] && fineAMRStructure[(i * 2) + 1]))
                    {
                        vector<double> leftFineConservedVariableVector = fineCurrentCells[(i * 2)].computeConservedVariableVector((i * 2) * (cellSpacing * 0.25), 0.0,
                                                                                                                                  0.0, materialParameters, blackHole);
                        vector<double> rightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1].computeConservedVariableVector(((i * 2) + 1) * (cellSpacing * 0.25),
                                                                                                                                       0.0, 0.0, materialParameters,
                                                                                                                                       blackHole);
                        
                        intermediateCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5,
                                                                                                             VectorAlgebra::addVectors(leftFineConservedVariableVector,
                                                                                                                                       rightFineConservedVariableVector)),
                                                                               i * (cellSpacing * 0.5), 0.0, 0.0, materialParameters, blackHole);
                    }
                    
                    if (!fineAMRStructure[(i * 2)] && intermediateAMRStructure[i])
                    {
                        fineCurrentCells[(i * 2)] = intermediateCurrentCells[i];
                    }
                    if (!fineAMRStructure[(i * 2) + 1] && intermediateAMRStructure[i])
                    {
                        fineCurrentCells[(i * 2) + 1] = intermediateCurrentCells[i];
                    }
                }
                
                vector<BlackHoleEulerStateVector> fineCurrentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(fineCurrentCells, 1);
                
#pragma omp parallel for
                for (int i = 0; i < (cellCount * 2); i++)
                {
                    vector<double> leftConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 1].computeConservedVariableVector((i * 2) * (cellSpacing * 0.25),
                                                                                                                                           0.0, 0.0, materialParameters,
                                                                                                                                           blackHole);
                    vector<double> conservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2].computeConservedVariableVector(((i * 2) + 1) * (cellSpacing * 0.25),
                                                                                                                                      0.0, 0.0, materialParameters,
                                                                                                                                      blackHole);
                    vector<double> rightConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 3].computeConservedVariableVector(((i * 2) + 2) *
                                                                                                                                           (cellSpacing * 0.25), 0.0, 0.0,
                                                                                                                                           materialParameters, blackHole);
                    
                    long conservedVariableCount = conservedVariableVector.size();
                    bool markedCell = false;
                    
                    for (int j = 0; j < conservedVariableCount; j++)
                    {
                        if (abs((conservedVariableVector[j] - leftConservedVariableVector[j]) / cellSpacing) >= AMRTolerance ||
                            abs((rightConservedVariableVector[j] - conservedVariableVector[j]) / cellSpacing) >= AMRTolerance)
                        {
                            markedCell = true;
                        }
                    }
                    
                    if (!markedCell && (fineAMRStructure[(i * 2)] && fineAMRStructure[(i * 2) + 1]))
                    {
                        fineAMRStructure[(i * 2)] = false;
                        fineAMRStructure[(i * 2) + 1] = false;
                        
                        intermediateAMRStructure[i] = true;
                        
                        vector<double> leftFineConservedVariableVector = fineCurrentCells[(i * 2)].computeConservedVariableVector((i * 2) * (cellSpacing * 0.25), 0.0,
                                                                                                                                  0.0, materialParameters, blackHole);
                        vector<double> rightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1].computeConservedVariableVector(((i * 2) + 1) * (cellSpacing * 0.25),
                                                                                                                                       0.0, 0.0, materialParameters,
                                                                                                                                       blackHole);
                        
                        intermediateCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5,
                                                                                                             VectorAlgebra::addVectors(leftFineConservedVariableVector,
                                                                                                                                       rightFineConservedVariableVector)),
                                                                               i * (cellSpacing * 0.5), 0.0, 0.0, materialParameters, blackHole);
                    }
                }
            }
        }
        else
        {
            fineTotalIteration += 1;
        }
        
        coarseCurrentTime += coarseTimeStep;
        coarseCurrentIteration += 1;
        
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            if (!coarseAMRStructure[i])
            {
                vector<double> leftIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2)].computeConservedVariableVector((i * 2) * (cellSpacing * 0.5),
                                                                                                                                          0.0, 0.0, materialParameters,
                                                                                                                                          blackHole);
                vector<double> rightIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2) + 1].computeConservedVariableVector(((i * 2) + 1) *
                                                                                                                                               (cellSpacing * 0.5), 0.0,
                                                                                                                                               0.0, materialParameters,
                                                                                                                                               blackHole);
                
                coarseCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5,
                                                                                               VectorAlgebra::addVectors(leftIntermediateConservedVariableVector,
                                                                                                                         rightIntermediateConservedVariableVector)),
                                                                 i * cellSpacing, 0.0, 0.0, materialParameters, blackHole);
            }
            
            if (!intermediateAMRStructure[(i * 2)])
            {
                intermediateCurrentCells[(i * 2)] = coarseCurrentCells[i];
            }
            if (!intermediateAMRStructure[(i * 2) + 1])
            {
                intermediateCurrentCells[(i * 2) + 1] = coarseCurrentCells[i];
            }
        }
        
        vector<BlackHoleEulerStateVector> intermediateCurrentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(intermediateCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            vector<double> leftConservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 1].computeConservedVariableVector((i * 2) *
                                                                                                                                          (cellSpacing * 0.5), 0.0, 0.0,
                                                                                                                                          materialParameters, blackHole);
            vector<double> conservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 2].computeConservedVariableVector(((i * 2) + 1) * (cellSpacing * 0.5),
                                                                                                                                      0.0, 0.0, materialParameters,
                                                                                                                                      blackHole);
            vector<double> rightConservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 3].computeConservedVariableVector(((i * 2) + 2) *
                                                                                                                                           (cellSpacing * 0.5), 0.0, 0.0,
                                                                                                                                           materialParameters, blackHole);
            
            long conservedVariableCount = conservedVariableVector.size();
            bool markedCell = false;
            
            for (int j = 0; j < conservedVariableCount; j++)
            {
                if (abs((conservedVariableVector[j] - leftConservedVariableVector[j]) / cellSpacing) >= (AMRTolerance * 0.5 * (1.0 / order)) ||
                    abs((rightConservedVariableVector[j] - conservedVariableVector[j]) / cellSpacing) >= (AMRTolerance * 0.5 * (1.0 / order)))
                {
                    markedCell = true;
                }
            }
            
            if (!markedCell && (intermediateAMRStructure[(i * 2)] && intermediateAMRStructure[(i * 2) + 1]))
            {
                intermediateAMRStructure[(i * 2)] = false;
                intermediateAMRStructure[(i * 2) + 1] = false;
                
                coarseAMRStructure[i] = true;
                
                vector<double> leftIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2)].computeConservedVariableVector((i * 2) * (cellSpacing * 0.5),
                                                                                                                                          0.0, 0.0, materialParameters,
                                                                                                                                          blackHole);
                vector<double> rightIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2) + 1].computeConservedVariableVector(((i * 2) + 1) *
                                                                                                                                               (cellSpacing * 0.5), 0.0,
                                                                                                                                               0.0, materialParameters,
                                                                                                                                               blackHole);
                
                coarseCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5,
                                                                                               VectorAlgebra::addVectors(leftIntermediateConservedVariableVector,
                                                                                                                         rightIntermediateConservedVariableVector)),
                                                                 i * cellSpacing, 0.0, 0.0, materialParameters, blackHole);
            }
        }
    }
    
    return tuple<vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<bool>, vector<bool>,
    vector<bool>>(coarseCurrentCells, intermediateCurrentCells, fineCurrentCells, coarseAMRStructure, intermediateAMRStructure, fineAMRStructure);
}

tuple<vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >,
vector<vector<bool> >, vector<vector<bool> >, vector<vector<bool> >>
BlackHoleEulerAMR::solve2DLevel2AMR(vector<vector<BlackHoleEulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                    double AMRTolerance, int order, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long rowCount = initialCells.size();
    long columnCount = initialCells[0].size();
    
    vector<vector<BlackHoleEulerStateVector> > coarseCurrentCells(rowCount, vector<BlackHoleEulerStateVector>(columnCount));
    vector<vector<BlackHoleEulerStateVector> > intermediateCurrentCells(rowCount * 2, vector<BlackHoleEulerStateVector>(columnCount * 2));
    vector<vector<BlackHoleEulerStateVector> > fineCurrentCells(rowCount * 4, vector<BlackHoleEulerStateVector>(columnCount * 4));
    
    vector<vector<bool> > coarseAMRStructure(rowCount, vector<bool>(columnCount));
    vector<vector<bool> > intermediateAMRStructure(rowCount * 2, vector<bool>(columnCount * 2));
    vector<vector<bool> > fineAMRStructure(rowCount * 4, vector<bool>(columnCount * 4));
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            coarseAMRStructure[i][j] = true;
            
            intermediateAMRStructure[(i * 2)][(j * 2)] = false;
            intermediateAMRStructure[(i * 2)][(j * 2) + 1] = false;
            intermediateAMRStructure[(i * 2) + 1][(j * 2)] = false;
            intermediateAMRStructure[(i * 2) + 1][(j * 2) + 1] = false;
            
            for (int k = 0; k < 4; k++)
            {
                for (int l = 0; l < 4; l++)
                {
                    fineAMRStructure[(i * 4) + k][(j * 4) + l] = false;
                }
            }
        }
    }
    
    double coarseCurrentTime = 0.0;
    int coarseCurrentIteration = 0;
    int fineTotalIteration = 0;
    coarseCurrentCells = initialCells;
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            intermediateCurrentCells[(i * 2)][(j * 2)] = coarseCurrentCells[i][j];
            intermediateCurrentCells[(i * 2)][(j * 2) + 1] = coarseCurrentCells[i][j];
            intermediateCurrentCells[(i * 2) + 1][(j * 2)] = coarseCurrentCells[i][j];
            intermediateCurrentCells[(i * 2) + 1][(j * 2) + 1] = coarseCurrentCells[i][j];
            
            for (int k = 0; k < 4; k++)
            {
                for (int l = 0; l < 4; l++)
                {
                    fineCurrentCells[(i * 4) + k][(j * 4) + l] = coarseCurrentCells[i][j];
                }
            }
        }
    }
    
    while (coarseCurrentTime < finalTime)
    {
        double coarseTimeStep = BlackHoleEulerSolvers::computeStableTimeStep2D(coarseCurrentCells, cellSpacing, CFLCoefficient, coarseCurrentTime, finalTime,
                                                                               fineTotalIteration, materialParameters, blackHole);
        
        vector<vector<BlackHoleEulerStateVector> > coarseCurrentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(coarseCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                vector<double> topConservedVariableVector = coarseCurrentCellsWithBoundary[i][j + 1].computeConservedVariableVector(j * cellSpacing, (i - 1) * cellSpacing,
                                                                                                                                    0.0, materialParameters, blackHole);
                vector<double> leftConservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j].computeConservedVariableVector((j - 1) * cellSpacing, i * cellSpacing,
                                                                                                                                     0.0, materialParameters, blackHole);
                vector<double> conservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j + 1].computeConservedVariableVector(j * cellSpacing, i * cellSpacing, 0.0,
                                                                                                                                     materialParameters, blackHole);
                vector<double> rightConservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j + 2].computeConservedVariableVector((j + 1) * cellSpacing,
                                                                                                                                          i * cellSpacing, 0.0,
                                                                                                                                          materialParameters, blackHole);
                vector<double> bottomConservedVariableVector = coarseCurrentCellsWithBoundary[i + 2][j + 1].computeConservedVariableVector(j * cellSpacing,
                                                                                                                                           (i + 1) * cellSpacing, 0.0,
                                                                                                                                           materialParameters, blackHole);
                
                long conservedVariableCount = conservedVariableVector.size();
                
                for (int k = 0; k < conservedVariableCount; k++)
                {
                    if (abs((conservedVariableVector[k] - topConservedVariableVector[k]) / cellSpacing) >= AMRTolerance ||
                        abs((bottomConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= AMRTolerance ||
                        abs((conservedVariableVector[k] - leftConservedVariableVector[k]) / cellSpacing) >= AMRTolerance ||
                        abs((rightConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= AMRTolerance)
                    {
                        coarseAMRStructure[i][j] = false;
                        
                        intermediateAMRStructure[(i * 2)][(j * 2)] = true;
                        intermediateAMRStructure[(i * 2)][(j * 2) + 1] = true;
                        intermediateAMRStructure[(i * 2) + 1][(j * 2)] = true;
                        intermediateAMRStructure[(i * 2) + 1][(j * 2) + 1] = true;
                    }
                }
            }
        }
        
        if (order == 1)
        {
            computeXFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, coarseAMRStructure, materialParameters, blackHole);
            computeYFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, coarseAMRStructure, materialParameters, blackHole);
            computeXFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, coarseAMRStructure, materialParameters, blackHole);
        }
        else
        {
            computeXSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, 0.0, 1, coarseAMRStructure, materialParameters, blackHole);
            computeYSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, 0.0, 1, coarseAMRStructure, materialParameters, blackHole);
            computeXSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, 0.0, 1, coarseAMRStructure, materialParameters, blackHole);
        }
        
        AMRHelper::outputCoarseStatus(coarseCurrentIteration + 1, coarseCurrentTime + coarseTimeStep, coarseTimeStep);
        
        bool refined = false;
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                if (intermediateAMRStructure[(i * 2)][(j * 2)] || intermediateAMRStructure[(i * 2)][(j * 2) + 1] ||
                    intermediateAMRStructure[(i * 2) + 1][(j * 2)] || intermediateAMRStructure[(i * 2) + 1][(j * 2) + 1])
                {
                    refined = true;
                }
            }
        }
        
        if (refined)
        {
            double intermediateCurrentTime = coarseCurrentTime;
            int intermediateCurrentIteration = 0;
            
            while (intermediateCurrentTime < (coarseCurrentTime + coarseTimeStep - pow(10.0, -4.0)))
            {
                double intermediateTimeStep = BlackHoleEulerSolvers::computeStableTimeStep2D(intermediateCurrentCells, cellSpacing * 0.5, CFLCoefficient,
                                                                                             intermediateCurrentTime, coarseCurrentTime + coarseTimeStep,
                                                                                             fineTotalIteration, materialParameters, blackHole);
                
                vector<vector<BlackHoleEulerStateVector> > intermediateCurrentCellsWithBoundary =
                BlackHoleEulerSolvers::insertBoundaryCells2D(intermediateCurrentCells, 1);
                
#pragma omp parallel for
                for (int i = 0; i < (rowCount * 2); i++)
                {
                    for (int j = 0; j < (columnCount * 2); j++)
                    {
                        vector<double> topConservedVariableVector = intermediateCurrentCellsWithBoundary[i][j + 1].
                        computeConservedVariableVector(j * (cellSpacing * 0.5), (i - 1) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                        vector<double> leftConservedVariableVector = intermediateCurrentCellsWithBoundary[i + 1][j].
                        computeConservedVariableVector((j - 1) * (cellSpacing * 0.5), i * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                        vector<double> conservedVariableVector = intermediateCurrentCellsWithBoundary[i + 1][j + 1].
                        computeConservedVariableVector(j * (cellSpacing * 0.5), i * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                        vector<double> rightConservedVariableVector = intermediateCurrentCellsWithBoundary[i + 1][j + 2].
                        computeConservedVariableVector((j + 1) * (cellSpacing * 0.5), i * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                        vector<double> bottomConservedVariableVector = intermediateCurrentCellsWithBoundary[i + 2][j + 1].
                        computeConservedVariableVector(j * (cellSpacing * 0.5), (i + 1) * (cellSpacing * 0.5), 0.0, materialParameters, blackHole);
                        
                        long conservedVariableCount = conservedVariableVector.size();
                        
                        for (int k = 0; k < conservedVariableCount; k++)
                        {
                            if ((abs((conservedVariableVector[k] - topConservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 2.0 * order) ||
                                 abs((bottomConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 2.0 * order) ||
                                 abs((conservedVariableVector[k] - leftConservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 2.0 * order) ||
                                 abs((rightConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 2.0 * order)) &&
                                intermediateAMRStructure[i][j])
                            {
                                intermediateAMRStructure[i][j] = false;
                                
                                fineAMRStructure[(i * 2)][(j * 2)] = true;
                                fineAMRStructure[(i * 2)][(j * 2) + 1] = true;
                                fineAMRStructure[(i * 2) + 1][(j * 2)] = true;
                                fineAMRStructure[(i * 2) + 1][(j * 2) + 1] = true;
                            }
                        }
                    }
                }
                
                if (order == 1)
                {
                    computeXFORCETimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep * 0.5, intermediateAMRStructure, materialParameters,
                                               blackHole);
                    computeYFORCETimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep, intermediateAMRStructure, materialParameters, blackHole);
                    computeXFORCETimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep * 0.5, intermediateAMRStructure, materialParameters,
                                               blackHole);
                }
                else
                {
                    computeXSLICTimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep * 0.5, 0.0, 1, intermediateAMRStructure,
                                              materialParameters, blackHole);
                    computeYSLICTimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep, 0.0, 1, intermediateAMRStructure, materialParameters,
                                              blackHole);
                    computeXSLICTimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep * 0.5, 0.0, 1, intermediateAMRStructure,
                                              materialParameters, blackHole);
                }
                
                AMRHelper::outputIntermediateStatus(intermediateCurrentIteration + 1, intermediateCurrentTime + intermediateTimeStep, intermediateTimeStep);
                
                bool refinedLevel2 = false;
#pragma omp parallel for
                for (int i = 0; i < (rowCount * 2); i++)
                {
                    for (int j = 0; j < (columnCount * 2); j++)
                    {
                        if (fineAMRStructure[(i * 2)][(j * 2)] || fineAMRStructure[(i * 2)][(j * 2) + 1] || fineAMRStructure[(i * 2) + 1][(j * 2)] ||
                            fineAMRStructure[(i * 2) + 1][(j * 2) + 1])
                        {
                            refinedLevel2 = true;
                        }
                    }
                }
                
                if (refinedLevel2)
                {
                    double fineCurrentTime = intermediateCurrentTime;
                    int fineCurrentIteration = 0;
                    
                    while (fineCurrentTime < (intermediateCurrentTime + intermediateTimeStep - pow(10.0, -4.0)))
                    {
                        double fineTimeStep = BlackHoleEulerSolvers::computeStableTimeStep2D(fineCurrentCells, cellSpacing * 0.25, CFLCoefficient, fineCurrentTime,
                                                                                             intermediateCurrentTime + intermediateTimeStep, fineTotalIteration,
                                                                                             materialParameters, blackHole);
                        
                        if (order == 1)
                        {
                            computeXFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep * 0.5, fineAMRStructure, materialParameters, blackHole);
                            computeYFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep, fineAMRStructure, materialParameters, blackHole);
                            computeXFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep * 0.5, fineAMRStructure, materialParameters, blackHole);
                        }
                        else
                        {
                            computeXSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep * 0.5, 0.0, 1, fineAMRStructure, materialParameters, blackHole);
                            computeYSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep, 0.0, 1, fineAMRStructure, materialParameters, blackHole);
                            computeXSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep * 0.5, 0.0, 1, fineAMRStructure, materialParameters, blackHole);
                        }
                        
                        fineCurrentTime += fineTimeStep;
                        fineCurrentIteration += 1;
                        fineTotalIteration += 1;
                        
                        AMRHelper::outputFineStatus(fineCurrentIteration, fineCurrentTime, fineTimeStep);
                    }
                }
                else
                {
                    fineTotalIteration += 1;
                }
                
                intermediateCurrentTime += intermediateTimeStep;
                intermediateCurrentIteration += 1;
                fineTotalIteration += 1;
                
#pragma omp parallel for
                for (int i = 0; i < (rowCount * 2); i++)
                {
                    for (int j = 0; j < (columnCount * 2); j++)
                    {
                        if (!intermediateAMRStructure[i][j])
                    }
                }
            }
        }
    }
    
    return tuple<vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >,
    vector<vector<bool> >, vector<vector<bool> >, vector<vector<bool> >>(coarseCurrentCells, intermediateCurrentCells, fineCurrentCells, coarseAMRStructure,
                                                                         intermediateAMRStructure, fineAMRStructure);
}
