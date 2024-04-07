#include "EulerAMR.hpp"

void EulerAMR::computeFORCETimeStepAMR(vector<EulerStateVector> & currentCells, double cellSpacing, double timeStep, vector<bool> AMRStructure,
                                       EulerMaterialParameters materialParameters)
{
    long cellCount = currentCells.size();
    vector<EulerStateVector> currentCellsWithBoundary = EulerSolvers::insertBoundaryCells(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        if (AMRStructure[i])
        {
            vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);
            
            vector<double> leftFluxVector = EulerFirstOrderSolver::computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing,
                                                                                     timeStep, materialParameters);
            vector<double> rightFluxVector = EulerFirstOrderSolver::computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing,
                                                                                      timeStep, materialParameters);
            
            currentCells[i].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                 cellSpacing, timeStep), materialParameters);
        }
    }
}

void EulerAMR::computeXFORCETimeStep2DAMR(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                          vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<EulerStateVector> > currentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (AMRStructure[i][j])
            {
                vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);
                
                vector<double> leftFluxVector = EulerFirstOrderSolver::computeXFORCEFlux(currentCellsWithBoundary[i + 1][j], currentCellsWithBoundary[i + 1][j + 1],
                                                                                         cellSpacing, timeStep, materialParameters);
                vector<double> rightFluxVector = EulerFirstOrderSolver::computeXFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 1][j + 2],
                                                                                          cellSpacing, timeStep, materialParameters);
                
                currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                        cellSpacing, timeStep), materialParameters);
            }
        }
    }
}

void EulerAMR::computeYFORCETimeStep2DAMR(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, vector<vector<bool> > AMRStructure,
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
            if (AMRStructure[i][j])
            {
                vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);
                
                vector<double> topFluxVector = EulerFirstOrderSolver::computeYFORCEFlux(currentCellsWithBoundary[i][j + 1], currentCellsWithBoundary[i + 1][j + 1],
                                                                                        cellSpacing, timeStep, materialParameters);
                vector<double> bottomFluxVector = EulerFirstOrderSolver::computeYFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 2][j + 1],
                                                                                           cellSpacing, timeStep, materialParameters);
                
                currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector,
                                                                                                        cellSpacing, timeStep), materialParameters);
            }
        }
    }
}

void EulerAMR::computeSLICTimeStepAMR(vector<EulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                      vector<bool> AMRStructure, EulerMaterialParameters materialParameters)
{
    long cellCount = currentCells.size();
    vector<EulerStateVector> currentCellsWithBoundary = EulerSolvers::insertBoundaryCells(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        if (AMRStructure[i])
        {
            vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);
            
            vector<double> leftFluxVector = EulerSecondOrderSolver::computeXSLICFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1],
                                                                                     currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3], cellSpacing,
                                                                                     timeStep, bias, slopeLimiter, materialParameters);
            vector<double> rightFluxVector = EulerSecondOrderSolver::computeXSLICFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2],
                                                                                      currentCellsWithBoundary[i + 3], currentCellsWithBoundary[i + 4], cellSpacing,
                                                                                      timeStep, bias, slopeLimiter, materialParameters);
            
            currentCells[i].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                 cellSpacing, timeStep), materialParameters);
        }
    }
}

void EulerAMR::computeXSLICTimeStep2DAMR(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                         vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<EulerStateVector> > currentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (AMRStructure[i][j])
            {
                vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);
                
                vector<double> leftFluxVector = EulerSecondOrderSolver::computeXSLICFlux(currentCellsWithBoundary[i + 2][j], currentCellsWithBoundary[i + 2][j + 1],
                                                                                         currentCellsWithBoundary[i + 2][j + 2], currentCellsWithBoundary[i + 2][j + 3],
                                                                                         cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
                vector<double> rightFluxVector = EulerSecondOrderSolver::computeXSLICFlux(currentCellsWithBoundary[i + 2][j + 1], currentCellsWithBoundary[i + 2][j + 2],
                                                                                          currentCellsWithBoundary[i + 2][j + 3], currentCellsWithBoundary[i + 2][j + 4],
                                                                                          cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
                
                currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                                        cellSpacing, timeStep), materialParameters);
            }
        }
    }
}

void EulerAMR::computeYSLICTimeStep2DAMR(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                         vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<EulerStateVector> > currentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(currentCells, 2);
    
#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (AMRStructure[i][j])
            {
                vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);
                
                vector<double> topFluxVector = EulerSecondOrderSolver::computeYSLICFlux(currentCellsWithBoundary[i][j + 2], currentCellsWithBoundary[i + 1][j + 2],
                                                                                        currentCellsWithBoundary[i + 2][j + 2], currentCellsWithBoundary[i + 3][j + 2],
                                                                                        cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
                vector<double> bottomFluxVector = EulerSecondOrderSolver::computeYSLICFlux(currentCellsWithBoundary[i + 1][j + 2], currentCellsWithBoundary[i + 2][j + 2],
                                                                                           currentCellsWithBoundary[i + 3][j + 2], currentCellsWithBoundary[i + 4][j + 2],
                                                                                           cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
                
                currentCells[i][j].setConservedVariableVector(EulerFirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector,
                                                                                                        cellSpacing, timeStep), materialParameters);
            }
        }
    }
}

tuple<vector<EulerStateVector>, vector<EulerStateVector>, vector<bool>, vector<bool>> EulerAMR::solveLevel1AMR(vector<EulerStateVector> initialCells,
                                                                                                               double cellSpacing, double CFLCoefficient,
                                                                                                               double finalTime, double AMRTolerance, int order,
                                                                                                               EulerMaterialParameters materialParameters)
{
    long cellCount = initialCells.size();
    
    vector<EulerStateVector> coarseCurrentCells(cellCount);
    vector<EulerStateVector> fineCurrentCells(cellCount * 2);
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
        double coarseTimeStep = EulerSolvers::computeStableTimeStep(coarseCurrentCells, cellSpacing, CFLCoefficient, coarseCurrentTime, finalTime,
                                                                    fineTotalIteration, materialParameters);
        
        vector<EulerStateVector> coarseCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells(coarseCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            vector<double> leftConservedVariableVector = coarseCurrentCellsWithBoundary[i].computeConservedVariableVector(materialParameters);
            vector<double> conservedVariableVector = coarseCurrentCellsWithBoundary[i + 1].computeConservedVariableVector(materialParameters);
            vector<double> rightConservedVariableVector = coarseCurrentCellsWithBoundary[i + 2].computeConservedVariableVector(materialParameters);
            
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
            computeFORCETimeStepAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, coarseAMRStructure, materialParameters);
        }
        else
        {
            computeSLICTimeStepAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, 0.0, 1, coarseAMRStructure, materialParameters);
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
                double fineTimeStep = EulerSolvers::computeStableTimeStep(fineCurrentCells, cellSpacing * 0.5, CFLCoefficient, fineCurrentTime,
                                                                          coarseCurrentTime + coarseTimeStep, fineTotalIteration, materialParameters);
                
                if (order == 1)
                {
                    computeFORCETimeStepAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep, fineAMRStructure, materialParameters);
                }
                else
                {
                    computeSLICTimeStepAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep, 0.0, 1, fineAMRStructure, materialParameters);
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
                vector<double> leftFineConservedVariableVector = fineCurrentCells[(i * 2)].computeConservedVariableVector(materialParameters);
                vector<double> rightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1].computeConservedVariableVector(materialParameters);
                
                coarseCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5,
                                                                                               VectorAlgebra::addVectors(leftFineConservedVariableVector,
                                                                                                                         rightFineConservedVariableVector)),
                                                                 materialParameters);
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
        
        vector<EulerStateVector> fineCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells(fineCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            vector<double> leftConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 1].computeConservedVariableVector(materialParameters);
            vector<double> conservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2].computeConservedVariableVector(materialParameters);
            vector<double> rightConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 3].computeConservedVariableVector(materialParameters);
            
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
                
                vector<double> leftFineConservedVariableVector = fineCurrentCells[(i * 2)].computeConservedVariableVector(materialParameters);
                vector<double> rightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1].computeConservedVariableVector(materialParameters);
                
                coarseCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5, VectorAlgebra::addVectors(leftFineConservedVariableVector,
                                                                                                                              rightFineConservedVariableVector)),
                                                                 materialParameters);
            }
        }
    }
    
    return tuple<vector<EulerStateVector>, vector<EulerStateVector>, vector<bool>, vector<bool>>(coarseCurrentCells, fineCurrentCells, coarseAMRStructure,
                                                                                                 fineAMRStructure);
}

tuple<vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<bool> >, vector<vector<bool> >>
EulerAMR::solve2DLevel1AMR(vector<vector<EulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance,
                           int order, EulerMaterialParameters materialParameters)
{
    long rowCount = initialCells.size();
    long columnCount = initialCells[0].size();
    
    vector<vector<EulerStateVector> > coarseCurrentCells(rowCount, vector<EulerStateVector>(columnCount));
    vector<vector<EulerStateVector> > fineCurrentCells(rowCount * 2, vector<EulerStateVector>(columnCount * 2));
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
        double coarseTimeStep = EulerSolvers::computeStableTimeStep2D(coarseCurrentCells, cellSpacing, CFLCoefficient, coarseCurrentTime, finalTime,
                                                                      fineTotalIteration, materialParameters);
        
        vector<vector<EulerStateVector> > coarseCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(coarseCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                vector<double> topConservedVariableVector = coarseCurrentCellsWithBoundary[i][j + 1].computeConservedVariableVector(materialParameters);
                vector<double> leftConservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j].computeConservedVariableVector(materialParameters);
                vector<double> conservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j + 1].computeConservedVariableVector(materialParameters);
                vector<double> rightConservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j + 2].computeConservedVariableVector(materialParameters);
                vector<double> bottomConservedVariableVector = coarseCurrentCellsWithBoundary[i + 2][j + 1].computeConservedVariableVector(materialParameters);
                
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
            computeXFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, coarseAMRStructure, materialParameters);
            computeYFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, coarseAMRStructure, materialParameters);
            computeXFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, coarseAMRStructure, materialParameters);
        }
        else
        {
            computeXSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, 0.0, 1, coarseAMRStructure, materialParameters);
            computeYSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, 0.0, 1, coarseAMRStructure, materialParameters);
            computeXSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, 0.0, 1, coarseAMRStructure, materialParameters);
        }
        
        AMRHelper::outputCoarseStatus(coarseCurrentIteration + 1, coarseCurrentTime + coarseTimeStep, coarseTimeStep);
        
        bool refined = false;
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
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
                double fineTimeStep = EulerSolvers::computeStableTimeStep2D(fineCurrentCells, cellSpacing * 0.5, CFLCoefficient, fineCurrentTime,
                                                                            coarseCurrentTime + coarseTimeStep, fineTotalIteration, materialParameters);
                
                if (order == 1)
                {
                    computeXFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep * 0.5, fineAMRStructure, materialParameters);
                    computeYFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep, fineAMRStructure, materialParameters);
                    computeXFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep * 0.5, fineAMRStructure, materialParameters);
                }
                else
                {
                    computeXSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep * 0.5, 0.0, 1, fineAMRStructure, materialParameters);
                    computeYSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep, 0.0, 1, fineAMRStructure, materialParameters);
                    computeXSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.5, fineTimeStep * 0.5, 0.0, 1, fineAMRStructure, materialParameters);
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
        
#pragma omp paralel for
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                if (!coarseAMRStructure[i][j])
                {
                    vector<double> topLeftFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2)].computeConservedVariableVector(materialParameters);
                    vector<double> topRightFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2) + 1].computeConservedVariableVector(materialParameters);
                    vector<double> bottomLeftFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2)].computeConservedVariableVector(materialParameters);
                    vector<double> bottomRightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2) + 1].computeConservedVariableVector(materialParameters);
                    
                    coarseCurrentCells[i][j].setConservedVariableVector(VectorAlgebra::multiplyVector(0.25, 
                        VectorAlgebra::addVectors(VectorAlgebra::addVectors(VectorAlgebra::addVectors(topLeftFineConservedVariableVector,
                                                                                                      topRightFineConservedVariableVector),
                                                                            bottomLeftFineConservedVariableVector), bottomRightFineConservedVariableVector)),
                                                                        materialParameters);
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
        
        vector<vector<EulerStateVector> > fineCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(fineCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j  < columnCount; j++)
            {
                vector<double> topConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 1][(j * 2) + 2].computeConservedVariableVector(materialParameters);
                vector<double> leftConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 1].computeConservedVariableVector(materialParameters);
                vector<double> conservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 2].computeConservedVariableVector(materialParameters);
                vector<double> rightConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 3].computeConservedVariableVector(materialParameters);
                vector<double> bottomConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 3][(j * 2) + 2].computeConservedVariableVector(materialParameters);
                
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
                    
                    vector<double> topLeftFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2)].computeConservedVariableVector(materialParameters);
                    vector<double> topRightFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2) + 1].computeConservedVariableVector(materialParameters);
                    vector<double> bottomLeftFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2)].computeConservedVariableVector(materialParameters);
                    vector<double> bottomRightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2) + 1].computeConservedVariableVector(materialParameters);
                    
                    coarseCurrentCells[i][j].setConservedVariableVector(VectorAlgebra::multiplyVector(0.25, 
                        VectorAlgebra::addVectors(VectorAlgebra::addVectors(VectorAlgebra::addVectors(topLeftFineConservedVariableVector,
                                                                                                      topRightFineConservedVariableVector),
                                                                            bottomLeftFineConservedVariableVector), bottomRightFineConservedVariableVector)),
                                                                        materialParameters);
                }
            }
        }
    }
    
    return tuple<vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<bool> >,
    vector<vector<bool> >>(coarseCurrentCells, fineCurrentCells, coarseAMRStructure, fineAMRStructure);
}

tuple<vector<EulerStateVector>, vector<EulerStateVector>, vector<EulerStateVector>, vector<bool>, vector<bool>, vector<bool>>
EulerAMR::solveLevel2AMR(vector<EulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance, int order,
                         EulerMaterialParameters materialParameters)
{
    long cellCount = initialCells.size();
    
    vector<EulerStateVector> coarseCurrentCells(cellCount);
    vector<EulerStateVector> intermediateCurrentCells(cellCount * 2);
    vector<EulerStateVector> fineCurrentCells(cellCount * 4);
    
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
        double coarseTimeStep = EulerSolvers::computeStableTimeStep(coarseCurrentCells, cellSpacing, CFLCoefficient, coarseCurrentTime, finalTime, fineTotalIteration,
                                                                    materialParameters);
        
        vector<EulerStateVector> coarseCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells(coarseCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            vector<double> leftConservedVariableVector = coarseCurrentCellsWithBoundary[i].computeConservedVariableVector(materialParameters);
            vector<double> conservedVariableVector = coarseCurrentCellsWithBoundary[i + 1].computeConservedVariableVector(materialParameters);
            vector<double> rightConservedVariableVector = coarseCurrentCellsWithBoundary[i + 2].computeConservedVariableVector(materialParameters);
            
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
            computeFORCETimeStepAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, coarseAMRStructure, materialParameters);
        }
        else
        {
            computeSLICTimeStepAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, 0.0, 1, coarseAMRStructure, materialParameters);
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
                double intermediateTimeStep = EulerSolvers::computeStableTimeStep(intermediateCurrentCells, cellSpacing * 0.5, CFLCoefficient,
                                                                                  intermediateCurrentTime, coarseCurrentTime + coarseTimeStep, fineTotalIteration,
                                                                                  materialParameters);
                
                vector<EulerStateVector> intermediateCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells(intermediateCurrentCells, 1);
                
#pragma omp parallel for
                for (int i = 0; i < (cellCount * 2); i++)
                {
                    vector<double> leftConservedVariableVector = intermediateCurrentCellsWithBoundary[i].computeConservedVariableVector(materialParameters);
                    vector<double> conservedVariableVector = intermediateCurrentCellsWithBoundary[i + 1].computeConservedVariableVector(materialParameters);
                    vector<double> rightConservedVariableVector = intermediateCurrentCellsWithBoundary[i + 2].computeConservedVariableVector(materialParameters);
                    
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
                    computeFORCETimeStepAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep, intermediateAMRStructure, materialParameters);
                }
                else
                {
                    computeSLICTimeStepAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep, 0.0, 1, intermediateAMRStructure, materialParameters);
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
                        double fineTimeStep = EulerSolvers::computeStableTimeStep(fineCurrentCells, cellSpacing * 0.25, CFLCoefficient, fineCurrentTime,
                                                                                  intermediateCurrentTime + intermediateTimeStep, fineTotalIteration, materialParameters);
                        
                        if (order == 1)
                        {
                            computeFORCETimeStepAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep, fineAMRStructure, materialParameters);
                        }
                        else
                        {
                            computeSLICTimeStepAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep, 0.0, 1, fineAMRStructure, materialParameters);
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
                        vector<double> leftFineConservedVariableVector = fineCurrentCells[(i * 2)].computeConservedVariableVector(materialParameters);
                        vector<double> rightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1].computeConservedVariableVector(materialParameters);
                            
                        intermediateCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5,
                                                                                                             VectorAlgebra::addVectors(leftFineConservedVariableVector,
                                                                                                                                       rightFineConservedVariableVector)),
                                                                               materialParameters);
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
                
                vector<EulerStateVector> fineCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells(fineCurrentCells, 1);
                
#pragma omp parallel for
                for (int i = 0; i < (cellCount * 2); i++)
                {
                    vector<double> leftConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 1].computeConservedVariableVector(materialParameters);
                    vector<double> conservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2].computeConservedVariableVector(materialParameters);
                    vector<double> rightConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 3].computeConservedVariableVector(materialParameters);
                    
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
                        
                        vector<double> leftFineConservedVariableVector = fineCurrentCells[(i * 2)].computeConservedVariableVector(materialParameters);
                        vector<double> rightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1].computeConservedVariableVector(materialParameters);
                        
                        intermediateCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5,
                                                                                                             VectorAlgebra::addVectors(leftFineConservedVariableVector,
                                                                                                                                       rightFineConservedVariableVector)),
                                                                               materialParameters);
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
                vector<double> leftIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2)].computeConservedVariableVector(materialParameters);
                vector<double> rightIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2) + 1].computeConservedVariableVector(materialParameters);
                    
                coarseCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5,
                                                                                               VectorAlgebra::addVectors(leftIntermediateConservedVariableVector,
                                                                                                                         rightIntermediateConservedVariableVector)),
                                                                 materialParameters);
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
        
        vector<EulerStateVector> intermediateCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells(intermediateCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < cellCount; i++)
        {
            vector<double> leftConservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 1].computeConservedVariableVector(materialParameters);
            vector<double> conservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 2].computeConservedVariableVector(materialParameters);
            vector<double> rightConservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 3].computeConservedVariableVector(materialParameters);
            
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
                
                vector<double> leftIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2)].computeConservedVariableVector(materialParameters);
                vector<double> rightIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2) + 1].computeConservedVariableVector(materialParameters);
                
                coarseCurrentCells[i].setConservedVariableVector(VectorAlgebra::multiplyVector(0.5,
                                                                                               VectorAlgebra::addVectors(leftIntermediateConservedVariableVector,
                                                                                                                         rightIntermediateConservedVariableVector)),
                                                                 materialParameters);
            }
        }
    }
    
    return tuple<vector<EulerStateVector>, vector<EulerStateVector>, vector<EulerStateVector>, vector<bool>, vector<bool>,
    vector<bool>>(coarseCurrentCells, intermediateCurrentCells, fineCurrentCells, coarseAMRStructure, intermediateAMRStructure, fineAMRStructure);
}

tuple<vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<bool> >, vector<vector<bool> >,
vector<vector<bool> >> EulerAMR::solve2DLevel2AMR(vector<vector<EulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                  double AMRTolerance, int order, EulerMaterialParameters materialParameters)
{
    long rowCount = initialCells.size();
    long columnCount = initialCells[0].size();
    
    vector<vector<EulerStateVector> > coarseCurrentCells(rowCount, vector<EulerStateVector>(columnCount));
    vector<vector<EulerStateVector> > intermediateCurrentCells(rowCount * 2, vector<EulerStateVector>(columnCount * 2));
    vector<vector<EulerStateVector> > fineCurrentCells(rowCount * 4, vector<EulerStateVector>(columnCount * 4));
    
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
        double coarseTimeStep = EulerSolvers::computeStableTimeStep2D(coarseCurrentCells, cellSpacing, CFLCoefficient, coarseCurrentTime, finalTime,
                                                                      fineTotalIteration, materialParameters);
        
        vector<vector<EulerStateVector> > coarseCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(coarseCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                vector<double> topConservedVariableVector = coarseCurrentCellsWithBoundary[i][j + 1].computeConservedVariableVector(materialParameters);
                vector<double> leftConservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j].computeConservedVariableVector(materialParameters);
                vector<double> conservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j + 1].computeConservedVariableVector(materialParameters);
                vector<double> rightConservedVariableVector = coarseCurrentCellsWithBoundary[i + 1][j + 2].computeConservedVariableVector(materialParameters);
                vector<double> bottomConservedVariableVector = coarseCurrentCellsWithBoundary[i + 2][j + 1].computeConservedVariableVector(materialParameters);
                
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
            computeXFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, coarseAMRStructure, materialParameters);
            computeYFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, coarseAMRStructure, materialParameters);
            computeXFORCETimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, coarseAMRStructure, materialParameters);
        }
        else
        {
            computeXSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, 0.0, 1, coarseAMRStructure, materialParameters);
            computeYSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep, 0.0, 1, coarseAMRStructure, materialParameters);
            computeXSLICTimeStep2DAMR(coarseCurrentCells, cellSpacing, coarseTimeStep * 0.5, 0.0, 1, coarseAMRStructure, materialParameters);
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
                double intermediateTimeStep = EulerSolvers::computeStableTimeStep2D(intermediateCurrentCells, cellSpacing * 0.5, CFLCoefficient,
                                                                                    intermediateCurrentTime, coarseCurrentTime + coarseTimeStep, fineTotalIteration,
                                                                                    materialParameters);
                
                vector<vector<EulerStateVector> > intermediateCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(intermediateCurrentCells, 1);
                
#pragma omp parallel for
                for (int i = 0; i < (rowCount * 2); i++)
                {
                    for (int j = 0; j < (columnCount * 2); j++)
                    {
                        vector<double> topConservedVariableVector = intermediateCurrentCellsWithBoundary[i][j + 1].computeConservedVariableVector(materialParameters);
                        vector<double> leftConservedVariableVector = intermediateCurrentCellsWithBoundary[i + 1][j].computeConservedVariableVector(materialParameters);
                        vector<double> conservedVariableVector = intermediateCurrentCellsWithBoundary[i + 1][j + 1].computeConservedVariableVector(materialParameters);
                        vector<double> rightConservedVariableVector = intermediateCurrentCellsWithBoundary[i + 1][j + 2].
                        computeConservedVariableVector(materialParameters);
                        vector<double> bottomConservedVariableVector = intermediateCurrentCellsWithBoundary[i + 2][j + 1].
                        computeConservedVariableVector(materialParameters);
                        
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
                    computeXFORCETimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep * 0.5, intermediateAMRStructure, materialParameters);
                    computeYFORCETimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep, intermediateAMRStructure, materialParameters);
                    computeXFORCETimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep * 0.5, intermediateAMRStructure, materialParameters);
                }
                else
                {
                    computeXSLICTimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep * 0.5, 0.0, 1, intermediateAMRStructure,
                                              materialParameters);
                    computeYSLICTimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep, 0.0, 1, intermediateAMRStructure,
                                              materialParameters);
                    computeXSLICTimeStep2DAMR(intermediateCurrentCells, cellSpacing * 0.5, intermediateTimeStep * 0.5, 0.0, 1, intermediateAMRStructure,
                                              materialParameters);
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
                        double fineTimeStep = EulerSolvers::computeStableTimeStep2D(fineCurrentCells, cellSpacing * 0.25, CFLCoefficient, fineCurrentTime,
                                                                                    intermediateCurrentTime + intermediateTimeStep, fineTotalIteration,
                                                                                    materialParameters);
                        
                        if (order == 1)
                        {
                            computeXFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep * 0.5, fineAMRStructure, materialParameters);
                            computeYFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep, fineAMRStructure, materialParameters);
                            computeXFORCETimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep * 0.5, fineAMRStructure, materialParameters);
                        }
                        else
                        {
                            computeXSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep * 0.5, 0.0, 1, fineAMRStructure, materialParameters);
                            computeYSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep, 0.0, 1, fineAMRStructure, materialParameters);
                            computeXSLICTimeStep2DAMR(fineCurrentCells, cellSpacing * 0.25, fineTimeStep * 0.5, 0.0, 1, fineAMRStructure, materialParameters);
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
                        if (!intermediateAMRStructure[i][j] && (fineAMRStructure[(i * 2)][(j * 2)] && fineAMRStructure[(i * 2)][(j * 2) + 1] &&
                                                                fineAMRStructure[(i * 2) + 1][(j * 2)] && fineAMRStructure[(i * 2) + 1][(j * 2) + 1]))
                        {
                            vector<double> topLeftFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2)].computeConservedVariableVector(materialParameters);
                            vector<double> topRightFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2) + 1].computeConservedVariableVector(materialParameters);
                            vector<double> bottomLeftFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2)].
                            computeConservedVariableVector(materialParameters);
                            vector<double> bottomRightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2) + 1].
                            computeConservedVariableVector(materialParameters);
                            
                            intermediateCurrentCells[i][j].setConservedVariableVector(VectorAlgebra::multiplyVector(0.25,
                                VectorAlgebra::addVectors(VectorAlgebra::addVectors(VectorAlgebra::addVectors(topLeftFineConservedVariableVector,
                                                                                                              topRightFineConservedVariableVector),
                                                                                    bottomLeftFineConservedVariableVector), bottomRightFineConservedVariableVector)),
                                                                                      materialParameters);
                        }
                        
                        if (!fineAMRStructure[(i * 2)][(j * 2)] && intermediateAMRStructure[i][j])
                        {
                            fineCurrentCells[(i * 2)][(j * 2)] = intermediateCurrentCells[i][j];
                        }
                        if (!fineAMRStructure[(i * 2)][(j * 2) + 1] && intermediateAMRStructure[i][j])
                        {
                            fineCurrentCells[(i * 2)][(j * 2) + 1] = intermediateCurrentCells[i][j];
                        }
                        if (!fineAMRStructure[(i * 2) + 1][(j * 2)] && intermediateAMRStructure[i][j])
                        {
                            fineCurrentCells[(i * 2) + 1][(j * 2)] = intermediateCurrentCells[i][j];
                        }
                        if (!fineAMRStructure[(i * 2) + 1][(j * 2) + 1] && intermediateAMRStructure[i][j])
                        {
                            fineCurrentCells[(i * 2) + 1][(j * 2) + 1] = intermediateCurrentCells[i][j];
                        }
                    }
                }
                
                vector<vector<EulerStateVector> > fineCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(fineCurrentCells, 1);
                
#pragma omp parallel for
                for (int i = 0; i < (rowCount * 2); i++)
                {
                    for (int j = 0; j < (columnCount * 2); j++)
                    {
                        vector<double> topConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 1][(j * 2) + 2].
                        computeConservedVariableVector(materialParameters);
                        vector<double> leftConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 1].
                        computeConservedVariableVector(materialParameters);
                        vector<double> conservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 2].
                        computeConservedVariableVector(materialParameters);
                        vector<double> rightConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 3].
                        computeConservedVariableVector(materialParameters);
                        vector<double> bottomConservedVariableVector = fineCurrentCellsWithBoundary[(i * 2) + 3][(j * 2) + 2].
                        computeConservedVariableVector(materialParameters);
                        
                        long conservedVariableCount = conservedVariableVector.size();
                        bool markedCell = false;
                        
                        for (int k = 0; k < conservedVariableCount; k++)
                        {
                            if (abs((conservedVariableVector[k] - topConservedVariableVector[k]) / cellSpacing) >= AMRTolerance ||
                                abs((bottomConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= AMRTolerance ||
                                abs((conservedVariableVector[k] - leftConservedVariableVector[k]) / cellSpacing) >= AMRTolerance ||
                                abs((rightConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= AMRTolerance)
                            {
                                markedCell = true;
                            }
                        }
                        
                        if (!markedCell && (fineAMRStructure[(i * 2)][(j * 2)] && fineAMRStructure[(i * 2)][(j * 2) + 1] &&
                                            fineAMRStructure[(i * 2) + 1][(j * 2)] && fineAMRStructure[(i * 2) + 1][(j * 2) + 1]))
                        {
                            fineAMRStructure[(i * 2)][(j * 2)] = false;
                            fineAMRStructure[(i * 2)][(j * 2) + 1] = false;
                            fineAMRStructure[(i * 2) + 1][(j * 2)] = false;
                            fineAMRStructure[(i * 2) + 1][(j * 2) + 1] = false;
                            
                            intermediateAMRStructure[i][j] = true;
                            
                            vector<double> topLeftFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2)].computeConservedVariableVector(materialParameters);
                            vector<double> topRightFineConservedVariableVector = fineCurrentCells[(i * 2)][(j * 2) + 1].
                            computeConservedVariableVector(materialParameters);
                            vector<double> bottomLeftFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2)].
                            computeConservedVariableVector(materialParameters);
                            vector<double> bottomRightFineConservedVariableVector = fineCurrentCells[(i * 2) + 1][(j * 2) + 1].
                            computeConservedVariableVector(materialParameters);
                            
                            intermediateCurrentCells[i][j].setConservedVariableVector(VectorAlgebra::multiplyVector(0.25,
                                VectorAlgebra::addVectors(VectorAlgebra::addVectors(VectorAlgebra::addVectors(topLeftFineConservedVariableVector,
                                                                                                              topRightFineConservedVariableVector),
                                                                                    bottomLeftFineConservedVariableVector), bottomRightFineConservedVariableVector)),
                                                                                      materialParameters);
                        }
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
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                if (!coarseAMRStructure[i][j])
                {
                    vector<double> topLeftIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2)][(j * 2)].
                    computeConservedVariableVector(materialParameters);
                    vector<double> topRightIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2)][(j * 2) + 1].
                    computeConservedVariableVector(materialParameters);
                    vector<double> bottomLeftIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2) + 1][(j * 2)].
                    computeConservedVariableVector(materialParameters);
                    vector<double> bottomRightIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2) + 1][(j * 2) + 1].
                    computeConservedVariableVector(materialParameters);
                    
                    coarseCurrentCells[i][j].setConservedVariableVector(VectorAlgebra::multiplyVector(0.25,
                        VectorAlgebra::addVectors(VectorAlgebra::addVectors(VectorAlgebra::addVectors(topLeftIntermediateConservedVariableVector,
                                                                                                      topRightIntermediateConservedVariableVector),
                                                                            bottomLeftIntermediateConservedVariableVector),
                                                  bottomRightIntermediateConservedVariableVector)), materialParameters);
                }
                
                if (!intermediateAMRStructure[(i * 2)][(j * 2)])
                {
                    intermediateCurrentCells[(i * 2)][(j * 2)] = coarseCurrentCells[i][j];
                }
                if (!intermediateAMRStructure[(i * 2)][(j * 2) + 1])
                {
                    intermediateCurrentCells[(i * 2)][(j * 2) + 1] = coarseCurrentCells[i][j];
                }
                if (!intermediateAMRStructure[(i * 2) + 1][(j * 2)])
                {
                    intermediateCurrentCells[(i * 2) + 1][(j * 2)] = coarseCurrentCells[i][j];
                }
                if (!intermediateAMRStructure[(i * 2) + 1][(j * 2) + 1])
                {
                    intermediateCurrentCells[(i * 2) + 1][(j * 2) + 1] = coarseCurrentCells[i][j];
                }
            }
        }
        
        vector<vector<EulerStateVector> > intermediateCurrentCellsWithBoundary = EulerSolvers::insertBoundaryCells2D(intermediateCurrentCells, 1);
        
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                vector<double> topConservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 1][(j * 2) + 2].
                computeConservedVariableVector(materialParameters);
                vector<double> leftConservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 1].
                computeConservedVariableVector(materialParameters);
                vector<double> conservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 2].
                computeConservedVariableVector(materialParameters);
                vector<double> rightConservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 2][(j * 2) + 3].
                computeConservedVariableVector(materialParameters);
                vector<double> bottomConservedVariableVector = intermediateCurrentCellsWithBoundary[(i * 2) + 3][(j * 2) + 2].
                computeConservedVariableVector(materialParameters);
                
                long conservedVariableCount = conservedVariableVector.size();
                bool markedCell = false;
                
                for (int k = 0; k < conservedVariableCount; k++)
                {
                    if (abs((conservedVariableVector[k] - topConservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 0.5 * (1.0 / order)) ||
                        abs((bottomConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 0.5 * (1.0 / order)) ||
                        abs((conservedVariableVector[k] - leftConservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 0.5 * (1.0 / order)) ||
                        abs((rightConservedVariableVector[k] - conservedVariableVector[k]) / cellSpacing) >= (AMRTolerance * 0.5 * (1.0 / order)))
                    {
                        markedCell = true;
                    }
                }
                
                if (!markedCell && (intermediateAMRStructure[(i * 2)][(j * 2)] && intermediateAMRStructure[(i * 2)][(j * 2) + 1] &&
                                    intermediateAMRStructure[(i * 2) + 1][(j * 2)] && intermediateAMRStructure[(i * 2) + 1][(j * 2) + 1]))
                {
                    intermediateAMRStructure[(i * 2)][(j * 2)] = false;
                    intermediateAMRStructure[(i * 2)][(j * 2) + 1] = false;
                    intermediateAMRStructure[(i * 2) + 1][(j * 2)] = false;
                    intermediateAMRStructure[(i * 2) + 1][(j * 2) + 1] = false;
                    
                    coarseAMRStructure[i][j] = true;
                    
                    vector<double> topLeftIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2)][(j * 2)].
                    computeConservedVariableVector(materialParameters);
                    vector<double> topRightIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2)][(j * 2) + 1].
                    computeConservedVariableVector(materialParameters);
                    vector<double> bottomLeftIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2) + 1][(j * 2)].
                    computeConservedVariableVector(materialParameters);
                    vector<double> bottomRightIntermediateConservedVariableVector = intermediateCurrentCells[(i * 2) + 1][(j * 2) + 1].
                    computeConservedVariableVector(materialParameters);
                    
                    coarseCurrentCells[i][j].setConservedVariableVector(VectorAlgebra::multiplyVector(0.25,
                        VectorAlgebra::addVectors(VectorAlgebra::addVectors(VectorAlgebra::addVectors(topLeftIntermediateConservedVariableVector,
                                                                                                      topRightIntermediateConservedVariableVector),
                                                                            bottomLeftIntermediateConservedVariableVector),
                                                  bottomRightIntermediateConservedVariableVector)), materialParameters);
                }
            }
        }
    }
    
    return tuple<vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<bool> >,
    vector<vector<bool> >, vector<vector<bool> >>(coarseCurrentCells, intermediateCurrentCells, fineCurrentCells, coarseAMRStructure, intermediateAMRStructure,
                                                  fineAMRStructure);
}
