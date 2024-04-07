#include "BlackHoleEulerForcingSolver.hpp"

vector<double> BlackHoleEulerForcingSolver::evolveConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                                          vector<double> rightConservedVariableVector, double xCoordinate, double yCoordinate,
                                                                          double zCoordinate, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                          EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    BlackHoleEulerStateVector leftStateVector;
    BlackHoleEulerStateVector middleStateVector;
    BlackHoleEulerStateVector rightStateVector;
    
    leftStateVector.setConservedVariableVector(leftConservedVariableVector, xCoordinate - cellSpacing, yCoordinate, zCoordinate, materialParameters, blackHole);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, xCoordinate + cellSpacing, yCoordinate, zCoordinate, materialParameters, blackHole);
    
    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, BlackHoleEulerStateVector::computeSourceTermVector(middleConservedVariableVector,
                                                                                                                          xCoordinate, yCoordinate, zCoordinate,
                                                                                                                          cellSpacing, materialParameters, blackHole));
    
    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    BlackHoleEulerStateVector middleStateVectorFirstStepEvolved = BlackHoleEulerSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector,
                                                                                                                          rightStateVector, cellSpacing, timeStep,
                                                                                                                          bias, slopeLimiter, xCoordinate,
                                                                                                                          yCoordinate, zCoordinate, materialParameters,
                                                                                                                          blackHole);
    
    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(xCoordinate, yCoordinate,
                                                                                                                                    zCoordinate, materialParameters,
                                                                                                                                    blackHole);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep,
                                                              BlackHoleEulerStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved,
                                                                                                                 xCoordinate, yCoordinate, zCoordinate, cellSpacing,
                                                                                                                 materialParameters, blackHole));
    
    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    BlackHoleEulerStateVector middleStateVectorSecondStepEvolved = BlackHoleEulerSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector,
                                                                                                                           rightStateVector, cellSpacing, timeStep,
                                                                                                                           bias, slopeLimiter, xCoordinate,
                                                                                                                           yCoordinate, zCoordinate, materialParameters,
                                                                                                                           blackHole);
    
    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(xCoordinate, yCoordinate,
                                                                                                                                      zCoordinate, materialParameters,
                                                                                                                                      blackHole);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep,
                                                             BlackHoleEulerStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved,
                                                                                                                xCoordinate, yCoordinate, zCoordinate, cellSpacing,
                                                                                                                materialParameters, blackHole));
    
    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    BlackHoleEulerStateVector middleStateVectorThirdStepEvolved = BlackHoleEulerSolvers::evolveStateByFractionalXTimeStep(1.0, leftStateVector, middleStateVector,
                                                                                                                          rightStateVector, cellSpacing, timeStep,
                                                                                                                          bias, slopeLimiter, xCoordinate,
                                                                                                                          yCoordinate, zCoordinate, materialParameters,
                                                                                                                          blackHole);
    
    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(xCoordinate, yCoordinate,
                                                                                                                                    zCoordinate, materialParameters,
                                                                                                                                    blackHole);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep,
                                                              BlackHoleEulerStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved,
                                                                                                                 xCoordinate, yCoordinate, zCoordinate, cellSpacing,
                                                                                                                 materialParameters, blackHole));
    
    return VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector((1.0 / 6.0),
        VectorAlgebra::addVectors(firstStep, VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(2.0, secondStep),
                                                                       VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(2.0, thirdStep), fourthStep)))));
}

vector<double> BlackHoleEulerForcingSolver::evolveConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                                            vector<double> rightConservedVariableVector, vector<double> topConservedVariableVector,
                                                                            vector<double> bottomConservedVariableVector, double xCoordinate, double yCoordinate,
                                                                            double zCoordinate, double cellSpacing, double timeStep, double bias,
                                                                            int slopeLimiter, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    BlackHoleEulerStateVector leftStateVector;
    BlackHoleEulerStateVector middleStateVector;
    BlackHoleEulerStateVector rightStateVector;
    
    BlackHoleEulerStateVector topStateVector;
    BlackHoleEulerStateVector bottomStateVector;
    
    leftStateVector.setConservedVariableVector(leftConservedVariableVector, xCoordinate - cellSpacing, yCoordinate, zCoordinate, materialParameters, blackHole);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, xCoordinate + cellSpacing, yCoordinate, zCoordinate, materialParameters, blackHole);
    
    topStateVector.setConservedVariableVector(topConservedVariableVector, xCoordinate, yCoordinate - cellSpacing, zCoordinate, materialParameters, blackHole);
    bottomStateVector.setConservedVariableVector(bottomConservedVariableVector, xCoordinate, yCoordinate + cellSpacing, zCoordinate, materialParameters, blackHole);
    
    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, BlackHoleEulerStateVector::computeSourceTermVector(middleConservedVariableVector,
                                                                                                                          xCoordinate, yCoordinate, zCoordinate,
                                                                                                                          cellSpacing, materialParameters, blackHole));
    
    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    
    BlackHoleEulerStateVector middleStateVectorFirstStepEvolved1 = BlackHoleEulerSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector,
                                                                                                                           rightStateVector, cellSpacing, timeStep,
                                                                                                                           bias, slopeLimiter, xCoordinate,
                                                                                                                           yCoordinate, zCoordinate,
                                                                                                                           materialParameters, blackHole);
    BlackHoleEulerStateVector middleStateVectorFirstStepEvolved2 = BlackHoleEulerSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector,
                                                                                                                           middleStateVectorFirstStepEvolved1,
                                                                                                                           bottomStateVector, cellSpacing, timeStep,
                                                                                                                           bias, slopeLimiter, xCoordinate,
                                                                                                                           yCoordinate, zCoordinate,
                                                                                                                           materialParameters, blackHole);
    BlackHoleEulerStateVector middleStateVectorFirstStepEvolved = BlackHoleEulerSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector,
                                                                                                                          middleStateVectorFirstStepEvolved2,
                                                                                                                          rightStateVector, cellSpacing, timeStep,
                                                                                                                          bias, slopeLimiter, xCoordinate,
                                                                                                                          yCoordinate, zCoordinate,
                                                                                                                          materialParameters, blackHole);
    
    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(xCoordinate, yCoordinate,
                                                                                                                                    zCoordinate, materialParameters,
                                                                                                                                    blackHole);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep,
                                                              BlackHoleEulerStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved,
                                                                                                                 xCoordinate, yCoordinate, zCoordinate, cellSpacing,
                                                                                                                 materialParameters, blackHole));
    
    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    
    BlackHoleEulerStateVector middleStateVectorSecondStepEvolved1 = BlackHoleEulerSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector,
                                                                                                                            rightStateVector, cellSpacing, timeStep,
                                                                                                                            bias, slopeLimiter, xCoordinate,
                                                                                                                            yCoordinate, zCoordinate,
                                                                                                                            materialParameters, blackHole);
    BlackHoleEulerStateVector middleStateVectorSecondStepEvolved2 = BlackHoleEulerSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector,
                                                                                                                            middleStateVectorSecondStepEvolved1,
                                                                                                                            bottomStateVector, cellSpacing, timeStep,
                                                                                                                            bias, slopeLimiter, xCoordinate,
                                                                                                                            yCoordinate, zCoordinate,
                                                                                                                            materialParameters, blackHole);
    BlackHoleEulerStateVector middleStateVectorSecondStepEvolved = BlackHoleEulerSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector,
                                                                                                                           middleStateVectorSecondStepEvolved2,
                                                                                                                           rightStateVector, cellSpacing, timeStep,
                                                                                                                           bias, slopeLimiter, xCoordinate,
                                                                                                                           yCoordinate, zCoordinate,
                                                                                                                           materialParameters, blackHole);
    
    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(xCoordinate, yCoordinate,
                                                                                                                                      zCoordinate, materialParameters,
                                                                                                                                      blackHole);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep,
                                                             BlackHoleEulerStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved,
                                                                                                                xCoordinate, yCoordinate, zCoordinate, cellSpacing,
                                                                                                                materialParameters, blackHole));
    
    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    
    BlackHoleEulerStateVector middleStateVectorThirdStepEvolved1 = BlackHoleEulerSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector,
                                                                                                                           rightStateVector, cellSpacing, timeStep,
                                                                                                                           bias, slopeLimiter, xCoordinate,
                                                                                                                           yCoordinate, zCoordinate,
                                                                                                                           materialParameters, blackHole);
    BlackHoleEulerStateVector middleStateVectorThirdStepEvolved2 = BlackHoleEulerSolvers::evolveStateByFractionalYTimeStep(1.0, topStateVector,
                                                                                                                           middleStateVectorThirdStepEvolved1,
                                                                                                                           bottomStateVector, cellSpacing, timeStep,
                                                                                                                           bias, slopeLimiter, xCoordinate,
                                                                                                                           yCoordinate, zCoordinate,
                                                                                                                           materialParameters, blackHole);
    BlackHoleEulerStateVector middleStateVectorThirdStepEvolved = BlackHoleEulerSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector,
                                                                                                                          middleStateVectorThirdStepEvolved2,
                                                                                                                          rightStateVector, cellSpacing, timeStep,
                                                                                                                          bias, slopeLimiter, xCoordinate,
                                                                                                                          yCoordinate, zCoordinate,
                                                                                                                          materialParameters, blackHole);
    
    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(xCoordinate, yCoordinate,
                                                                                                                                    zCoordinate, materialParameters,
                                                                                                                                    blackHole);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep,
                                                              BlackHoleEulerStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved,
                                                                                                                 xCoordinate, yCoordinate, zCoordinate, cellSpacing,
                                                                                                                 materialParameters, blackHole));
    
    return VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector((1.0 / 6.0),
        VectorAlgebra::addVectors(firstStep, VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(2.0, secondStep),
                                                                       VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(2.0, thirdStep), fourthStep)))));
}

void BlackHoleEulerForcingSolver::computeRungeKuttaTimeStep(vector<BlackHoleEulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias,
                                                            int slopeLimiter, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long cellCount = currentCells.size();
    vector<BlackHoleEulerStateVector> currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        if (!blackHole.inExcisionRegion(i * cellSpacing, 0.0, 0.0) && !blackHole.inExcisionRegion((i - 1) * cellSpacing, 0.0, 0.0) &&
            !blackHole.inExcisionRegion((i + 1) * cellSpacing, 0.0, 0.0))
        {
            vector<double> leftConservedVariableVector = currentCellsWithBoundary[i].computeConservedVariableVector((i - 1) * cellSpacing, 0.0, 0.0,
                                                                                                                    materialParameters, blackHole);
            vector<double> middleConservedVariableVector = currentCellsWithBoundary[i + 1].computeConservedVariableVector(i * cellSpacing, 0.0, 0.0,
                                                                                                                          materialParameters, blackHole);
            vector<double> rightConservedVariableVector = currentCellsWithBoundary[i + 2].computeConservedVariableVector((i + 1) * cellSpacing, 0.0, 0.0,
                                                                                                                         materialParameters, blackHole);
            
            currentCells[i].setConservedVariableVector(evolveConservedVariableVector(leftConservedVariableVector, middleConservedVariableVector,
                                                                                     rightConservedVariableVector, i * cellSpacing, 0.0, 0.0, cellSpacing, timeStep,
                                                                                     bias, slopeLimiter, materialParameters, blackHole), i * cellSpacing, 0.0, 0.0,
                                                       materialParameters, blackHole);
        }
    }
}

void BlackHoleEulerForcingSolver::computeRungeKuttaTimeStep2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                                              double bias, int slopeLimiter, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    long rowCount = currentCells.size();
    long columnCount = currentCells[0].size();
    vector<vector<BlackHoleEulerStateVector> > currentCellsWithBoundary = BlackHoleEulerSolvers::insertBoundaryCells2D(currentCells, 1);
    
#pragma omp parallel for
    for (int i = 0 ; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (!blackHole.inExcisionRegion(j * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion((j - 1) * cellSpacing, i * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion((j + 1) * cellSpacing, i * cellSpacing, 0.0) && !blackHole.inExcisionRegion(j * cellSpacing, (i - 1) * cellSpacing, 0.0) &&
                !blackHole.inExcisionRegion(j * cellSpacing, (i + 1) * cellSpacing, 0.0))
            {
                vector<double> leftConservedVariableVector = currentCellsWithBoundary[i + 1][j].computeConservedVariableVector((j - 1) * cellSpacing, i * cellSpacing,
                                                                                                                               0.0, materialParameters, blackHole);
                vector<double> middleConservedVariableVector = currentCellsWithBoundary[i + 1][j + 1].computeConservedVariableVector(j * cellSpacing, i * cellSpacing,
                                                                                                                                     0.0, materialParameters, blackHole);
                vector<double> rightConservedVariableVector = currentCellsWithBoundary[i + 1][j + 2].computeConservedVariableVector((j + 1) * cellSpacing, i * cellSpacing,
                                                                                                                                    0.0, materialParameters, blackHole);
                
                vector<double> topConservedVariableVector = currentCellsWithBoundary[i][j + 1].computeConservedVariableVector(j * cellSpacing, (i - 1) * cellSpacing,
                                                                                                                              0.0, materialParameters, blackHole);
                vector<double> bottomConservedVariableVector = currentCellsWithBoundary[i + 2][j + 1].computeConservedVariableVector(j * cellSpacing, (i + 1) * cellSpacing,
                                                                                                                                     0.0, materialParameters, blackHole);
                
                currentCells[i][j].setConservedVariableVector(evolveConservedVariableVector2D(leftConservedVariableVector, middleConservedVariableVector,
                                                                                              rightConservedVariableVector, topConservedVariableVector,
                                                                                              bottomConservedVariableVector, j * cellSpacing, i * cellSpacing, 0.0,
                                                                                              cellSpacing, timeStep, bias, slopeLimiter, materialParameters, blackHole),
                                                              j * cellSpacing, i * cellSpacing, 0.0, materialParameters, blackHole);
            }
        }
    }
}
