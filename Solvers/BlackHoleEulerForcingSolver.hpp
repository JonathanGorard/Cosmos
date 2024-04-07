#ifndef BlackHoleEulerForcingSolver_hpp
#define BlackHoleEulerForcingSolver_hpp

#include "BlackHoleEulerSolvers.hpp"

class BlackHoleEulerForcingSolver
{
public:
    static vector<double> evolveConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                        vector<double> rightConservedVaraibleVector, double xCoordinate, double yCoordinate, double zCoordinate,
                                                        double cellSpacing, double timeStep, double bias, int slopeLimiter, EulerMaterialParameters materialParameters,
                                                        BlackHoleSpacetime blackHole);
    static vector<double> evolveConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                          vector<double> rightConservedVariableVector, vector<double> topConservedVariableVector,
                                                          vector<double> bottomConservedVariableVector, double xCoordinate, double yCoordinate, double zCoordinate,
                                                          double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                          EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static void computeRungeKuttaTimeStep(vector<BlackHoleEulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                          EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static void computeRungeKuttaTimeStep2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                            int slopeLimiter, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
};

#endif
