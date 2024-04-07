#ifndef BlackHoleEulerSecondOrderSolver_hpp
#define BlackHoleEulerSecondOrderSolver_hpp

#include "BlackHoleEulerFirstOrderSolver.hpp"

class BlackHoleEulerSecondOrderSolver
{
public:
    static vector<double> computeXSLICFlux(BlackHoleEulerStateVector leftLeftStateVector, BlackHoleEulerStateVector leftStateVector,
                                           BlackHoleEulerStateVector rightStateVector, BlackHoleEulerStateVector rightRightStateVector, double cellSpacing,
                                           double timeStep, double bias, int slopeLimiter, double xCoordinate, double yCoordinate, double zCoordinate,
                                           EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static vector<double> computeYSLICFlux(BlackHoleEulerStateVector topTopStateVector, BlackHoleEulerStateVector topStateVector,
                                           BlackHoleEulerStateVector bottomStateVector, BlackHoleEulerStateVector bottomBottomStateVector, double cellSpacing,
                                           double timeStep, double bias, int slopeLimiter, double xCoordinate, double yCoordinate, double zCoordinate,
                                           EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static void computeSLICTimeStep(vector<BlackHoleEulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                    EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static void computeXSLICTimeStep2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                       int slopeLimiter, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static void computeYSLICTimeStep2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                       int slopeLimiter, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static vector<BlackHoleEulerStateVector> solve(vector<BlackHoleEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                   double bias, int slopeLimiter, int subcyclingIterations, EulerMaterialParameters materialParameters,
                                                   BlackHoleSpacetime blackHole);
    
    static vector<vector<BlackHoleEulerStateVector> > solve2D(vector<vector<BlackHoleEulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient,
                                                              double finalTime, double bias, int slopeLimiter, int subcyclingIterations,
                                                              EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
};

#endif
