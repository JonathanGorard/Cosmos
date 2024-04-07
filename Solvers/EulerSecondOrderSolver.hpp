#ifndef EulerSecondOrderSolver_hpp
#define EulerSecondOrderSolver_hpp

#include "EulerFirstOrderSolver.hpp"

class EulerSecondOrderSolver
{
public:
    static vector<double> computeXSLICFlux(EulerStateVector leftLeftStateVector, EulerStateVector leftStateVector, EulerStateVector rightStateVector,
                                           EulerStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           EulerMaterialParameters materialParameters);
    static vector<double> computeYSLICFlux(EulerStateVector topTopStateVector, EulerStateVector topStateVector, EulerStateVector bottomStateVector,
                                           EulerStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           EulerMaterialParameters materialParameters);
    
    static void computeSLICTimeStep(vector<EulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                    EulerMaterialParameters materialParameters);
    
    static void computeXSLICTimeStep2D(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                       EulerMaterialParameters materialParameters);
    static void computeYSLICTimeStep2D(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                       EulerMaterialParameters materialParameters);
    
    static vector<EulerStateVector> solve(vector<EulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                          int slopeLimiter, EulerMaterialParameters materialParameters);
    
    static vector<vector<EulerStateVector> > solve2D(vector<vector<EulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                     double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
};

#endif
