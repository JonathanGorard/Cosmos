#ifndef RelativisticEulerSecondOrderSolver_hpp
#define RelativisticEulerSecondOrderSolver_hpp

#include "RelativisticEulerFirstOrderSolver.hpp"

class RelativisticEulerSecondOrderSolver
{
public:
    static vector<double> computeXSLICFlux(RelativisticEulerStateVector leftLeftStateVector, RelativisticEulerStateVector leftStateVector,
                                           RelativisticEulerStateVector rightStateVector, RelativisticEulerStateVector rightRightStateVector, double cellSpacing,
                                           double timeStep, double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
    static vector<double> computeYSLICFlux(RelativisticEulerStateVector topTopStateVector, RelativisticEulerStateVector topStateVector,
                                           RelativisticEulerStateVector bottomStateVector, RelativisticEulerStateVector bottomBottomStateVector, double cellSpacing,
                                           double timeStep, double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
    
    static void computeSLICTimeStep(vector<RelativisticEulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                    EulerMaterialParameters materialParameters);
    
    static void computeXSLICTimeStep2D(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                       int slopeLimiter, EulerMaterialParameters materialParameters);
    static void computeYSLICTimeStep2D(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                       int slopeLimiter, EulerMaterialParameters materialParameters);
    
    static vector<RelativisticEulerStateVector> solve(vector<RelativisticEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                      double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
    
    static vector<vector<RelativisticEulerStateVector> > solve2D(vector<vector<RelativisticEulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient,
                                                                 double finalTime, double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
};

#endif
