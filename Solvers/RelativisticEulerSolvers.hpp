#ifndef RelativisticEulerSolvers_hpp
#define RelativisticEulerSolvers_hpp

#include "../RelativisticEuler/RelativisticEulerStateVector.hpp"
#include "EulerSolvers.hpp"

class RelativisticEulerSolvers
{
public:
    static vector<RelativisticEulerStateVector> insertBoundaryCells(vector<RelativisticEulerStateVector> & currentCells, int boundarySize);
    static vector<vector<RelativisticEulerStateVector> > insertBoundaryCells2D(vector<vector<RelativisticEulerStateVector> > & currentCells, int boundarySize);
    
    static double computeMaximumWaveSpeed(vector<RelativisticEulerStateVector> & currentCells, EulerMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed2D(vector<vector<RelativisticEulerStateVector> > & currentCells, EulerMaterialParameters materialParameters);
    
    static double computeStableTimeStep(vector<RelativisticEulerStateVector> currentCells, double cellSpacing, double CFLCoefficient, double currentTime,
                                        double finalTime, int currentIteration, EulerMaterialParameters materialParameters);
    static double computeStableTimeStep2D(vector<vector<RelativisticEulerStateVector> > currentCells, double cellSpacing, double CFLCoefficient, double currentTime,
                                          double finalTime, int currentIteration, EulerMaterialParameters materialParameters);
    
    static RelativisticEulerStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue,
                                                                  vector<double> evolutionVector, int side, EulerMaterialParameters materialParameters);
    static RelativisticEulerStateVector evolveStateByHalfXTimeStep(RelativisticEulerStateVector leftStateVector, RelativisticEulerStateVector middleStateVector,
                                                                   RelativisticEulerStateVector rightStateVector, double cellSpacing, double timeStep, double bias,
                                                                   int slopeLimiter, int side, EulerMaterialParameters materialParameters);
    static RelativisticEulerStateVector evolveStateByHalfYTimeStep(RelativisticEulerStateVector topStateVector, RelativisticEulerStateVector middleStateVector,
                                                                   RelativisticEulerStateVector bottomStateVector, double cellSpacing, double timeStep, double bias,
                                                                   int slopeLimiter, int side, EulerMaterialParameters materialParameters);
};

#endif
