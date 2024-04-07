#ifndef EulerSolvers_hpp
#define EulerSolvers_hpp

#include "../Euler/EulerStateVector.hpp"
#include "SlopeLimiters.hpp"
#include "/usr/local/include/omp.h"
#include <iostream>
using namespace std;

class EulerSolvers
{
public:
    static vector<EulerStateVector> insertBoundaryCells(vector<EulerStateVector> & currentCells, int boundarySize);
    static vector<vector<EulerStateVector> > insertBoundaryCells2D(vector<vector<EulerStateVector> > & currentCells, int boundarySize);
    
    static double computeMaximumWaveSpeed(vector<EulerStateVector> & currentCells, EulerMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed2D(vector<vector<EulerStateVector> > & currentCells, EulerMaterialParameters materialParameters);
    
    static double computeStableTimeStep(double timeStep, double currentTime, double finalTime, int curentIteration);
    static double computeStableTimeStep(vector<EulerStateVector> currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                        int currentIteration, EulerMaterialParameters materialParameters);
    static double computeStableTimeStep2D(vector<vector<EulerStateVector> > currentCells, double cellSpacing, double CFLCoefficient, double currentTime,
                                          double finalTime, int currentIteration, EulerMaterialParameters materialParameters);
    
    static vector<double> computeFractionalEvolutionVector(double stepFraction, vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing,
                                                           double timeStep);
    static vector<double> computeEvolutionVector(vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing, double timeStep);
    
    static EulerStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector,
                                                      int side, EulerMaterialParameters materialParameters);
    static EulerStateVector evolveStateByHalfXTimeStep(EulerStateVector leftStateVector, EulerStateVector middleStateVector, EulerStateVector rightStateVector,
                                                       double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                       EulerMaterialParameters materialParameters);
    static EulerStateVector evolveStateByHalfYTimeStep(EulerStateVector topStateVector, EulerStateVector middleStateVector, EulerStateVector bottomStateVector,
                                                       double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                       EulerMaterialParameters materialParameters);
    
    static void outputStatus(int currentIteration, double currentTime, double timeStep);
};

#endif
