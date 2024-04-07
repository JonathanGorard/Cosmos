#ifndef BlackHoleEulerSolvers_hpp
#define BlackHoleEulerSolvers_hpp

#include "../BlackHoleEuler/BlackHoleEulerStateVector.hpp"
#include "EulerSolvers.hpp"

class BlackHoleEulerSolvers
{
public:
    static vector<BlackHoleEulerStateVector> insertBoundaryCells(vector<BlackHoleEulerStateVector> & currentCells, int boundarySize);
    static vector<vector<BlackHoleEulerStateVector> > insertBoundaryCells2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, int boundarySize);
    
    static double computeMaximumWaveSpeed(vector<BlackHoleEulerStateVector> & currentCells, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static double computeMaximumWaveSpeed2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, EulerMaterialParameters materialParameters,
                                            BlackHoleSpacetime blackHole);
    
    static double computeStableTimeStep(vector<BlackHoleEulerStateVector> currentCells, double cellSpacing, double CFLCoefficient, double currentTime,
                                        double finalTime, int currentIteration, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static double computeStableTimeStep2D(vector<vector<BlackHoleEulerStateVector> > currentCells, double cellSpacing, double CFLCoefficient, double currentTime,
                                          double finalTime, int currentIteration, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static BlackHoleEulerStateVector evolveStateByFractionalTimeStep(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution,
                                                                     double xCoordinate, double yCoordinate, double zCoordinate,
                                                                     EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static BlackHoleEulerStateVector evolveStateByFractionalXTimeStep(double stepFraction, BlackHoleEulerStateVector leftStateVector,
                                                                      BlackHoleEulerStateVector middleStateVector, BlackHoleEulerStateVector rightStateVector,
                                                                      double cellSpacing, double timeStep, double bias, int slopeLimiter, double xCoordinate,
                                                                      double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                                      BlackHoleSpacetime blackHole);
    static BlackHoleEulerStateVector evolveStateByFractionalYTimeStep(double stepFraction, BlackHoleEulerStateVector topStateVector,
                                                                      BlackHoleEulerStateVector middleStateVector, BlackHoleEulerStateVector bottomStateVector,
                                                                      double cellSpacing, double timeStep, double bias, int slopeLimiter, double xCoordinate,
                                                                      double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                                      BlackHoleSpacetime blackHole);
    
    static BlackHoleEulerStateVector evolveStateByHalfXTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue,
                                                                vector<double> evolutionVector, int side, double cellSpacing, double xCoordinate, double yCoordinate,
                                                                double zCoordinate, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static BlackHoleEulerStateVector evolveStateByHalfXTimeStep(BlackHoleEulerStateVector leftStateVector, BlackHoleEulerStateVector middleStateVector,
                                                                BlackHoleEulerStateVector rightStateVector, double cellSpacing, double timeStep, double bias,
                                                                int slopeLimiter, int side, double xCoordinate, double yCoordinate, double zCoordinate,
                                                                EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static BlackHoleEulerStateVector evolveStateByHalfYTimeStep(vector<double> topExtrapolatedValue, vector<double> bottomExtrapolatedValue,
                                                                vector<double> evolutionVector, int side, double cellSpacing, double xCoordinate, double yCoordinate,
                                                                double zCoordinate, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static BlackHoleEulerStateVector evolveStateByHalfYTimeStep(BlackHoleEulerStateVector topStateVector, BlackHoleEulerStateVector middleStateVector,
                                                                BlackHoleEulerStateVector bottomStateVector, double cellSpacing, double timeStep, double bias,
                                                                int slopeLimiter, int side, double xCoordinate, double yCoordinate, double zCoordinate,
                                                                EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
};

#endif
