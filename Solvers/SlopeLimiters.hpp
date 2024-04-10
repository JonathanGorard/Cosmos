#ifndef SlopeLimiters_hpp
#define SlopeLimiters_hpp

#include "../Euler/EulerStateVector.hpp"
#include "../RelativisticEuler/RelativisticEulerStateVector.hpp"
#include "../BlackHoleEuler/BlackHoleEulerStateVector.hpp"
#include "../Mathematics/VectorAlgebra.hpp"
#include <omp.h>

class SlopeLimiters
{
public:
    
    static double computeGradientRatio(double steepness, double bias);
    
    static double computeSuperBeeSlopeLimiter(double steepness, double bias);
    static double computeVanLeerSlopeLimiter(double steepness, double bias);
    static double computeMinBeeSlopeLimiter(double steepness, double bias);
    
    static double computeSlopeLimiter(double steepness, double bias, int slopeLimiter);
    
    static vector<double> computeSlopeVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                             vector<double> rightConservedVariableVector, double bias, int slopeLimiter);
    static vector<double> computeSlopeVector(EulerStateVector leftStateVector, EulerStateVector middleSttaeVector, EulerStateVector rightStateVector, double bias,
                                             int slopeLimiter, EulerMaterialParameters materialParameters);
    static vector<double> computeSlopeVector(RelativisticEulerStateVector leftStateVector, RelativisticEulerStateVector middleStateVector,
                                             RelativisticEulerStateVector rightStateVector, double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
    
    static vector<double> computeSlopeVectorX(BlackHoleEulerStateVector leftStateVector, BlackHoleEulerStateVector middleStateVector,
                                              BlackHoleEulerStateVector rightStateVector, double cellSpacing, double bias, int slopeLimiter, double xCoordinate,
                                              double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static vector<double> computeSlopeVectorY(BlackHoleEulerStateVector topStateVector, BlackHoleEulerStateVector middleStateVector,
                                              BlackHoleEulerStateVector bottomStateVector, double cellSpacing, double bias, int slopeLimiter, double xCoordinate,
                                              double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
};

#endif
