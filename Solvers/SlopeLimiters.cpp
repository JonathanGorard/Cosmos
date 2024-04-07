#include "SlopeLimiters.hpp"

double SlopeLimiters::computeGradientRatio(double steepness, double bias)
{
    return 2.0 / ((1.0 - bias) + ((1.0 + bias) * steepness));
}

double SlopeLimiters::computeSuperBeeSlopeLimiter(double steepness, double bias)
{
    if (steepness < 0.0)
    {
        return 0.0;
    }
    else if (steepness>= 0.0 && steepness < 0.5)
    {
        return 2.0 * steepness;
    }
    else if (steepness >= 0.5 && steepness < 1.0)
    {
        return 1.0;
    }
    else{
        return min(min(steepness, computeGradientRatio(steepness, bias)), 2.0);
    }
}

double SlopeLimiters::computeVanLeerSlopeLimiter(double steepness, double bias)
{
    if (steepness < 0.0)
    {
        return 0.0;
    }
    else
    {
        return min((2.0 * steepness) / (1.0 + steepness), computeGradientRatio(steepness, bias));
    }
}

double SlopeLimiters::computeMinBeeSlopeLimiter(double steepness, double bias)
{
    if (steepness < 0.0)
    {
        return 0.0;
    }
    else if (steepness >= 0.0 && steepness < 1.0)
    {
        return steepness;
    }
    else
    {
        return min(1.0, computeGradientRatio(steepness, bias));
    }
}

double SlopeLimiters::computeSlopeLimiter(double steepness, double bias, int slopeLimiter)
{
    if (slopeLimiter == 0)
    {
        return computeSuperBeeSlopeLimiter(steepness, bias);
    }
    else if (slopeLimiter == 1)
    {
        return computeVanLeerSlopeLimiter(steepness, bias);
    }
    else
    {
        return computeMinBeeSlopeLimiter(steepness, bias);
    }
}

vector<double> SlopeLimiters::computeSlopeVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                 vector<double> rightConservedVariableVector, double bias, int slopeLimiter)
{
    vector<double> leftConservedVariableVectorDifference = VectorAlgebra::subtractVectors(middleConservedVariableVector, leftConservedVariableVector);
    vector<double> rightConservedVariableVectorDifference = VectorAlgebra::subtractVectors(rightConservedVariableVector, middleConservedVariableVector);
    
    vector<double> slopeVector = VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(0.5 * (1.0 + bias), leftConservedVariableVectorDifference),
                                                           VectorAlgebra::multiplyVector(0.5 * (1.0 - bias), rightConservedVariableVectorDifference));
    long conservedVariableCount = slopeVector.size();
    
#pragma omp parallel for
    for (int i = 0; i < conservedVariableCount; i++)
    {
        double numerator = leftConservedVariableVectorDifference[i];
        double denominator = rightConservedVariableVectorDifference[i];
        
        if (abs(numerator) < pow(10.0, -5.0))
        {
            numerator = pow(10.0, -5.0);
        }
        if (abs(denominator) < pow(10.0, -5.0))
        {
            denominator = pow(10.0, -5.0);
        }
        
        double steepness = numerator / denominator;
        slopeVector[i] *= computeSlopeLimiter(steepness, bias, slopeLimiter);
    }
    
    return slopeVector;
}

vector<double> SlopeLimiters::computeSlopeVector(EulerStateVector leftStateVector, EulerStateVector middleStateVector, EulerStateVector rightStateVector,
                                                 double bias, int slopeLimiter, EulerMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);
    
    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(RelativisticEulerStateVector leftStateVector, RelativisticEulerStateVector middleStateVector,
                                                 RelativisticEulerStateVector rightStateVector, double bias, int slopeLimiter, EulerMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);
    
    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVectorX(BlackHoleEulerStateVector leftStateVector, BlackHoleEulerStateVector middleStateVector,
                                                  BlackHoleEulerStateVector rightStateVector, double cellSpacing, double bias, int slopeLimiter, double xCoordinate,
                                                  double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(xCoordinate - cellSpacing, yCoordinate, zCoordinate,
                                                                                                materialParameters, blackHole);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(xCoordinate + cellSpacing, yCoordinate, zCoordinate,
                                                                                                  materialParameters, blackHole);
    
    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVectorY(BlackHoleEulerStateVector topStateVector, BlackHoleEulerStateVector middleStateVector,
                                                  BlackHoleEulerStateVector bottomStateVector, double cellSpacing, double bias, int slopeLimiter, double xCoordinate,
                                                  double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(xCoordinate, yCoordinate - cellSpacing, zCoordinate,
                                                                                              materialParameters, blackHole);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(xCoordinate, yCoordinate, zCoordinate, materialParameters, blackHole);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(xCoordinate, yCoordinate + cellSpacing, zCoordinate,
                                                                                                    materialParameters, blackHole);
    
    return computeSlopeVector(topConservedVariableVector, middleConservedVariableVector, bottomConservedVariableVector, bias, slopeLimiter);
}
