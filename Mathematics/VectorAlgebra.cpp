#include "VectorAlgebra.hpp"

vector<double> VectorAlgebra::addVectors(vector<double> vector1, vector<double> vector2)
{
    long componentCount = vector1.size();
    vector<double> sumVector(componentCount);
    
#pragma omp parallel for
    for (int i = 0; i < componentCount; i++)
    {
        sumVector[i] = vector1[i] + vector2[i];
    }
    
    return sumVector;
}

vector<double> VectorAlgebra::subtractVectors(vector<double> vector1, vector<double> vector2)
{
    return addVectors(vector1, multiplyVector(-1.0, vector2));
}

vector<double> VectorAlgebra::multiplyVector(double scalar1, vector<double> vector1)
{
    long componentCount = vector1.size();
    vector<double> productVector(componentCount);
    
#pragma omp parallel for
    for (int i = 0; i < componentCount; i++)
    {
        productVector[i] = scalar1 * vector1[i];
    }
    
    return productVector;
}
