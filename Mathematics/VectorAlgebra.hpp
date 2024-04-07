#ifndef VectorAlgebra_hpp
#define VectorAlgebra_hpp

#include <cmath>
#include <vector>
using namespace std;

class VectorAlgebra
{
public:
    static vector<double> addVectors(vector<double> vector1, vector<double> vector2);
    static vector<double> subtractVectors(vector<double> vector1, vector<double> vector2);
    
    static vector<double> multiplyVector(double scalar1, vector<double> vector1);
};

#endif
