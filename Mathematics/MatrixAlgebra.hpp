#ifndef MatrixAlgebra_hpp
#define MatrixAlgebra_hpp

#include <cmath>
#include <vector>
using namespace std;

class MatrixAlgebra
{
public:
    static vector<vector<double> > computeIdentityMatrix(int matrixDimension);
    
    static vector<vector<double> > addMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
    static vector<vector<double> > subtractMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
    static vector<vector<double> > multiplyMatrix(double scalar1, vector<vector<double> > matrix1);
    static vector<vector<double> > multiplyMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
    
    static double computeMatrixTrace(vector<vector<double> > matrix1);
    static double computeMatrixDeterminant(vector<vector<double> > matrix1);
    
    static vector<vector<double> > computeMatrixInverse(vector<vector<double> > matrix1);
};

#endif
