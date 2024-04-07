#include "MatrixAlgebra.hpp"

vector<vector<double> > MatrixAlgebra::computeIdentityMatrix(int matrixDimension)
{
    vector<vector<double> > identityMatrix(matrixDimension, vector<double>(matrixDimension));
    
    for (int i = 0; i < matrixDimension; i++)
    {
        for (int j = 0; j < matrixDimension; j++)
        {
            if (i == j)
            {
                identityMatrix[i][j] = 1.0;
            }
            else
            {
                identityMatrix[i][j] = 0.0;
            }
        }
    }
    
    return identityMatrix;
}

vector<vector<double> > MatrixAlgebra::addMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    long rowCount = matrix1.size();
    long columnCount = matrix1[0].size();
    vector<vector<double> > sumMatrix(rowCount, vector<double>(columnCount));
    
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            sumMatrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    
    return sumMatrix;
}

vector<vector<double> > MatrixAlgebra::subtractMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    return addMatrices(matrix1, multiplyMatrix(-1.0, matrix2));
}

vector<vector<double> > MatrixAlgebra::multiplyMatrix(double scalar1, vector<vector<double> > matrix1)
{
    long rowCount = matrix1.size();
    long columnCount = matrix1[0].size();
    vector<vector<double> > productMatrix(rowCount, vector<double>(columnCount));
    
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            productMatrix[i][j] = scalar1 * matrix1[i][j];
        }
    }
    
    return productMatrix;
}

vector<vector<double> > MatrixAlgebra::multiplyMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    long matrix1RowCount = matrix1.size();
    long matrix1ColumnCount = matrix1[0].size();
    long matrix2ColumnCount = matrix2[0].size();
    vector<vector<double> > productMatrix(matrix1RowCount, vector<double>(matrix2ColumnCount));
    
    for (int i = 0; i < matrix1RowCount; i++)
    {
        for (int j = 0; j < matrix2ColumnCount; j++)
        {
            productMatrix[i][j] = 0.0;
        }
    }
    
    for (int i = 0; i < matrix1RowCount; i++)
    {
        for (int j = 0; j < matrix2ColumnCount; j++)
        {
            for (int k = 0; k < matrix1ColumnCount; k++)
            {
                productMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    
    return productMatrix;
}

double MatrixAlgebra::computeMatrixTrace(vector<vector<double> > matrix1)
{
    long rowCount = matrix1.size();
    double trace = 0.0;
    
    for (int i = 0; i < rowCount; i++)
    {
        trace += matrix1[i][i];
    }
    
    return trace;
}

double MatrixAlgebra::computeMatrixDeterminant(vector<vector<double> > matrix1)
{
    long matrixDimension = matrix1.size();
    
    if (matrixDimension == 2)
    {
        return (matrix1[0][0] * matrix1[1][1]) - (matrix1[1][0] * matrix1[0][1]);
    }
    else
    {
        double matrixDeterminant = 0.0;
        
        for (int i = 0; i < matrixDimension; i++)
        {
            vector<vector<double> > cofactorMatrix(matrixDimension - 1, vector<double>(matrixDimension - 1));
            
            for (int j = 0; j < matrixDimension - 1; j++)
            {
                for (int k = 0; k < matrixDimension - 1; k++)
                {
                    if (k < i)
                    {
                        cofactorMatrix[j][k] = matrix1[j + 1][k];
                    }
                    else
                    {
                        cofactorMatrix[j][k] = matrix1[j + 1][k + 1];
                    }
                }
            }
            
            if (i % 2 == 0)
            {
                matrixDeterminant += matrix1[0][i] * computeMatrixDeterminant(cofactorMatrix);
            }
            else
            {
                matrixDeterminant -= matrix1[0][i] * computeMatrixDeterminant(cofactorMatrix);
            }
        }
        
        return matrixDeterminant;
    }
}

vector<vector<double> > MatrixAlgebra::computeMatrixInverse(vector<vector<double> > matrix1)
{
    long matrixDimension = matrix1.size();
    
    if (matrixDimension == 2)
    {
        double matrixDeterminant = computeMatrixDeterminant(matrix1);
        double matrixTrace = computeMatrixTrace(matrix1);
        vector<vector<double> > identityMatrix = computeIdentityMatrix(2);
        
        return multiplyMatrix((1.0 / matrixDeterminant), subtractMatrices(multiplyMatrix(matrixTrace, identityMatrix), matrix1));
    }
    else if (matrixDimension == 3)
    {
        double matrixDeterminant = computeMatrixDeterminant(matrix1);
        double matrixTrace = computeMatrixTrace(matrix1);
        double matrixSquaredTrace = computeMatrixTrace(multiplyMatrices(matrix1, matrix1));
        vector<vector<double> > identityMatrix = computeIdentityMatrix(3);
        
        return multiplyMatrix((1.0 / matrixDeterminant),
                              addMatrices(subtractMatrices(multiplyMatrix(0.5 * ((matrixTrace * matrixTrace) - matrixSquaredTrace), identityMatrix),
                                                           multiplyMatrix(matrixTrace, matrix1)), multiplyMatrices(matrix1, matrix1)));
    }
    else if (matrixDimension == 4)
    {
        double matrixDeterminant = computeMatrixDeterminant(matrix1);
        double matrixTrace = computeMatrixTrace(matrix1);
        double matrixSquaredTrace = computeMatrixTrace(multiplyMatrices(matrix1, matrix1));
        double matrixCubedTrace = computeMatrixTrace(multiplyMatrices(multiplyMatrices(matrix1, matrix1), matrix1));
        vector<vector<double> > identityMatrix = computeIdentityMatrix(4);
        
        return multiplyMatrix((1.0 / matrixDeterminant),
                              subtractMatrices(addMatrices(subtractMatrices(multiplyMatrix((1.0 / 6.0) * ((matrixTrace * matrixTrace * matrixTrace) -
                                                                                                          (3.0 * matrixTrace * matrixSquaredTrace) +
                                                                                                          (2.0 * matrixCubedTrace)), identityMatrix),
                                                                            multiplyMatrix(0.5 * ((matrixTrace * matrixTrace) - matrixSquaredTrace), matrix1)),
                                                           multiplyMatrix(matrixTrace, multiplyMatrices(matrix1, matrix1))),
                                               multiplyMatrices(multiplyMatrices(matrix1, matrix1), matrix1)));
    }
    else
    {
        return matrix1;
    }
}
