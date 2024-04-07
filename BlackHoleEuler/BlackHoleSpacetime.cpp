#include "BlackHoleSpacetime.hpp"

BlackHoleSpacetime::BlackHoleSpacetime()
{
    blackHoleMass = 1.0;
    blackHoleSpin = 0.0;
    
    blackHoleXPosition = 0.5;
    blackHoleYPosition = 0.5;
    blackHoleZPosition = 0.0;
}

BlackHoleSpacetime::BlackHoleSpacetime(double newBlackHoleMass)
{
    blackHoleMass = newBlackHoleMass;
    blackHoleSpin = 0.0;
    
    blackHoleXPosition = 0.5;
    blackHoleYPosition = 0.5;
    blackHoleZPosition = 0.0;
}

BlackHoleSpacetime::BlackHoleSpacetime(double newBlackHoleMass, double newBlackHoleXPosition, double newBlackHoleYPosition, double newBlackHoleZPosition)
{
    blackHoleMass = newBlackHoleMass;
    blackHoleSpin = 0.0;
    
    blackHoleXPosition = newBlackHoleXPosition;
    blackHoleYPosition = newBlackHoleYPosition;
    blackHoleZPosition = newBlackHoleZPosition;
}

BlackHoleSpacetime::BlackHoleSpacetime(double newBlackHoleMass, double newBlackHoleSpin)
{
    blackHoleMass = newBlackHoleMass;
    blackHoleSpin = newBlackHoleSpin;
    
    blackHoleXPosition = 0.5;
    blackHoleYPosition = 0.5;
    blackHoleZPosition = 0.0;
}

BlackHoleSpacetime::BlackHoleSpacetime(double newBlackHoleMass, double newBlackHoleSpin, double newBlackHoleXPosition, double newBlackHoleYPosition,
                                       double newBlackHoleZPosition)
{
    blackHoleMass = newBlackHoleMass;
    blackHoleSpin = newBlackHoleSpin;
    
    blackHoleXPosition = newBlackHoleXPosition;
    blackHoleYPosition = newBlackHoleYPosition;
    blackHoleZPosition = newBlackHoleZPosition;
}

double BlackHoleSpacetime::computeKerrSchildScalar(double xCoordinate, double yCoordinate, double zCoordinate)
{
    double normSquared = ((xCoordinate - blackHoleXPosition) * (xCoordinate - blackHoleXPosition)) +
    ((yCoordinate - blackHoleYPosition) * (yCoordinate - blackHoleYPosition)) + ((zCoordinate - blackHoleZPosition) * (zCoordinate - blackHoleZPosition));
    
    double radialDistance = sqrt(0.5 * (normSquared - ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin)) +
                                        sqrt(((normSquared - ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin))) *
                                              (normSquared - ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin)))) +
                                             (4.0 * ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin) *
                                                     ((zCoordinate - blackHoleZPosition) * (zCoordinate - blackHoleZPosition)))))));
    
    return - (blackHoleMass * (radialDistance * radialDistance * radialDistance)) /
    ((radialDistance * radialDistance * radialDistance * radialDistance) + ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin) *
                                                                            ((zCoordinate - blackHoleZPosition) * (zCoordinate - blackHoleZPosition))));
}

vector<double> BlackHoleSpacetime::computeKerrSchildVector(double xCoordinate, double yCoordinate, double zCoordinate)
{
    double normSquared = ((xCoordinate - blackHoleXPosition) * (xCoordinate - blackHoleXPosition)) +
    ((yCoordinate - blackHoleYPosition) * (yCoordinate - blackHoleYPosition)) + ((zCoordinate - blackHoleZPosition) * (zCoordinate - blackHoleZPosition));
    
    double radialDistance = sqrt(0.5 * (normSquared - ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin)) +
                                        sqrt(((normSquared - ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin))) *
                                              (normSquared - ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin)))) +
                                             (4.0 * ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin) *
                                                     ((zCoordinate - blackHoleZPosition) * (zCoordinate - blackHoleZPosition)))))));
    
    vector<double> kerrSchildVector(3);
    
    kerrSchildVector[0] = - ((radialDistance * (xCoordinate - blackHoleXPosition)) + (blackHoleSpin * blackHoleMass * (yCoordinate - blackHoleYPosition))) /
    ((radialDistance * radialDistance) + ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin)));
    kerrSchildVector[1] = - ((radialDistance * (yCoordinate - blackHoleYPosition)) - (blackHoleSpin * blackHoleMass * (xCoordinate - blackHoleXPosition))) /
    ((radialDistance * radialDistance) + ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin)));
    kerrSchildVector[2] = - (zCoordinate - blackHoleZPosition) / radialDistance;
    
    return kerrSchildVector;
}

vector<double> BlackHoleSpacetime::computeKerrSchildSpacetimeVector(double xCoordinate, double yCoordinate, double zCoordinate)
{
    double normSquared = ((xCoordinate - blackHoleXPosition) * (xCoordinate - blackHoleXPosition)) +
    ((yCoordinate - blackHoleYPosition) * (yCoordinate - blackHoleYPosition)) + ((zCoordinate - blackHoleZPosition) * (zCoordinate - blackHoleZPosition));
    
    double radialDistance = sqrt(0.5 * (normSquared - ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin)) +
                                        sqrt(((normSquared - ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin))) *
                                              (normSquared - ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin)))) +
                                             (4.0 * ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin) *
                                                     ((zCoordinate - blackHoleZPosition) * (zCoordinate - blackHoleZPosition)))))));
    
    vector<double> kerrSchildSpacetimeVector(4);
    
    kerrSchildSpacetimeVector[0] = -1.0;
    kerrSchildSpacetimeVector[1] = - ((radialDistance * (xCoordinate - blackHoleXPosition)) + (blackHoleSpin * blackHoleMass * (yCoordinate - blackHoleYPosition))) /
    ((radialDistance * radialDistance) + ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin)));
    kerrSchildSpacetimeVector[2] = - ((radialDistance + (yCoordinate - blackHoleYPosition)) - (blackHoleSpin * blackHoleMass * (xCoordinate - blackHoleXPosition))) /
    ((radialDistance * radialDistance) + ((blackHoleMass * blackHoleMass) * (blackHoleSpin * blackHoleSpin)));
    kerrSchildSpacetimeVector[3] = - (zCoordinate - blackHoleZPosition) / radialDistance;
    
    return kerrSchildSpacetimeVector;
}

vector<double> BlackHoleSpacetime::computeKerrSchildScalarDerivative(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing)
{
    vector<double> kerrSchildScalarDerivative(3);
    
    kerrSchildScalarDerivative[0] = (computeKerrSchildScalar(xCoordinate + (0.5 * cellSpacing), yCoordinate, zCoordinate) -
                                     computeKerrSchildScalar(xCoordinate - (0.5 * cellSpacing), yCoordinate, zCoordinate)) / cellSpacing;
    kerrSchildScalarDerivative[1] = (computeKerrSchildScalar(xCoordinate, yCoordinate + (0.5 * cellSpacing), zCoordinate) -
                                     computeKerrSchildScalar(xCoordinate, yCoordinate - (0.5 * cellSpacing), zCoordinate)) / cellSpacing;
    kerrSchildScalarDerivative[2] = (computeKerrSchildScalar(xCoordinate, yCoordinate, zCoordinate + (0.5 * cellSpacing)) -
                                     computeKerrSchildScalar(xCoordinate, yCoordinate, zCoordinate - (0.5 * cellSpacing))) / cellSpacing;
    
    return kerrSchildScalarDerivative;
}

vector<vector<double> > BlackHoleSpacetime::computeKerrSchildVectorDerivative(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing)
{
    vector<vector<double> > kerrSchildVectorDerivative(3, vector<double>(3));
    
    vector<double> kerrSchildVectorXDerivative = VectorAlgebra::multiplyVector((1.0 / cellSpacing),
        VectorAlgebra::subtractVectors(computeKerrSchildVector(xCoordinate + (0.5 * cellSpacing), yCoordinate, zCoordinate),
                                       computeKerrSchildVector(xCoordinate - (0.5 * cellSpacing), yCoordinate, zCoordinate)));
    vector<double> kerrSchildVectorYDerivative = VectorAlgebra::multiplyVector((1.0 / cellSpacing),
        VectorAlgebra::subtractVectors(computeKerrSchildVector(xCoordinate, yCoordinate + (0.5 * cellSpacing), zCoordinate),
                                       computeKerrSchildVector(xCoordinate, yCoordinate - (0.5 * cellSpacing), zCoordinate)));
    vector<double> kerrSchildVectorZDerivative = VectorAlgebra::multiplyVector((1.0 / cellSpacing),
        VectorAlgebra::subtractVectors(computeKerrSchildVector(xCoordinate, yCoordinate, zCoordinate + (0.5 * cellSpacing)),
                                       computeKerrSchildVector(xCoordinate, yCoordinate, zCoordinate - (0.5 * cellSpacing))));
    
    for (int i = 0; i < 3; i++)
    {
        kerrSchildVectorDerivative[0][i] = kerrSchildVectorXDerivative[i];
        kerrSchildVectorDerivative[1][i] = kerrSchildVectorYDerivative[i];
        kerrSchildVectorDerivative[2][i] = kerrSchildVectorZDerivative[i];
    }
    
    return kerrSchildVectorDerivative;
}

vector<vector<double> > BlackHoleSpacetime::computeSpatialMetricTensor(double xCoordinate, double yCoordinate, double zCoordinate)
{
    vector<vector<double> > euclideanMetricTensor(3, vector<double>(3));
    
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i == j)
            {
                euclideanMetricTensor[i][j] = 1.0;
            }
            else
            {
                euclideanMetricTensor[i][j] = 0.0;
            }
        }
    }
    
    double kerrSchildScalar = computeKerrSchildScalar(xCoordinate, yCoordinate, zCoordinate);
    vector<double> kerrSchildVector = computeKerrSchildVector(xCoordinate, yCoordinate, zCoordinate);
    
    vector<vector<double> > spatialMetricTensor(3, vector<double>(3));
    
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            spatialMetricTensor[i][j] = euclideanMetricTensor[i][j] - (2.0 * kerrSchildScalar * kerrSchildVector[i] * kerrSchildVector[j]);
        }
    }
    
    return spatialMetricTensor;
}

vector<vector<double> > BlackHoleSpacetime::computeSpacetimeMetricTensor(double xCoordinate, double yCoordinate, double zCoordinate)
{
    vector<vector<double> > minkowskiMetricTensor(4, vector<double>(4));
    
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                if (i == 0 && j == 0)
                {
                    minkowskiMetricTensor[i][j] = -1.0;
                }
                else
                {
                    minkowskiMetricTensor[i][j] = 1.0;
                }
            }
            else
            {
                minkowskiMetricTensor[i][j] = 0.0;
            }
        }
    }
    
    double kerrSchildScalar = computeKerrSchildScalar(xCoordinate, yCoordinate, zCoordinate);
    vector<double> kerrSchildSpacetimeVector = computeKerrSchildSpacetimeVector(xCoordinate, yCoordinate, zCoordinate);
    
    vector<vector<double> > spacetimeMetricTensor(4, vector<double>(4));
    
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            spacetimeMetricTensor[i][j] = minkowskiMetricTensor[i][j] - (2.0 * kerrSchildScalar * kerrSchildSpacetimeVector[i] * kerrSchildSpacetimeVector[j]);
        }
    }
    
    return spacetimeMetricTensor;
}

vector<vector<double> > BlackHoleSpacetime::computeInverseSpatialMetricTensor(double xCoordinate, double yCoordinate, double zCoordinate)
{
    return MatrixAlgebra::computeMatrixInverse(computeSpatialMetricTensor(xCoordinate, yCoordinate, zCoordinate));
}

vector<vector<double> > BlackHoleSpacetime::computeInverseSpacetimeMetricTensor(double xCoordinate, double yCoordinate, double zCoordinate)
{
    return MatrixAlgebra::computeMatrixInverse(computeSpacetimeMetricTensor(xCoordinate, yCoordinate, zCoordinate));
}

double BlackHoleSpacetime::computeSpatialMetricDeterminant(double xCoordinate, double yCoordinate, double zCoordinate)
{
    return MatrixAlgebra::computeMatrixDeterminant(computeSpatialMetricTensor(xCoordinate, yCoordinate, zCoordinate));
}

double BlackHoleSpacetime::computeSpacetimeMetricDeterminant(double xCoordinate, double yCoordinate, double zCoordinate)
{
    return MatrixAlgebra::computeMatrixDeterminant(computeSpacetimeMetricTensor(xCoordinate, yCoordinate, zCoordinate));
}

vector<vector<vector<double> > > BlackHoleSpacetime::computeSpatialMetricTensorDerivative(double xCoordinate, double yCoordinate, double zCoordinate,
                                                                                          double cellSpacing)
{
    vector<vector<vector<double> > > spatialMetricTensorDerivative(3, vector<vector<double> >(3, vector<double>(3)));
    
    vector<vector<double> > spatialMetricTensorXDerivative = MatrixAlgebra::multiplyMatrix((1.0 / cellSpacing),
        MatrixAlgebra::subtractMatrices(computeSpatialMetricTensor(xCoordinate + (0.5 * cellSpacing), yCoordinate, zCoordinate),
                                        computeSpatialMetricTensor(xCoordinate - (0.5 * cellSpacing), yCoordinate, zCoordinate)));
    vector<vector<double> > spatialMetricTensorYDerivative = MatrixAlgebra::multiplyMatrix((1.0 / cellSpacing),
        MatrixAlgebra::subtractMatrices(computeSpatialMetricTensor(xCoordinate, yCoordinate + (0.5 * cellSpacing), zCoordinate),
                                        computeSpatialMetricTensor(xCoordinate, yCoordinate - (0.5 * cellSpacing), zCoordinate)));
    vector<vector<double> > spatialMetricTensorZDerivative = MatrixAlgebra::multiplyMatrix((1.0 / cellSpacing),
        MatrixAlgebra::subtractMatrices(computeSpatialMetricTensor(xCoordinate, yCoordinate, zCoordinate + (0.5 * cellSpacing)),
                                        computeSpatialMetricTensor(xCoordinate, yCoordinate, zCoordinate - (0.5 * cellSpacing))));
    
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            spatialMetricTensorDerivative[0][i][j] = spatialMetricTensorXDerivative[i][j];
            spatialMetricTensorDerivative[1][i][j] = spatialMetricTensorYDerivative[i][j];
            spatialMetricTensorDerivative[2][i][j] = spatialMetricTensorZDerivative[i][j];
        }
    }
    
    return spatialMetricTensorDerivative;
}

double BlackHoleSpacetime::computeLapseFunction(double xCoordinate, double yCoordinate, double zCoordinate)
{
    return 1.0 / sqrt(1.0 - (2.0 * computeKerrSchildScalar(xCoordinate, yCoordinate, zCoordinate)));
}

vector<double> BlackHoleSpacetime::computeShiftVector(double xCoordinate, double yCoordinate, double zCoordinate)
{
    double kerrSchildScalar = computeKerrSchildScalar(xCoordinate, yCoordinate, zCoordinate);
    vector<double> kerrSchildVector = computeKerrSchildVector(xCoordinate, yCoordinate, zCoordinate);
    
    vector<double> shiftVector(3);
    
    for (int i = 0; i < 3; i++)
    {
        shiftVector[i] = ((2.0 * kerrSchildScalar) / (1.0 - (2.0 * kerrSchildScalar))) * kerrSchildVector[i];
    }
    
    return shiftVector;
}

vector<double> BlackHoleSpacetime::computeLapseFunctionDerivative(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing)
{
    vector<double> lapseFunctionDerivative(3);
    
    lapseFunctionDerivative[0] = (computeLapseFunction(xCoordinate + (0.5 * cellSpacing), yCoordinate, zCoordinate) -
                                  computeLapseFunction(xCoordinate - (0.5 * cellSpacing), yCoordinate, zCoordinate)) / cellSpacing;
    lapseFunctionDerivative[1] = (computeLapseFunction(xCoordinate, yCoordinate + (0.5 * cellSpacing), zCoordinate) -
                                  computeLapseFunction(xCoordinate, yCoordinate - (0.5 * cellSpacing), zCoordinate)) / cellSpacing;
    lapseFunctionDerivative[2] = (computeLapseFunction(xCoordinate, yCoordinate, zCoordinate + (0.5 * cellSpacing)) -
                                  computeLapseFunction(xCoordinate, yCoordinate, zCoordinate - (0.5 * cellSpacing))) / cellSpacing;
    
    return lapseFunctionDerivative;
}

vector<vector<double> > BlackHoleSpacetime::computeShiftVectorDerivative(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing)
{
    vector<vector<double> > shiftVectorDerivative(3, vector<double>(3));
    
    vector<double> shiftVectorXDerivative = VectorAlgebra::multiplyVector((1.0 / cellSpacing),
        VectorAlgebra::subtractVectors(computeShiftVector(xCoordinate + (0.5 * cellSpacing), yCoordinate, zCoordinate),
                                       computeShiftVector(xCoordinate - (0.5 * cellSpacing), yCoordinate, zCoordinate)));
    vector<double> shiftVectorYDerivative = VectorAlgebra::multiplyVector((1.0 / cellSpacing),
        VectorAlgebra::subtractVectors(computeShiftVector(xCoordinate, yCoordinate + (0.5 * cellSpacing), zCoordinate),
                                       computeShiftVector(xCoordinate, yCoordinate - (0.5 * cellSpacing), zCoordinate)));
    vector<double> shiftVectorZDerivative = VectorAlgebra::multiplyVector((1.0 / cellSpacing),
        VectorAlgebra::subtractVectors(computeShiftVector(xCoordinate, yCoordinate, zCoordinate + (0.5 * cellSpacing)),
                                       computeShiftVector(xCoordinate, yCoordinate, zCoordinate - (0.5 * cellSpacing))));
    
    for (int i = 0; i < 3; i++)
    {
        shiftVectorDerivative[0][i] = shiftVectorXDerivative[i];
        shiftVectorDerivative[1][i] = shiftVectorYDerivative[i];
        shiftVectorDerivative[2][i] = shiftVectorZDerivative[i];
    }
    
    return shiftVectorDerivative;
}

vector<vector<double> > BlackHoleSpacetime::computeExtrinsicCurvatureTensor(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing)
{
    vector<vector<double> > extrinsicCurvatureTensor(3, vector<double>(3));
    
    double lapseFunction = computeLapseFunction(xCoordinate, yCoordinate, zCoordinate);
    
    double kerrSchildScalar = computeKerrSchildScalar(xCoordinate, yCoordinate, zCoordinate);
    vector<double> kerrSchildVector = computeKerrSchildVector(xCoordinate, yCoordinate, zCoordinate);
    
    vector<double> kerrSchildScalarDerivative = computeKerrSchildScalarDerivative(xCoordinate, yCoordinate, zCoordinate, cellSpacing);
    vector<vector<double> > kerrSchildVectorDerivative = computeKerrSchildVectorDerivative(xCoordinate, yCoordinate, zCoordinate, cellSpacing);
    
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            extrinsicCurvatureTensor[i][j] = lapseFunction * (- (kerrSchildVector[i] * kerrSchildScalarDerivative[j]) -
                                                              (kerrSchildVector[j] * kerrSchildScalarDerivative[i]) -
                                                              (kerrSchildScalar * kerrSchildVectorDerivative[j][i]) -
                                                              (kerrSchildScalar * kerrSchildVectorDerivative[i][j]));
            
            for (int k = 0; k < 3; k++)
            {
                extrinsicCurvatureTensor[i][j] += lapseFunction * 2.0 * (kerrSchildScalar * kerrSchildScalar) *
                ((kerrSchildVector[i] * kerrSchildVector[k] * kerrSchildVectorDerivative[k][j]) +
                 (kerrSchildVector[j] * kerrSchildVector[k] * kerrSchildVectorDerivative[k][i]));
                
                extrinsicCurvatureTensor[i][j] += lapseFunction * 2.0 * kerrSchildScalar * kerrSchildVector[i] * kerrSchildVector[j] * kerrSchildVector[k] *
                kerrSchildScalarDerivative[k];
            }
        }
    }
    
    return extrinsicCurvatureTensor;
}

bool BlackHoleSpacetime::inExcisionRegion(double xCoordinate, double yCoordinate, double zCoordinate)
{
    double radialCoordinate = sqrt(((xCoordinate - blackHoleXPosition) * (xCoordinate - blackHoleXPosition)) +
                                   ((yCoordinate - blackHoleYPosition) * (yCoordinate - blackHoleYPosition)) +
                                   ((zCoordinate - blackHoleZPosition) * (zCoordinate - blackHoleZPosition)));
    
    if (radialCoordinate <= (blackHoleMass * (1.0 + sqrt(1.0 - (blackHoleSpin * blackHoleSpin)))))
    {
        return true;
    }
    else
    {
        return false;
    }
}

void BlackHoleSpacetime::setBlackHoleMass(double newBlackHoleMass)
{
    blackHoleMass = newBlackHoleMass;
}

void BlackHoleSpacetime::setBlackHoleSpin(double newBlackHoleSpin)
{
    blackHoleSpin = newBlackHoleSpin;
}

void BlackHoleSpacetime::setBlackHoleXPosition(double newBlackHoleXPosition)
{
    blackHoleXPosition = newBlackHoleXPosition;
}

void BlackHoleSpacetime::setBlackHoleYPosition(double newBlackHoleYPosition)
{
    blackHoleYPosition = newBlackHoleYPosition;
}

void BlackHoleSpacetime::setBlackHoleZPosition(double newBlackHoleZPosition)
{
    blackHoleZPosition = newBlackHoleZPosition;
}

double BlackHoleSpacetime::getBlackHoleMass()
{
    return blackHoleMass;
}

double BlackHoleSpacetime::getBlackHoleSpin()
{
    return blackHoleSpin;
}

double BlackHoleSpacetime::getBlackHoleXPosition()
{
    return blackHoleXPosition;
}

double BlackHoleSpacetime::getBlackHoleYPosition()
{
    return blackHoleYPosition;
}

double BlackHoleSpacetime::getBlackHoleZPosition()
{
    return blackHoleZPosition;
}
