#ifndef BlackHoleSpacetime_hpp
#define BlackHoleSpacetime_hpp

#include "../Mathematics/VectorAlgebra.hpp"
#include "../Mathematics/MatrixAlgebra.hpp"

class BlackHoleSpacetime
{
public:
    BlackHoleSpacetime();
    BlackHoleSpacetime(double newBlackHoleMass);
    BlackHoleSpacetime(double newBlackHoleMass, double newBlackHoleXPosition, double newBlackHoleYPosition, double newBlackHoleZPosition);
    BlackHoleSpacetime(double newBlackHoleMass, double newBlackHoleSpin);
    BlackHoleSpacetime(double newBlackHoleMass, double newBlackHoleSpin, double newBlackHoleXPosition, double newBlackHoleYPosition, double newBlackHoleZPosition);
    
    double computeKerrSchildScalar(double xCoordinate, double yCoordinate, double zCoordinate);
    vector<double> computeKerrSchildVector(double xCoordinate, double yCoordinate, double zCoordinate);
    vector<double> computeKerrSchildSpacetimeVector(double xCoordinate, double yCoordinate, double zCoordinate);
    
    vector<double> computeKerrSchildScalarDerivative(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing);
    vector<vector<double> > computeKerrSchildVectorDerivative(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing);
    
    vector<vector<double> > computeSpatialMetricTensor(double xCoordinate, double yCoordinate, double zCoordinate);
    vector<vector<double> > computeSpacetimeMetricTensor(double xCoordinate, double yCoordinate, double zCoordinate);
    
    vector<vector<double> > computeInverseSpatialMetricTensor(double xCoordinate, double yCoordinate, double zCoordinate);
    vector<vector<double> > computeInverseSpacetimeMetricTensor(double xCoordiante, double yCoordinate, double zCoordinate);
    
    double computeSpatialMetricDeterminant(double xCoordinate, double yCoordinate, double zCoordinate);
    double computeSpacetimeMetricDeterminant(double xCoordinate, double yCoordinate, double zCoordinate);
    
    vector<vector<vector<double> > > computeSpatialMetricTensorDerivative(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing);
    
    double computeLapseFunction(double xCoordinate, double yCoordinate, double zCoordinate);
    vector<double> computeShiftVector(double xCoordinate, double yCoordinate, double zCoordinate);
    
    vector<double> computeLapseFunctionDerivative(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing);
    vector<vector<double> > computeShiftVectorDerivative(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing);
    
    vector<vector<double> > computeExtrinsicCurvatureTensor(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing);
    
    bool inExcisionRegion(double xCoordinate, double yCoordinate, double zCoordinate);
    
    void setBlackHoleMass(double newBlackHoleMass);
    void setBlackHoleSpin(double newBlackHoleSpin);
    void setBlackHoleXPosition(double newBlackHoleXPosition);
    void setBlackHoleYPosition(double newBlackHoleYPosition);
    void setBlackHoleZPosition(double newBlackHoleZPosition);
    
    double getBlackHoleMass();
    double getBlackHoleSpin();
    double getBlackHoleXPosition();
    double getBlackHoleYPosition();
    double getBlackHoleZPosition();
    
private:
    double blackHoleMass;
    double blackHoleSpin;
    double blackHoleXPosition;
    double blackHoleYPosition;
    double blackHoleZPosition;
};

#endif
