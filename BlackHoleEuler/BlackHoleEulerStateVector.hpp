#ifndef BlackHoleEulerStateVector_hpp
#define BlackHoleEulerStateVector_hpp

#include "../RelativisticEuler/RelativisticEulerEquationOfState.hpp"
#include "BlackHoleSpacetime.hpp"
#include "BlackHoleEulerEigenvalues.hpp"

class BlackHoleEulerStateVector
{
public:
    BlackHoleEulerStateVector();
    BlackHoleEulerStateVector(double newDensity, double newXVelocity, double newYVelocity, double newZVelocity, double newPressure);
    
    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector);
    void setConservedVariableVector(vector<double> newConservedVariableVector, double xCoordinate, double yCoordinate, double zCoordinate,
                                    EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector(double xCoordinate, double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                                  BlackHoleSpacetime blackHole);
    
    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, double xCoordinate, double yCoordinate, double zCoordinate,
                                             EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    vector<double> computeXFluxVector(double xCoordinate, double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                      BlackHoleSpacetime blackHole);
    
    static vector<double> computeYFluxVector(vector<double> conservedVariableVector, double xCoordinate, double yCoordinate, double zCoordinate,
                                             EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    vector<double> computeYFluxVector(double xCoordinate, double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters,
                                      BlackHoleSpacetime blackHole);
    
    static vector<double> computeSourceTermVector(vector<double> conservedVariableVector, double xCoordinate, double yCoordinate, double zCoordinate,
                                                  double cellSpacing, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    vector<double> computeSourceTermVector(double xCoordinate, double yCoordinate, double zCoordinate, double cellSpacing,
                                           EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    double computeSoundSpeed(double xCoordinate, double yCoordinate, double zCoordinate, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    void setDensity(double newDensity);
    void setXVelocity(double newXVelocity);
    void setYVelocity(double newYVelocity);
    void setZVelocity(double newZVelocity);
    void setPressure(double newPressure);
    
    double getDensity();
    double getXVelocity();
    double getYVelocity();
    double getZVelocity();
    double getPressure();
    
private:
    double density;
    double xVelocity;
    double yVelocity;
    double zVelocity;
    double pressure;
};

#endif
