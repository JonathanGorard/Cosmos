#ifndef EulerStateVector_hpp
#define EulerStateVector_hpp

#include "EulerEquationOfState.hpp"
#include <vector>
using namespace std;

class EulerStateVector
{
public:
    EulerStateVector();
    EulerStateVector(double newDensity, double newXVelocity, double newYVelocity, double newZVelocity, double newPressure);
    
    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector);
    void setConservedVariableVector(vector<double> newConservedVariableVector, EulerMaterialParameters materialParameters);
    
    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector(EulerMaterialParameters materialParameters);
    
    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters materialParameters);
    vector<double> computeXFluxVector(EulerMaterialParameters materialParameters);
    
    static vector<double> computeYFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters materialParameters);
    vector<double> computeYFluxVector(EulerMaterialParameters materialParameters);
    
    double computeSpecificInternalEnergy(EulerMaterialParameters materialParameters);
    double computeTotalEnergy(EulerMaterialParameters materialParameters);
    double computeSoundSpeed(EulerMaterialParameters materialParameters);
    
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
