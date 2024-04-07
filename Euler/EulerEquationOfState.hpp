#ifndef EulerEquationOfState_hpp
#define EulerEquationOfState_hpp

#include "EulerMaterialParameters.hpp"
#include <cmath>
using namespace std;

class EulerEquationOfState
{
public:
    static double computeSpecificInternalEnergy(double density, double pressure, EulerMaterialParameters materialParameters);
    static double computePressure(double density, double specificInternalEnergy, EulerMaterialParameters materialParameters);
    static double computeSoundSpeed(double density, double pressure, EulerMaterialParameters materialParameters);
};

#endif
