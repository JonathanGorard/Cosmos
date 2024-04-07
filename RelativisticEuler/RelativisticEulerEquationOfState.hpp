#ifndef RelativisticEulerEquationOfState_hpp
#define RelativisticEulerEquationOfState_hpp

#include "../Euler/EulerMaterialParameters.hpp"
#include <cmath>
using namespace std;

class RelativisticEulerEquationOfState
{
public:
    static double computeSpecificEnthalpy(double density, double pressure, EulerMaterialParameters materialParameters);
    static double computeSoundSpeed(double density, double pressure, EulerMaterialParameters materialParameters);
};

#endif
