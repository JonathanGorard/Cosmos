#ifndef BlackHoleEulerEigenvalues_hpp
#define BlackHoleEulerEigenvalues_hpp

#include "BlackHoleSpacetime.hpp"
#include "../RelativisticEuler/RelativisticEulerEquationOfState.hpp"
#include "../Mathematics/VectorAlgebra.hpp"

class BlackHoleEulerEigenvalues
{
public:
    static vector<double> computeMaterialWaveEigenvalues(double xCoordinate, double yCoordinate, double zCoordinate, double xVelocity, double yVelocity,
                                                         double zVelocity, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static vector<double> computeFastAcousticWaveEigenvalues(double xCoordinate, double yCoordinate, double zCoordinate, double density, double pressure,
                                                             double xVelocity, double yVelocity, double zVelocity, EulerMaterialParameters materialParameters,
                                                             BlackHoleSpacetime blackHole);
    static vector<double> computeSlowAcousticWaveEigenvalues(double xCoordinate, double yCoordinate, double zCoordinate, double density, double pressure,
                                                             double xVelocity, double yVelocity, double zVelocity, EulerMaterialParameters materialParameters,
                                                             BlackHoleSpacetime blackHole);
};

#endif
