#ifndef BlackHoleEulerFirstOrderSolver_hpp
#define BlackHoleEulerFirstOrderSolver_hpp

#include "../Mathematics/VectorAlgebra.hpp"
#include "BlackHoleEulerForcingSolver.hpp"
#include "EulerFirstOrderSolver.hpp"

class BlackHoleEulerFirstOrderSolver
{
public:
    static vector<double> computeXLaxFriedrichsFlux(BlackHoleEulerStateVector leftStateVector, BlackHoleEulerStateVector rightStateVector, double xCoordinate,
                                                    double yCoordinate, double zCoordinate, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static vector<double> computeXRichtmyerFlux(BlackHoleEulerStateVector leftStateVector, BlackHoleEulerStateVector rightStateVector, double xCoordinate,
                                                double yCoordinate, double zCoordinate, double cellSpacing, double timeSep,
                                                EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static vector<double> computeXFORCEFlux(BlackHoleEulerStateVector leftStateVector, BlackHoleEulerStateVector rightStateVector, double xCoordinate,
                                            double yCoordinate, double zCoordinate, double cellSpacing, double timeStep,
                                            EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static vector<double> computeYLaxFriedrichsFlux(BlackHoleEulerStateVector topStateVector, BlackHoleEulerStateVector bottomStateVector, double xCoordinate,
                                                    double yCoordinate, double zCoordinate, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static vector<double> computeYRichtmyerFlux(BlackHoleEulerStateVector topStateVector, BlackHoleEulerStateVector bottomStateVector, double xCoordinate,
                                                double yCoordinate, double zCoordinate, double cellSpacing, double timeStep,
                                                EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static vector<double> computeYFORCEFlux(BlackHoleEulerStateVector topStateVector, BlackHoleEulerStateVector bottomStateVector, double xCoordinate,
                                            double yCoordinate, double zCoordinate, double cellSpacing, double timeStep,
                                            EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static void computeFORCETimeStep(vector<BlackHoleEulerStateVector> & currentCells, double cellSpacing, double timeStep, EulerMaterialParameters materialParameters,
                                     BlackHoleSpacetime blackHole);
    
    static void computeXFORCETimeStep2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                        EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static void computeYFORCETimeStep2D(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                        EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static vector<BlackHoleEulerStateVector> solve(vector<BlackHoleEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                   int subcyclingIterations, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static vector<vector<BlackHoleEulerStateVector> > solve2D(vector<vector<BlackHoleEulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient,
                                                              double finalTime, int subcyclingIterations, EulerMaterialParameters materialParameters,
                                                              BlackHoleSpacetime blackHole);
};

#endif
