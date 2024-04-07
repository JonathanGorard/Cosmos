#ifndef RelativisticEulerFirstOrderSolver_hpp
#define RelativisticEulerFirstOrderSolver_hpp

#include "../Mathematics/VectorAlgebra.hpp"
#include "RelativisticEulerSolvers.hpp"
#include "EulerFirstOrderSolver.hpp"

class RelativisticEulerFirstOrderSolver
{
public:
    static vector<double> computeXLaxFriedrichsFlux(RelativisticEulerStateVector leftStateVector, RelativisticEulerStateVector rightStateVector, double cellSpacing,
                                                    double timeStep, EulerMaterialParameters materialParameters);
    static vector<double> computeXRichtmyerFlux(RelativisticEulerStateVector leftStateVector, RelativisticEulerStateVector rightStateVector, double cellSpacing,
                                                double timeStep, EulerMaterialParameters materialParameters);
    static vector<double> computeXFORCEFlux(RelativisticEulerStateVector leftStateVector, RelativisticEulerStateVector rightStateVector, double cellSpacing,
                                            double timeStep, EulerMaterialParameters materialParameters);
    
    static vector<double> computeYLaxFriedrichsFlux(RelativisticEulerStateVector topStateVector, RelativisticEulerStateVector bottomStateVector, double cellSpacing,
                                                    double timeStep, EulerMaterialParameters materialParameters);
    static vector<double> computeYRichtmyerFlux(RelativisticEulerStateVector topStateVector, RelativisticEulerStateVector bottomStateVector, double cellSpacing,
                                                double timeStep, EulerMaterialParameters materialParameters);
    static vector<double> computeYFORCEFlux(RelativisticEulerStateVector topStateVector, RelativisticEulerStateVector bottomStateVector, double cellSpacing,
                                            double timeStep, EulerMaterialParameters materialParameters);
    
    static void computeFORCETimeStep(vector<RelativisticEulerStateVector> & currentCells, double cellSpacing, double timeStep,
                                     EulerMaterialParameters materialParameters);
    
    static void computeXFORCETimeStep2D(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                        EulerMaterialParameters materialParameters);
    static void computeYFORCETimeStep2D(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                        EulerMaterialParameters materialParameters);
    
    static vector<RelativisticEulerStateVector> solve(vector<RelativisticEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                      EulerMaterialParameters materialParameters);
    
    static vector<vector<RelativisticEulerStateVector> > solve2D(vector<vector<RelativisticEulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient,
                                                                 double finalTime, EulerMaterialParameters materialParameters);
};

#endif
