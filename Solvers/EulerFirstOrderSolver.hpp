#ifndef EulerFirstOrderSolver_hpp
#define EulerFirstOrderSolver_hpp

#include "../Euler/EulerStateVector.hpp"
#include "../Mathematics/VectorAlgebra.hpp"
#include "EulerSolvers.hpp"

class EulerFirstOrderSolver
{
public:
    static vector<double> computeLaxFriedrichsFlux(vector<double> leftConservedVariableVector, vector<double> rightConservedVariableVector, vector<double> leftFluxVector,
                                                   vector<double> rightFluxVector, double cellSpacing, double timeStep);
    static vector<double> computeRichtmyerFlux(vector<double> leftConservedVariableVector, vector<double> rightConservedVariableVector, vector<double> leftFluxVector,
                                               vector<double> rightFluxVector, double cellSpacing, double timeStep);
    static vector<double> computeFORCEFlux(vector<double> laxFriedrichsFlux, vector<double> richtmyerFlux);
    
    static vector<double> computeXLaxFriedrichsFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters materialParameters);
    static vector<double> computeXRichtmyerFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters materialParameters);
    static vector<double> computeXFORCEFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                            EulerMaterialParameters materialParameters);
    
    static vector<double> computeYLaxFriedrichsFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters materialParameters);
    static vector<double> computeYRichtmyerFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters materialParameters);
    static vector<double> computeYFORCEFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            EulerMaterialParameters materialParameters);
    
    static vector<double> computeFORCEUpdate(vector<double> conservedVariableVector, vector<double> leftFluxVector, vector<double> rightFluxVector,
                                             double cellSpacing, double timeStep);
    
    static void computeFORCETimeStep(vector<EulerStateVector> & currentCells, double cellSpacing, double timeStep,
                                     EulerMaterialParameters materialParameters);
    
    static void computeXFORCETimeStep2D(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                        EulerMaterialParameters materialParameters);
    static void computeYFORCETimeStep2D(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                        EulerMaterialParameters materialParameters);

    static vector<EulerStateVector> solve(vector<EulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                          EulerMaterialParameters materialParameters);
    
    static vector<vector<EulerStateVector> > solve2D(vector<vector<EulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                     EulerMaterialParameters materialParameters);
};

#endif
