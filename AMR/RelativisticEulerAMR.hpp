#ifndef RelativisticEulerAMR_hpp
#define RelativisticEulerAMR_hpp

#include "../Solvers/RelativisticEulerFirstOrderSolver.hpp"
#include "../Solvers/RelativisticEulerSecondOrderSolver.hpp"
#include "AMRHelper.hpp"

class RelativisticEulerAMR
{
public:
    static void computeFORCETimeStepAMR(vector<RelativisticEulerStateVector> & currentCells, double cellSpacing, double timeStsep, vector<bool> AMRStructure,
                                        EulerMaterialParameters materialParameters);
    
    static void computeXFORCETimeStep2DAMR(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                           vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters);
    static void computeYFORCETimeStep2DAMR(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                           vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters);
    
    static void computeSLICTimeStepAMR(vector<RelativisticEulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                       vector<bool> AMRStructure, EulerMaterialParameters materialParameters);
    
    static void computeXSLICTimeStep2DAMR(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                          int slopeLimiter, vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters);
    static void computeYSLICTimeStep2DAMR(vector<vector<RelativisticEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                          int slopeLimiter, vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters);
    
    static tuple<vector<RelativisticEulerStateVector>, vector<RelativisticEulerStateVector>, vector<bool>, vector<bool>>
    solveLevel1AMR(vector<RelativisticEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance, int order,
                   EulerMaterialParameters materialParameters);
    static tuple<vector<vector<RelativisticEulerStateVector> >, vector<vector<RelativisticEulerStateVector> >, vector<vector<bool> >, vector<vector<bool> >>
    solve2DLevel1AMR(vector<vector<RelativisticEulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance,
                   int order, EulerMaterialParameters materialParameters);
    
    static tuple<vector<RelativisticEulerStateVector>, vector<RelativisticEulerStateVector>, vector<RelativisticEulerStateVector>, vector<bool>, vector<bool>,
    vector<bool>> solveLevel2AMR(vector<RelativisticEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance,
                                 int order, EulerMaterialParameters materialParameters);
    static tuple<vector<vector<RelativisticEulerStateVector> >, vector<vector<RelativisticEulerStateVector> >, vector<vector<RelativisticEulerStateVector> >,
    vector<vector<bool> >, vector<vector<bool> >, vector<vector<bool> >> solve2DLevel2AMR(vector<vector<RelativisticEulerStateVector> > initialCells,
                                                                                          double cellSpacing, double CFLCoefficient, double finalTime,
                                                                                          double AMRTolerance, int order, EulerMaterialParameters materialParameters);
};

#endif
