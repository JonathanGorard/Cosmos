#ifndef EulerAMR_hpp
#define EulerAMR_hpp

#include "../Solvers/EulerFirstOrderSolver.hpp"
#include "../Solvers/EulerSecondOrderSolver.hpp"
#include "AMRHelper.hpp"

class EulerAMR
{
public:
    static void computeFORCETimeStepAMR(vector<EulerStateVector> & currentCells, double cellSpacing, double timeStep, vector<bool> AMRStructure,
                                        EulerMaterialParameters materialParameters);
    
    static void computeXFORCETimeStep2DAMR(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, vector<vector<bool> > AMRStructure,
                                           EulerMaterialParameters materialParameters);
    static void computeYFORCETimeStep2DAMR(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, vector<vector<bool> > AMRStructure,
                                           EulerMaterialParameters materialParameters);
    
    static void computeSLICTimeStepAMR(vector<EulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                       vector<bool> AMRStructure, EulerMaterialParameters materialParameters);
    
    static void computeXSLICTimeStep2DAMR(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                          vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters);
    static void computeYSLICTimeStep2DAMR(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                          vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters);
    
    static tuple<vector<EulerStateVector>, vector<EulerStateVector>, vector<bool>, vector<bool>> solveLevel1AMR(vector<EulerStateVector> initialCells,
                                                                                                                double cellSpacing, double CFLCoefficient,
                                                                                                                double finalTime, double AMRTolerance, int order,
                                                                                                                EulerMaterialParameters materialParameters);
    static tuple<vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<bool> >, vector<vector<bool> >>
    solve2DLevel1AMR(vector<vector<EulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance, int order,
                     EulerMaterialParameters materialParameters);
    
    static tuple<vector<EulerStateVector>, vector<EulerStateVector>, vector<EulerStateVector>, vector<bool>, vector<bool>, vector<bool>>
    solveLevel2AMR(vector<EulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance, int order,
                   EulerMaterialParameters materialParameters);
    static tuple<vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<bool> >, vector<vector<bool> >,
    vector<vector<bool> >> solve2DLevel2AMR(vector<vector<EulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                            double AMRTolerance, int order, EulerMaterialParameters materialParameters);
    
};

#endif
