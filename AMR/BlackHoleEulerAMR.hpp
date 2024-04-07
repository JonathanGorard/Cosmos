#ifndef BlackHoleEulerAMR_hpp
#define BlackHoleEulerAMR_hpp

#include "../Solvers/BlackHoleEulerFirstOrderSolver.hpp"
#include "../Solvers/BlackHoleEulerSecondOrderSolver.hpp"
#include "AMRHelper.hpp"

class BlackHoleEulerAMR
{
public:
    static void computeFORCETimeStepAMR(vector<BlackHoleEulerStateVector> & currentCells, double cellSpacing, double timeStep, vector<bool> AMRStructure,
                                        EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static void computeXFORCETimeStep2DAMR(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                           vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static void computeYFORCETimeStep2DAMR(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep,
                                           vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static void computeSLICTimeStepAMR(vector<BlackHoleEulerStateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                       vector<bool> AMRStructure, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static void computeXSLICTimeStep2DAMR(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                          int slopeLimiter, vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static void computeYSLICTimeStep2DAMR(vector<vector<BlackHoleEulerStateVector> > & currentCells, double cellSpacing, double timeStep, double bias,
                                          int slopeLimiter, vector<vector<bool> > AMRStructure, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static tuple<vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<bool>, vector<bool>>
    solveLevel1AMR(vector<BlackHoleEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance, int order,
                   EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static tuple<vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >, vector<vector<bool> >, vector<vector<bool> >>
    solve2DLevel1AMR(vector<vector<BlackHoleEulerStateVector> > initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance,
                     int order, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    
    static tuple<vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<bool>, vector<bool>,
    vector<bool>> solveLevel2AMR(vector<BlackHoleEulerStateVector> initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double AMRTolerance,
                                 int order, EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
    static tuple<vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >,
    vector<vector<bool> >, vector<vector<bool> >, vector<vector<bool> >> solve2DLevel2AMR(vector<vector<BlackHoleEulerStateVector> > initialCells, double cellSpacing,
                                                                                          double CFLCoefficient, double finalTime, double AMRTolerance, int order,
                                                                                          EulerMaterialParameters materialParameters, BlackHoleSpacetime blackHole);
};

#endif
