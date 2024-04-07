#ifndef RelativisticEulerAMRTests_hpp
#define RelativisticEulerAMRTests_hpp

#include "../AMR/RelativisticEulerAMR.hpp"
#include <fstream>
using namespace std;

class RelativisticEulerAMRTests
{
public:
    static void solveMildlyRelativisticShock(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solveStronglyRelativisticBlast(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solvePerturbedDensityTest(int cellCount, int order, double AMRTolerance, int AMRLevel);
    
    static void solve2DMildlyRelativisticShock(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solve2DStronglyRelativisticBlast(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solve2DPerturbedDensityTest(int cellCount, int order, double AMRTolerance, int AMRLevel);
    
    static void outputSolutionLevel1AMR(tuple<vector<RelativisticEulerStateVector>, vector<RelativisticEulerStateVector>, vector<bool>, vector<bool>> solution);
    static void outputSolutionLevel2AMR(tuple<vector<RelativisticEulerStateVector>, vector<RelativisticEulerStateVector>, vector<RelativisticEulerStateVector>,
                                        vector<bool>, vector<bool>, vector<bool>> solution);
    
    static void outputSolution2DLevel1AMR(tuple<vector<vector<RelativisticEulerStateVector> >, vector<vector<RelativisticEulerStateVector> >, vector<vector<bool> >,
                                          vector<vector<bool> >> solution);
    static void outputSolution2DLevel2AMR(tuple<vector<vector<RelativisticEulerStateVector> >, vector<vector<RelativisticEulerStateVector> >,
                                          vector<vector<RelativisticEulerStateVector> >, vector<vector<bool> >, vector<vector<bool> >, vector<vector<bool> >> solution);
};

#endif
