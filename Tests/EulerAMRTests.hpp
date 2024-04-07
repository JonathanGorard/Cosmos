#ifndef EulerAMRTests_hpp
#define EulerAMRTests_hpp

#include "../AMR/EulerAMR.hpp"
#include <fstream>
using namespace std;

class EulerAMRTests
{
public:
    static void solveToroTest1AMR(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solveToroTest2AMR(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solveToroTest3AMR(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solveToroTest4AMR(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solveToroTest5AMR(int cellCount, int order, double AMRTolerance, int AMRLevel);
    
    static void solve2DToroTest1AMR(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solve2DToroTest2AMR(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solve2DToroTest3AMR(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solve2DToroTest4AMR(int cellCount, int order, double AMRTolerance, int AMRLevel);
    static void solve2DToroTest5AMR(int cellCount, int order, double AMRTolerance, int AMRLevel);
    
    static void outputSolutionLevel1AMR(tuple<vector<EulerStateVector>, vector<EulerStateVector>, vector<bool>, vector<bool>> solution);
    static void outputSolutionLevel2AMR(tuple<vector<EulerStateVector>, vector<EulerStateVector>, vector<EulerStateVector>, vector<bool>, vector<bool>,
                                        vector<bool>> solution);
    
    static void outputSolution2DLevel1AMR(tuple<vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<bool> >, vector<vector<bool> >> solution);
    static void outputSolution2DLevel2AMR(tuple<vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >,
                                          vector<vector<bool> >, vector<vector<bool> >, vector<vector<bool> >> solution);
};

#endif
