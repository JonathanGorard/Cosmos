#ifndef BlackHoleEulerAMRTests_hpp
#define BlackHoleEulerAMRTests_hpp

#include "../AMR/BlackHoleEulerAMR.hpp"
#include <fstream>
using namespace std;

class BlackHoleEulerAMRTests
{
public:
    static void solve1DSphericalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel);
    static void solve1DSpinningSphericalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel);
    
    static void solve2DSphericalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel);
    static void solve2DSpinningSphericalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel);
    
    static void solve2DSpheroidalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel);
    static void solve2DSpinningSpheroidalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel);
    
    static void outputSolutionLevel1AMR(tuple<vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<bool>, vector<bool>> solution,
                                        BlackHoleSpacetime blackHole);
    static void outputSolutionLevel2AMR(tuple<vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>,
                                        vector<bool>, vector<bool>, vector<bool>> solution, BlackHoleSpacetime blackHole);
    
    static void outputSolution2DLevel1AMR(tuple<vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >, vector<vector<bool> >,
                                          vector<vector<bool> >> solution, BlackHoleSpacetime spacetime);
};

#endif
