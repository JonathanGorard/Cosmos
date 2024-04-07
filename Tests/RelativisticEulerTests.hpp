#ifndef RelativisticEulerTests_hpp
#define RelativisticEulerTests_hpp

#include "../Solvers/RelativisticEulerFirstOrderSolver.hpp"
#include "../Solvers/RelativisticEulerSecondOrderSolver.hpp"
#include <fstream>
using namespace std;

class RelativisticEulerTests
{
public:
    static void solveMildlyRelativisticShock(int cellCount, int order);
    static void solveStronglyRelativisticBlast(int cellCount, int order);
    static void solvePerturbedDensityTest(int cellCount, int order);
    
    static void solve2DMildlyRelativisticShock(int cellCount, int order);
    static void solve2DStronglyRelativisticBlast(int cellCount, int order);
    static void solve2DPerturbedDensityTest(int cellCount, int order);
    
    static void outputSolution(vector<RelativisticEulerStateVector> solution);
    static void outputSolution2D(vector<vector<RelativisticEulerStateVector> > solution);
};

#endif
