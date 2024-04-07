#ifndef BlackHoleEulerTests_hpp
#define BlackHoleEulerTests_hpp

#include "../Solvers/BlackHoleEulerFirstOrderSolver.hpp"
#include "../Solvers/BlackHoleEulerSecondOrderSolver.hpp"

#include <fstream>
using namespace std;

class BlackHoleEulerTests
{
public:
    static void solve1DSphericalAccretionTest(int cellCount, int order, int subcyclingIterations);
    static void solve1DSpinningSphericalAccretionTest(int cellCount, int order, int subcyclingIterations);
    
    static void solve2DSphericalAccretionTest(int cellCount, int order, int subcyclingIterations);
    static void solve2DSpinningSphericalAccretionTest(int cellCount, int order, int subcyclingIterations);
    
    static void solve2DSpheroidalAccretionTest(int cellCount, int order, int subcyclingIterations);
    static void solve2DSpinningSpheroidalAccretionTest(int cellCount, int order, int subcyclingIterations);
    
    static void outputSolution(vector<BlackHoleEulerStateVector> solution, BlackHoleSpacetime blackHole);
    static void outputSolution2D(vector<vector<BlackHoleEulerStateVector> > solution, BlackHoleSpacetime blackHole);
};

#endif
