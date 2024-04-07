#ifndef EulerTests_hpp
#define EulerTests_hpp

#include "../Solvers/EulerFirstOrderSolver.hpp"
#include "../Solvers/EulerSecondOrderSolver.hpp"
#include <fstream>
using namespace std;

class EulerTests
{
public:
    static void solveToroTest1(int cellCount, int order);
    static void solveToroTest2(int cellCount, int order);
    static void solveToroTest3(int cellCount, int order);
    static void solveToroTest4(int cellCount, int order);
    static void solveToroTest5(int cellCount, int order);
    
    static void solve2DToroTest1(int cellCount, int order);
    static void solve2DToroTest2(int cellCount, int order);
    static void solve2DToroTest3(int cellCount, int order);
    static void solve2DToroTest4(int cellCount, int order);
    static void solve2DToroTest5(int cellCount, int order);
    
    static void outputSolution(vector<EulerStateVector> solution);
    static void outputSolution2D(vector<vector<EulerStateVector> > solution);
};

#endif
