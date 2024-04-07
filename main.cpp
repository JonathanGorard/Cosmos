#include <iostream>
#include <fstream>
#include "Tests/EulerTests.hpp"
#include "Tests/EulerAMRTests.hpp"
#include "Tests/RelativisticEulerTests.hpp"
#include "Tests/RelativisticEulerAMRTests.hpp"
#include "Tests/BlackHoleEulerTests.hpp"
#include "Tests/BlackHoleEulerAMRTests.hpp"

using namespace std;

int main()
{
    //EulerTests::solveToroTest1(800, 2);
    //EulerTests::solveToroTest2(800, 2);
    //EulerTests::solveToroTest3(2000, 2);
    //EulerTests::solve2DToroTest1(200, 2);
    //EulerTests::solve2DToroTest2(200, 2);
    //EulerTests::solve2DToroTest3(200, 2);
    //EulerAMRTests::solveToroTest1AMR(200, 2, 1.0, 2);
    //EulerAMRTests::solve2DToroTest1AMR(100, 2, 1.0, 2);
    
    RelativisticEulerTests::solveMildlyRelativisticShock(800, 2);
    //RelativisticEulerTests::solveStronglyRelativisticBlast(400, 2);
    //RelativisticEulerTests::solvePerturbedDensityTest(400, 2);
    //RelativisticEulerTests::solve2DMildlyRelativisticShock(400, 2);
    //RelativisticEulerTests::solve2DStronglyRelativisticBlast(200, 2);
    //RelativisticEulerTests::solve2DPerturbedDensityTest(400, 2);
    //RelativisticEulerAMRTests::solveMildlyRelativisticShock(100, 2, 5.0, 2);
    //RelativisticEulerAMRTests::solveStronglyRelativisticBlast(100, 2, 100.0, 2);
    //RelativisticEulerAMRTests::solvePerturbedDensityTest(200, 2, 10.0, 2);
    //RelativisticEulerAMRTests::solve2DMildlyRelativisticShock(100, 2, 10.0, 2);
    //RelativisticEulerAMRTests::solve2DStronglyRelativisticBlast(100, 1, 1000.0, 2);
    //RelativisticEulerAMRTests::solve2DPerturbedDensityTest(100, 2, 20.0, 2);
    
    //BlackHoleEulerTests::solve1DSphericalAccretionTest(400, 2, 0);
    //BlackHoleEulerTests::solve1DSpinningSphericalAccretionTest(400, 2, 0);
    //BlackHoleEulerTests::solve2DSphericalAccretionTest(100, 2, 0);
    //BlackHoleEulerTests::solve2DSpinningSphericalAccretionTest(100, 2, 0);
    //BlackHoleEulerTests::solve2DSpheroidalAccretionTest(100, 2, 0);
    //BlackHoleEulerTests::solve2DSpinningSpheroidalAccretionTest(200, 2, 0);
    //BlackHoleEulerAMRTests::solve1DSphericalAccretionTest(200, 2, 0, 5.0, 2);
    //BlackHoleEulerAMRTests::solve1DSpinningSphericalAccretionTest(200, 2, 0, 5.0, 2);
    //BlackHoleEulerAMRTests::solve2DSphericalAccretionTest(100, 2, 0, 5.0, 1);
    //BlackHoleEulerAMRTests::solve2DSpinningSphericalAccretionTest(100, 2, 0, 5.0, 1);
    //BlackHoleEulerAMRTests::solve2DSpinningSpheroidalAccretionTest(100, 2, 0, 5.0, 1);
    
    return 0;
}
