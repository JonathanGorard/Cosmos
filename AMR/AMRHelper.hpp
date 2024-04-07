#ifndef AMRHelper_hpp
#define AMRHelper_hpp

#include "../Solvers/EulerFirstOrderSolver.hpp"

class AMRHelper
{
public:
    static void outputCoarseStatus(int coarseCurrentIteration, double coarseCurrentTime, double coarseTimeStep);
    static void outputIntermediateStatus(int intermediateCurrentIteration, double intermediateCurrentTime, double intermediateTimeStep);
    static void outputFineStatus(int fineCurrentIteration, double fineCurrentTime, double fineTimeStep);
};

#endif
