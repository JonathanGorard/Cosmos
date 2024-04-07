#include "AMRHelper.hpp"

void AMRHelper::outputCoarseStatus(int coarseCurrentIteration, double coarseCurrentTime, double timeStep)
{
    cout << "Coarse Grid Iteration = " << coarseCurrentIteration << "; Time = " << coarseCurrentTime << "; Timestep = " << timeStep << endl;
}

void AMRHelper::outputIntermediateStatus(int intermediateCurrentIteration, double intermediateCurrentTime, double intermediateTimeStep)
{
    cout << "     Refinement Level 1 Iteration = " << intermediateCurrentIteration << "; Time = " << intermediateCurrentTime
    << "; Timestep = " << intermediateTimeStep << endl;
}

void AMRHelper::outputFineStatus(int fineCurrentIteration, double fineCurrentTime, double fineTimeStep)
{
    cout << "          Refinement Level 2 Iteration = " << fineCurrentIteration << "; Time = " << fineCurrentTime << "; Timestep = " << fineTimeStep << endl;
}
