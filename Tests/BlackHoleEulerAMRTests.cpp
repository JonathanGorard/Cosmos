#include "BlackHoleEulerAMRTests.hpp"

void BlackHoleEulerAMRTests::solve1DSphericalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<BlackHoleEulerStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    BlackHoleSpacetime blackHole(0.05, 0.5, 0.0, 0.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        if (i >= (0.1 * cellCount) && i <= (0.9 * cellCount))
        {
            initialCells[i] = BlackHoleEulerStateVector(10.0, 0.0, 0.0, 0.0, 0.1);
        }
        else
        {
            initialCells[i] = BlackHoleEulerStateVector(0.1, 0.0, 0.0, 0.0, 0.1);
        }
    }
    
    if (AMRLevel == 1)
    {
        outputSolutionLevel1AMR(BlackHoleEulerAMR::solveLevel1AMR(initialCells, cellSpacing, 0.95, 0.45, AMRTolerance, order, materialParameters, blackHole),
                                blackHole);
    }
    else
    {
        outputSolutionLevel2AMR(BlackHoleEulerAMR::solveLevel2AMR(initialCells, cellSpacing, 0.95, 0.45, AMRTolerance, order, materialParameters, blackHole),
                                blackHole);
    }
}

void BlackHoleEulerAMRTests::solve1DSpinningSphericalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<BlackHoleEulerStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    BlackHoleSpacetime blackHole(0.1, 0.99, 0.5, 0.0, 0.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        if (i >= (0.1 * cellCount) && i <= (0.9 * cellCount))
        {
            initialCells[i] = BlackHoleEulerStateVector(10.0, 0.0, 0.0, 0.0, 0.1);
        }
        else
        {
            initialCells[i] = BlackHoleEulerStateVector(0.1, 0.0, 0.0, 0.0, 0.1);
        }
    }
    
    if (AMRLevel == 1)
    {
        outputSolutionLevel1AMR(BlackHoleEulerAMR::solveLevel1AMR(initialCells, cellSpacing, 0.95, 0.45, AMRTolerance, order, materialParameters, blackHole),
                                blackHole);
    }
    else
    {
        outputSolutionLevel2AMR(BlackHoleEulerAMR::solveLevel2AMR(initialCells, cellSpacing, 0.95, 0.45, AMRTolerance, order, materialParameters, blackHole),
                                blackHole);
    }
}

void BlackHoleEulerAMRTests::solve2DSphericalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<BlackHoleEulerStateVector> > initialCells(cellCount, vector<BlackHoleEulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    BlackHoleSpacetime blackHole(0.05, 0.5, 0.5, 0.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j< cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.4 * cellCount))
            {
                initialCells[i][j] = BlackHoleEulerStateVector(10.0, 0.0, 0.0, 0.0, 0.1);
            }
            else
            {
                initialCells[i][j] = BlackHoleEulerStateVector(0.1, 0.0, 0.0, 0.0, 0.1);
            }
        }
    }
    
    if (AMRLevel == 1)
    {
        outputSolution2DLevel1AMR(BlackHoleEulerAMR::solve2DLevel1AMR(initialCells, cellSpacing, 0.95, 0.45, AMRTolerance, order, materialParameters, blackHole),
                                  blackHole);
    }
}

void BlackHoleEulerAMRTests::solve2DSpinningSphericalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<BlackHoleEulerStateVector> > initialCells(cellCount, vector<BlackHoleEulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    BlackHoleSpacetime blackHole(0.1, 0.99, 0.5, 0.5, 0.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i -  (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.4 * cellCount))
            {
                initialCells[i][j] = BlackHoleEulerStateVector(10.0, 0.0, 0.0, 0.0, 0.1);
            }
            else
            {
                initialCells[i][j] = BlackHoleEulerStateVector(0.1, 0.0, 0.0, 0.0, 0.1);
            }
        }
    }
    
    if (AMRLevel == 1)
    {
        outputSolution2DLevel1AMR(BlackHoleEulerAMR::solve2DLevel1AMR(initialCells, cellSpacing, 0.95, 0.45, AMRTolerance, order, materialParameters, blackHole),
                                  blackHole);
    }
}

void BlackHoleEulerAMRTests::solve2DSpinningSpheroidalAccretionTest(int cellCount, int order, int subcyclingIterations, double AMRTolerance, int AMRLevel)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<BlackHoleEulerStateVector> > initialCells(cellCount, vector<BlackHoleEulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    BlackHoleSpacetime blackHole(0.1, 0.99, 0.5, 0.5, 0.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount)) * 0.8) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)) * 1.2)) <= (0.4 * cellCount))
            {
                initialCells[i][j] = BlackHoleEulerStateVector(10.0, 0.0, 0.0, 0.0, 0.1);
            }
            else
            {
                initialCells[i][j] = BlackHoleEulerStateVector(0.1, 0.0, 0.0, 0.0, 0.1);
            }
        }
    }
    
    if (AMRLevel == 1)
    {
        outputSolution2DLevel1AMR(BlackHoleEulerAMR::solve2DLevel1AMR(initialCells, cellSpacing, 0.95, 0.45, AMRTolerance, order, materialParameters, blackHole),
                                  blackHole);
    }
}

void BlackHoleEulerAMRTests::outputSolutionLevel1AMR(tuple<vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>, vector<bool>, vector<bool>> solution,
                                                     BlackHoleSpacetime blackHole)
{
    vector<BlackHoleEulerStateVector> coarseSolution = get<0>(solution);
    vector<BlackHoleEulerStateVector> fineSolution = get<1>(solution);
    vector<bool> coarseAMRStructure = get<2>(solution);
    vector<bool> fineAMRStructure = get<3>(solution);
    
    long cellCount = coarseSolution.size();
    double cellSpacing = 1.0 / cellCount;
    
    ofstream densityFile("/Users/jonathangorard/Documents/Cosmos/density.dat");
    ofstream xVelocityFile("/Users/jonathangorard/Documents/Cosmos/xVelocity.dat");
    ofstream yVelocityFile("/Users/jonathangorard/Documents/Cosmos/yVelocity.dat");
    
    for (int i = 0; i < cellCount; i++)
    {
        if (coarseAMRStructure[i])
        {
            if (!blackHole.inExcisionRegion(i * cellSpacing, 0.0, 0.0))
            {
                densityFile << (cellSpacing * i) << " " << coarseSolution[i].getDensity() << endl;
                xVelocityFile << (cellSpacing * i) << " " << coarseSolution[i].getXVelocity() << endl;
                yVelocityFile << (cellSpacing * i) << " " << coarseSolution[i].getYVelocity() << endl;
            }
            else
            {
                densityFile << (cellSpacing * i) << " " << 0.0 << endl;
                xVelocityFile << (cellSpacing * i) << " " << 0.0 << endl;
                yVelocityFile << (cellSpacing * i) << " " << 0.0 << endl;
            }
        }
        else if (fineAMRStructure[(i * 2)] && fineAMRStructure[(i * 2) + 1])
        {
            for (int j = 0; j < 2; j++)
            {
                if (!blackHole.inExcisionRegion((cellSpacing * 0.5) * ((i * 2) + j), 0.0, 0.0))
                {
                    densityFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << fineSolution[(i * 2) + j].getDensity() << endl;
                    xVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << fineSolution[(i * 2) + j].getXVelocity() << endl;
                    yVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << fineSolution[(i * 2) + j].getYVelocity() << endl;
                }
                else
                {
                    densityFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << 0.0 << endl;
                    xVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << 0.0 << endl;
                    yVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << 0.0 << endl;
                }
            }
        }
    }
    
    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
}

void BlackHoleEulerAMRTests::outputSolutionLevel2AMR(tuple<vector<BlackHoleEulerStateVector>, vector<BlackHoleEulerStateVector>,
                                                     vector<BlackHoleEulerStateVector>, vector<bool>, vector<bool>, vector<bool>> solution, BlackHoleSpacetime blackHole)
{
    vector<BlackHoleEulerStateVector> coarseSolution = get<0>(solution);
    vector<BlackHoleEulerStateVector> intermediateSolution = get<1>(solution);
    vector<BlackHoleEulerStateVector> fineSolution = get<2>(solution);
    
    vector<bool> coarseAMRStructure = get<3>(solution);
    vector<bool> intermediateAMRStructure = get<4>(solution);
    vector<bool> fineAMRStructure = get<5>(solution);
    
    long cellCount = coarseSolution.size();
    double cellSpacing = 1.0 / cellCount;
    
    ofstream densityFile("/Users/jonathangorard/Documents/Cosmos/density.dat");
    ofstream xVelocityFile("/Users/jonathangorard/Documents/Cosmos/xVelocity.dat");
    ofstream yVelocityFile("/Users/jonathangorard/Documents/Cosmos/yVelocity.dat");
    
    for (int i = 0; i < cellCount; i++)
    {
        if (coarseAMRStructure[i])
        {
            if (!blackHole.inExcisionRegion(cellSpacing * i, 0.0, 0.0))
            {
                densityFile << (cellSpacing * i) << " " << coarseSolution[i].getDensity() << endl;
                xVelocityFile << (cellSpacing * i) << " " << coarseSolution[i].getXVelocity() << endl;
                yVelocityFile << (cellSpacing * i) << " " << coarseSolution[i].getYVelocity() << endl;
            }
            else
            {
                densityFile << (cellSpacing * i) << " " << 0.0 << endl;
                xVelocityFile << (cellSpacing * i) << " " << 0.0 << endl;
                yVelocityFile << (cellSpacing * i) << " " << 0.0 << endl;
            }
        }
        else if (intermediateAMRStructure[(i * 2)])
        {
            if (!blackHole.inExcisionRegion((cellSpacing * 0.5) * (i * 2), 0.0, 0.0))
            {
                densityFile << ((cellSpacing * 0.5) * (i * 2)) << " " << intermediateSolution[(i * 2)].getDensity() << endl;
                xVelocityFile << ((cellSpacing * 0.5) * (i * 2)) << " " << intermediateSolution[(i * 2)].getXVelocity() << endl;
                yVelocityFile << ((cellSpacing * 0.5) * (i * 2)) << " " << intermediateSolution[(i * 2)].getYVelocity() << endl;
            }
            else
            {
                densityFile << ((cellSpacing * 0.5) * (i * 2)) << " " << 0.0 << endl;
                xVelocityFile << ((cellSpacing * 0.5) * (i * 2)) << " " << 0.0 << endl;
                yVelocityFile << ((cellSpacing * 0.5) * (i * 2)) << " " << 0.0 << endl;
            }
            
            if (intermediateAMRStructure[(i * 2) + 1])
            {
                if (!blackHole.inExcisionRegion((cellSpacing * 0.5) * ((i * 2) + 1), 0.0, 0.0))
                {
                    densityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getDensity() << endl;
                    xVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getXVelocity() << endl;
                    yVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getYVelocity() << endl;
                }
                else
                {
                    densityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << 0.0 << endl;
                    xVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << 0.0 << endl;
                    yVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << 0.0 << endl;
                }
            }
            else
            {
                for (int j = 0; j < 2; j++)
                {
                    if (!blackHole.inExcisionRegion((cellSpacing * 0.25) * ((i * 4) + 2 + j), 0.0, 0.0))
                    {
                        densityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getDensity() << endl;
                        xVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getXVelocity() << endl;
                        yVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getYVelocity() << endl;
                    }
                    else
                    {
                        densityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << 0.0 << endl;
                        xVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << 0.0 << endl;
                        yVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << 0.0 << endl;
                    }
                }
            }
        }
        else if (fineAMRStructure[(i * 4)] && fineAMRStructure[(i * 4) + 1])
        {
            for (int j = 0; j < 2; j++)
            {
                if (!blackHole.inExcisionRegion((cellSpacing * 0.25) * ((i * 4) + j), 0.0, 0.0))
                {
                    densityFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << fineSolution[(i * 4) + j].getDensity() << endl;
                    xVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << fineSolution[(i * 4) + j].getXVelocity() << endl;
                    yVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << fineSolution[(i * 4) + j].getYVelocity() << endl;
                }
                else
                {
                    densityFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << 0.0 << endl;
                    xVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << 0.0 << endl;
                    yVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << 0.0 << endl;
                }
            }
            
            if (fineAMRStructure[(i * 4) + 2] && fineAMRStructure[(i * 4) + 3])
            {
                for (int j = 0; j < 2; j++)
                {
                    if (!blackHole.inExcisionRegion((cellSpacing * 0.25) * ((i * 4) + 2 + j), 0.0, 0.0))
                    {
                        densityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getDensity() << endl;
                        xVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getXVelocity() << endl;
                        yVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getYVelocity() << endl;
                    }
                    else
                    {
                        densityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << 0.0 << endl;
                        xVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << 0.0 << endl;
                        yVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << 0.0 << endl;
                    }
                }
            }
            else
            {
                if (!blackHole.inExcisionRegion((cellSpacing * 0.5) * ((i * 2) + 1), 0.0, 0.0))
                {
                    densityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getDensity() << endl;
                    xVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getXVelocity() << endl;
                    yVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getYVelocity() << endl;
                }
                else
                {
                    densityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << 0.0 << endl;
                    xVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << 0.0 << endl;
                    yVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << 0.0 << endl;
                }
            }
        }
    }
    
    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
}

void BlackHoleEulerAMRTests::outputSolution2DLevel1AMR(tuple<vector<vector<BlackHoleEulerStateVector> >, vector<vector<BlackHoleEulerStateVector> >,
                                                       vector<vector<bool> >, vector<vector<bool> >> solution, BlackHoleSpacetime blackHole)
{
    vector<vector<BlackHoleEulerStateVector> > coarseSolution = get<0>(solution);
    vector<vector<BlackHoleEulerStateVector> > fineSolution = get<1>(solution);
    vector<vector<bool> > coarseAMRStructure = get<2>(solution);
    vector<vector<bool> > fineAMRStructure = get<3>(solution);
    
    long rowCount = coarseSolution.size();
    long columnCount = coarseSolution[0].size();
    double rowCellSpacing = 1.0 / rowCount;
    double columnCellSpacing = 1.0 / columnCount;
    
    ofstream densityFile("/Users/jonathangorard/Documents/Cosmos/density.dat");
    
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (coarseAMRStructure[i][j])
            {
                if (!blackHole.inExcisionRegion(columnCellSpacing * j, rowCellSpacing * i, 0.0))
                {
                    densityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << rowCellSpacing << " " << columnCellSpacing << " "
                    << coarseSolution[i][j].getDensity() << endl;
                }
                else
                {
                    densityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << rowCellSpacing << " " << columnCellSpacing << " " << 0.0 << endl;
                }
            }
            else if (fineAMRStructure[(i * 2)][(j * 2)] && fineAMRStructure[(i * 2)][(j * 2) + 1] && fineAMRStructure[(i * 2) + 1][(j * 2)] &&
                     fineAMRStructure[(i * 2) + 1][(j * 2)])
            {
                for (int k = 0; k < 2; k++)
                {
                    for (int l = 0; l < 2; l++)
                    {
                        if (!blackHole.inExcisionRegion((columnCellSpacing * 0.5) * ((j * 2) + l), (rowCellSpacing * 0.5) * ((i * 2) + k), 0.0))
                        {
                            densityFile << ((rowCellSpacing * 0.5) * ((i * 2) + k)) << " " << ((columnCellSpacing * 0.5) * ((j * 2) + l)) << " " << (rowCellSpacing * 0.5)
                            << " " << (columnCellSpacing * 0.5) << " " << fineSolution[(i * 2) + k][(j * 2) + l].getDensity() << endl;
                        }
                        else
                        {
                            densityFile << ((rowCellSpacing * 0.5) * ((i * 2) + k)) << " " << ((columnCellSpacing * 0.5) * ((j * 2) + l)) << " " << (rowCellSpacing * 0.5)
                            << " " << (columnCellSpacing * 0.5) << " " << 0.0 << endl;
                        }
                    }
                }
            }
        }
    }
    
    densityFile.close();
}
