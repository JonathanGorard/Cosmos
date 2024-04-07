#include "BlackHoleEulerTests.hpp"

void BlackHoleEulerTests::solve1DSphericalAccretionTest(int cellCount, int order, int subcyclingIterations)
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
    
    if (order == 1)
    {
        outputSolution(BlackHoleEulerFirstOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.45, subcyclingIterations, materialParameters, blackHole), blackHole);
    }
    else
    {
        outputSolution(BlackHoleEulerSecondOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.45, 0.0, 1, subcyclingIterations, materialParameters, blackHole),
                       blackHole);
    }
}

void BlackHoleEulerTests::solve1DSpinningSphericalAccretionTest(int cellCount, int order, int subcyclingIterations)
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
    
    if (order == 1)
    {
        outputSolution(BlackHoleEulerFirstOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.45, subcyclingIterations, materialParameters, blackHole), blackHole);
    }
    else
    {
        outputSolution(BlackHoleEulerSecondOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.45, 0.0, 1, subcyclingIterations, materialParameters, blackHole),
                       blackHole);
    }
}

void BlackHoleEulerTests::solve2DSphericalAccretionTest(int cellCount, int order, int subcyclingIterations)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<BlackHoleEulerStateVector> > initialCells(cellCount, vector<BlackHoleEulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    BlackHoleSpacetime blackHole(0.05, 0.5, 0.5, 0.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= 0.4 * cellCount)
            {
                initialCells[i][j] = BlackHoleEulerStateVector(10.0, 0.0, 0.0, 0.0, 0.1);
            }
            else
            {
                initialCells[i][j] = BlackHoleEulerStateVector(0.1, 0.0, 0.0, 0.0, 0.1);
            }
        }
    }
    
    if (order == 1)
    {
        outputSolution2D(BlackHoleEulerFirstOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.45, subcyclingIterations, materialParameters, blackHole), blackHole);
    }
    else
    {
        outputSolution2D(BlackHoleEulerSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.45, 0.0, 1, subcyclingIterations, materialParameters, blackHole),
                         blackHole);
    }
}

void BlackHoleEulerTests::solve2DSpinningSphericalAccretionTest(int cellCount, int order, int subcyclingIterations)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<BlackHoleEulerStateVector> > initialCells(cellCount, vector<BlackHoleEulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    BlackHoleSpacetime blackHole(0.1, 0.99, 0.5, 0.5, 0.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= 0.4 * cellCount)
            {
                initialCells[i][j] = BlackHoleEulerStateVector(10.0, 0.0, 0.0, 0.0, 0.1);
            }
            else
            {
                initialCells[i][j] = BlackHoleEulerStateVector(0.1, 0.0, 0.0, 0.0, 0.1);
            }
        }
    }
    
    if (order == 1)
    {
        outputSolution2D(BlackHoleEulerFirstOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.45, subcyclingIterations, materialParameters, blackHole), blackHole);
    }
    else
    {
        outputSolution2D(BlackHoleEulerSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.45, 0.0, 1, subcyclingIterations, materialParameters, blackHole),
                         blackHole);
    }
}

void BlackHoleEulerTests::solve2DSpheroidalAccretionTest(int cellCount, int order, int subcyclingIterations)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<BlackHoleEulerStateVector> > initialCells(cellCount, vector<BlackHoleEulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    BlackHoleSpacetime blackHole(0.05, 0.5, 0.5, 0.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount)) * 0.8) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)) * 1.2)) <= 0.4 * cellCount)
            {
                initialCells[i][j] = BlackHoleEulerStateVector(10.0, 0.0, 0.0, 0.0, 0.1);
            }
            else
            {
                initialCells[i][j] = BlackHoleEulerStateVector(0.1, 0.0, 0.0, 0.0, 0.1);
            }
        }
    }
    
    if (order == 1)
    {
        outputSolution2D(BlackHoleEulerFirstOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.45, subcyclingIterations, materialParameters, blackHole), blackHole);
    }
    else
    {
        outputSolution2D(BlackHoleEulerSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.45, 0.0, 1, subcyclingIterations, materialParameters, blackHole),
                         blackHole);
    }
}

void BlackHoleEulerTests::solve2DSpinningSpheroidalAccretionTest(int cellCount, int order, int subcyclingIterations)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<BlackHoleEulerStateVector> > initialCells(cellCount, vector<BlackHoleEulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    BlackHoleSpacetime blackHole(0.1, 0.99, 0.5, 0.5, 0.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount)) * 0.8) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)) * 1.2)) <= 0.4 * cellCount)
            {
                initialCells[i][j] = BlackHoleEulerStateVector(10.0, 0.0, 0.0, 0.0, 0.1);
            }
            else
            {
                initialCells[i][j] = BlackHoleEulerStateVector(0.1, 0.0, 0.0, 0.0, 0.1);
            }
        }
    }
    
    if (order == 1)
    {
        outputSolution2D(BlackHoleEulerFirstOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.45, subcyclingIterations, materialParameters, blackHole), blackHole);
    }
    else
    {
        outputSolution2D(BlackHoleEulerSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.45, 0.0, 1, subcyclingIterations, materialParameters, blackHole),
                         blackHole);
    }
}

void BlackHoleEulerTests::outputSolution(vector<BlackHoleEulerStateVector> solution, BlackHoleSpacetime blackHole)
{
    long cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;
    
    ofstream densityFile("/Users/jonathangorard/Documents/Cosmos/density.dat");
    ofstream xVelocityFile("/Users/jonathangorard/Documents/Cosmos/xVelocity.dat");
    ofstream yVelocityFile("/Users/jonathangorard/Documents/Cosmos/yVelocity.dat");
    ofstream zVelocityFile("/Users/jonathangorard/Documents/Cosmos/zVelocity.dat");
    ofstream pressureFile("/Users/jonathangorard/Documents/Cosmos/pressure.dat");
    
    for (int i = 0; i < cellCount; i++)
    {
        if (!blackHole.inExcisionRegion(i * cellSpacing, 0.0, 0.0))
        {
            densityFile << (cellSpacing * i) << " " << solution[i].getDensity() << endl;
            xVelocityFile << (cellSpacing * i) << " " << solution[i].getXVelocity() << endl;
            yVelocityFile << (cellSpacing * i) << " " << solution[i].getYVelocity() << endl;
            zVelocityFile << (cellSpacing * i) << " " << solution[i].getZVelocity() << endl;
            pressureFile << (cellSpacing * i) << " " << solution[i].getPressure() << endl;
        }
        else
        {
            densityFile << (cellSpacing * i) << " " << 0.0 << endl;
            xVelocityFile << (cellSpacing * i) << " " << 0.0 << endl;
            yVelocityFile << (cellSpacing * i) << " " << 0.0 << endl;
            zVelocityFile << (cellSpacing * i) << " " << 0.0 << endl;
            pressureFile << (cellSpacing * i) << " " << 0.0 << endl;
        }
    }
    
    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
    zVelocityFile.close();
    pressureFile.close();
}

void BlackHoleEulerTests::outputSolution2D(vector<vector<BlackHoleEulerStateVector> > solution, BlackHoleSpacetime blackHole)
{
    long rowCount = solution.size();
    long columnCount = solution[0].size();
    double rowCellSpacing = 1.0 / rowCount;
    double columnCellSpacing = 1.0 / columnCount;
    
    ofstream densityFile("/Users/jonathangorard/Documents/Cosmos/density.dat");
    ofstream xVelocityFile("/Users/jonathangorard/Documents/Cosmos/xVelocity.dat");
    ofstream yVelocityFile("/Users/jonathangorard/Documents/Cosmos/yVelocity.dat");
    ofstream zVelocityFile("/Users/jonathangorard/Documents/Cosmos/zVelocity.dat");
    ofstream pressureFile("/Users/jonathangorard/Documents/Cosmos/pressure.dat");
    
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            if (!blackHole.inExcisionRegion(j * columnCellSpacing, i * rowCellSpacing, 0.0))
            {
                densityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << solution[i][j].getDensity() << endl;
                xVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << solution[i][j].getXVelocity() << endl;
                yVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << solution[i][j].getYVelocity() << endl;
                zVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << solution[i][j].getZVelocity() << endl;
                pressureFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << solution[i][j].getPressure() << endl;
            }
            else
            {
                densityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << 0.0 << endl;
                xVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << 0.0 << endl;
                yVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << 0.0 << endl;
                zVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << 0.0 << endl;
                pressureFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << 0.0 << endl;
            }
        }
        
        densityFile << endl;
        xVelocityFile << endl;
        yVelocityFile << endl;
        zVelocityFile << endl;
        pressureFile << endl;
    }
    
    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
    zVelocityFile.close();
    pressureFile.close();
}
