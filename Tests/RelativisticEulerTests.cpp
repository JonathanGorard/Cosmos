#include "RelativisticEulerTests.hpp"

void RelativisticEulerTests::solveMildlyRelativisticShock(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<RelativisticEulerStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = RelativisticEulerStateVector(10.0, 0.0, 0.0, 0.0, 13.3);
        }
        else
        {
            initialCells[i] = RelativisticEulerStateVector(1.0, 0.0, 0.0, 0.0, 0.0);
        }
    }
    
    if (order == 1)
    {
        outputSolution(RelativisticEulerFirstOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.4, materialParameters));
    }
    else
    {
        outputSolution(RelativisticEulerSecondOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.01, 0.0, 1, materialParameters));
    }
}

void RelativisticEulerTests::solveStronglyRelativisticBlast(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<RelativisticEulerStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = RelativisticEulerStateVector(1.0, 0.0, 0.0, 0.0, 1000.0);
        }
        else
        {
            initialCells[i] = RelativisticEulerStateVector(1.0, 0.0, 0.0, 0.0, 0.01);
        }
    }
    
    if (order == 1)
    {
        outputSolution(RelativisticEulerFirstOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.4, materialParameters));
    }
    else
    {
        outputSolution(RelativisticEulerSecondOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.4, 0.0, 1, materialParameters));
    }
}

void RelativisticEulerTests::solvePerturbedDensityTest(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<RelativisticEulerStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = RelativisticEulerStateVector(5.0, 0.0, 0.0, 0.0, 50.0);
        }
        else
        {
            initialCells[i] = RelativisticEulerStateVector(2.0 + 0.3 * sin(50.0 * (cellSpacing * i)), 0.0, 0.0, 0.0, 5.0);
        }
    }
    
    if (order == 1)
    {
        outputSolution(RelativisticEulerFirstOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.35, materialParameters));
    }
    else
    {
        outputSolution(RelativisticEulerSecondOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.35, 0.0, 1, materialParameters));
    }
}

void RelativisticEulerTests::solve2DMildlyRelativisticShock(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<RelativisticEulerStateVector> > initialCells(cellCount, vector<RelativisticEulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = RelativisticEulerStateVector(10.0, 0.0, 0.0, 0.0, 13.3);
            }
            else
            {
                initialCells[i][j] = RelativisticEulerStateVector(1.0, 0.0, 0.0, 0.0, 0.0);
            }
        }
    }
    
    if (order == 1)
    {
        outputSolution2D(RelativisticEulerFirstOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.2, materialParameters));
    }
    else
    {
        outputSolution2D(RelativisticEulerSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.2, 0.0, 1, materialParameters));
    }
}

void RelativisticEulerTests::solve2DStronglyRelativisticBlast(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<RelativisticEulerStateVector> > initialCells(cellCount, vector<RelativisticEulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = RelativisticEulerStateVector(1.0, 0.0, 0.0, 0.0, 1000.0);
            }
            else
            {
                initialCells[i][j] = RelativisticEulerStateVector(1.0, 0.0, 0.0, 0.0, 0.01);
            }
        }
    }
    
    if (order == 1)
    {
        outputSolution2D(RelativisticEulerFirstOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.15, materialParameters));
    }
    else
    {
        outputSolution2D(RelativisticEulerSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.15, 0.0, 1, materialParameters));
    }
}

void RelativisticEulerTests::solve2DPerturbedDensityTest(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<RelativisticEulerStateVector> > initialCells(cellCount, vector<RelativisticEulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(5.0 / 3.0);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = RelativisticEulerStateVector(5.0, 0.0, 0.0, 0.0, 50.0);
            }
            else
            {
                initialCells[i][j] = RelativisticEulerStateVector(2.0 + 0.3 * sin(50.0 * sqrt(((cellSpacing * (i - (0.5 * cellCount))) *
                                                                                               (cellSpacing * (i - (0.5 * cellCount)))) +
                                                                                              ((cellSpacing * (j - (0.5 * cellCount))) *
                                                                                               (cellSpacing * (j - (0.5 * cellCount)))))), 0.0, 0.0, 0.0, 5.0);
            }
        }
    }
    
    if (order == 1)
    {
        outputSolution2D(RelativisticEulerFirstOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.175, materialParameters));
    }
    else
    {
        outputSolution2D(RelativisticEulerSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.175, 0.0, 1, materialParameters));
    }
}

void RelativisticEulerTests::outputSolution(vector<RelativisticEulerStateVector> solution)
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
        densityFile << (cellSpacing * i) << " " << solution[i].getDensity() << endl;
        xVelocityFile << (cellSpacing * i) << " " << solution[i].getXVelocity() << endl;
        yVelocityFile << (cellSpacing * i) << " " << solution[i].getYVelocity() << endl;
        zVelocityFile << (cellSpacing * i) << " " << solution[i].getZVelocity() << endl;
        pressureFile << (cellSpacing * i) << " " << solution[i].getPressure() << endl;
    }
    
    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
    zVelocityFile.close();
    pressureFile.close();
}

void RelativisticEulerTests::outputSolution2D(vector<vector<RelativisticEulerStateVector> > solution)
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
            densityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << solution[i][j].getDensity() << endl;
            xVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << solution[i][j].getXVelocity() << endl;
            yVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << solution[i][j].getYVelocity() << endl;
            zVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << solution[i][j].getZVelocity() << endl;
            pressureFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << solution[i][j].getPressure() << endl;
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
