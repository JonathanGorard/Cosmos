#include "EulerTests.hpp"

void EulerTests::solveToroTest1(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<EulerStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters(1.4);
    
    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0);
        }
        else
        {
            initialCells[i] = EulerStateVector(0.125, 0.0, 0.0, 0.0, 0.1);
        }
    }
    
    if (order == 1)
    {
        outputSolution(EulerFirstOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.25, materialParameters));
    }
    else
    {
        outputSolution(EulerSecondOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.25, 0.0, 1, materialParameters));
    }
}

void EulerTests::solveToroTest2(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<EulerStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters(1.4);
    
    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = EulerStateVector(1.0, -2.0, 0.0, 0.0, 0.4);
        }
        else
        {
            initialCells[i] = EulerStateVector(1.0, 2.0, 0.0, 0.0, 0.4);
        }
    }
    
    if (order == 1)
    {
        outputSolution(EulerFirstOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.15, materialParameters));
    }
    else
    {
        outputSolution(EulerSecondOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.15, 0.0, 1, materialParameters));
    }
}

void EulerTests::solveToroTest3(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<EulerStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters(1.4);
    
    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1000.0);
        }
        else
        {
            initialCells[i] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 0.01);
        }
    }
    
    if (order == 1)
    {
        outputSolution(EulerFirstOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.012, materialParameters));
    }
    else
    {
        outputSolution(EulerSecondOrderSolver::solve(initialCells, cellSpacing, 0.95, 0.012, 0.0, 1, materialParameters));
    }
}

void EulerTests::solve2DToroTest1(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<EulerStateVector> > initialCells(cellCount, vector<EulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(1.4);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0);
            }
            else
            {
                initialCells[i][j] = EulerStateVector(0.125, 0.0, 0.0, 0.0, 0.1);
            }
        }
    }
    
    if (order == 1)
    {
        outputSolution2D(EulerFirstOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.125, materialParameters));
    }
    else
    {
        outputSolution2D(EulerSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.125, 0.0, 1, materialParameters));
    }
}

void EulerTests::solve2DToroTest2(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<EulerStateVector> > initialCells(cellCount, vector<EulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(1.4);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = EulerStateVector(1.0, -2.0 * sin(atan2(j - (0.5 * cellCount), i - (0.5 * cellCount))),
                                                      -2.0 * cos(atan2(j - (0.5 * cellCount), i - (0.5 * cellCount))), 0.0, 0.4);
            }
            else
            {
                initialCells[i][j] = EulerStateVector(1.0, 2.0 * sin(atan2(j - (0.5 * cellCount), i - (0.5 * cellCount))),
                                                      2.0 * cos(atan2(j - (0.5 * cellCount), i - (0.5 * cellCount))), 0.0, 0.4);
            }
        }
    }
    
    if (order == 1)
    {
        outputSolution2D(EulerFirstOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.075, materialParameters));
    }
    else
    {
        outputSolution2D(EulerSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.075, 0.0, 1, materialParameters));
    }
}

void EulerTests::solve2DToroTest3(int cellCount, int order)
{
    double cellSpacing = 1.0 / cellCount;
    
    vector<vector<EulerStateVector> > initialCells(cellCount, vector<EulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(1.4);
    
    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1000.0);
            }
            else
            {
                initialCells[i][j] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 0.01);
            }
        }
    }
    
    if (order == 1)
    {
        outputSolution2D(EulerFirstOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.006, materialParameters));
    }
    else
    {
        outputSolution2D(EulerSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.95, 0.006, 0.0, 1, materialParameters));
    }
}


void EulerTests::outputSolution(vector<EulerStateVector> solution)
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

void EulerTests::outputSolution2D(vector<vector<EulerStateVector> > solution)
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
