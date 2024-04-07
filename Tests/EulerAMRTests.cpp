#include "EulerAMRTests.hpp"

void EulerAMRTests::solveToroTest1AMR(int cellCount, int order, double AMRTolerance, int AMRLevel)
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
    
    if (AMRLevel == 1)
    {
        outputSolutionLevel1AMR(EulerAMR::solveLevel1AMR(initialCells, cellSpacing, 0.95, 0.25, AMRTolerance, order, materialParameters));
    }
    else
    {
        outputSolutionLevel2AMR(EulerAMR::solveLevel2AMR(initialCells, cellSpacing, 0.95, 0.25, AMRTolerance, order, materialParameters));
    }
}

void EulerAMRTests::solve2DToroTest1AMR(int cellCount, int order, double AMRTolerance, int AMRLevel)
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
    
    if (AMRLevel == 1)
    {
        outputSolution2DLevel1AMR(EulerAMR::solve2DLevel1AMR(initialCells, cellSpacing, 0.95, 0.125, AMRTolerance, order, materialParameters));
    }
    else
    {
        outputSolution2DLevel2AMR(EulerAMR::solve2DLevel2AMR(initialCells, cellSpacing, 0.95, 0.125, AMRTolerance, order, materialParameters));
    }
}

void EulerAMRTests::outputSolutionLevel1AMR(tuple<vector<EulerStateVector>, vector<EulerStateVector>, vector<bool>, vector<bool>> solution)
{
    vector<EulerStateVector> coarseSolution = get<0>(solution);
    vector<EulerStateVector> fineSolution = get<1>(solution);
    vector<bool> coarseAMRStructure = get<2>(solution);
    vector<bool> fineAMRStructure = get<3>(solution);
    
    long cellCount = coarseSolution.size();
    double cellSpacing = 1.0 / cellCount;
    
    ofstream densityFile("/Users/jonathangorard/Documents/Cosmos/density.dat");
    ofstream xVelocityFile("/Users/jonathangorard/Documents/Cosmos/xVelocity.dat");
    ofstream yVelocityFile("/Users/jonathangorard/Documents/Cosmos/yVelocity.dat");
    ofstream zVelocityFile("/Users/jonathangorard/Documents/Cosmos/zVelocity.dat");
    ofstream pressureFile("/Users/jonathangorard/Documents/Cosmos/pressure.dat");
    
    for (int i = 0; i < cellCount; i++)
    {
        if (coarseAMRStructure[i])
        {
            densityFile << (cellSpacing * i) << " " << coarseSolution[i].getDensity() << endl;
            xVelocityFile << (cellSpacing * i) << " " << coarseSolution[i].getXVelocity() << endl;
            yVelocityFile << (cellSpacing * i) << " " << coarseSolution[i].getYVelocity() << endl;
            zVelocityFile << (cellSpacing * i) << " " << coarseSolution[i].getZVelocity() << endl;
            pressureFile << (cellSpacing * i) << " " << coarseSolution[i].getPressure() << endl;
        }
        else if (fineAMRStructure[(i * 2)] && fineAMRStructure[(i * 2) + 1])
        {
            for (int j = 0; j < 2; j++)
            {
                densityFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << fineSolution[(i * 2) + j].getDensity() << endl;
                xVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << fineSolution[(i * 2) + j].getXVelocity() << endl;
                yVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << fineSolution[(i * 2) + j].getYVelocity() << endl;
                zVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << fineSolution[(i * 2) + j].getZVelocity() << endl;
                pressureFile << ((cellSpacing * 0.5) * ((i * 2) + j)) << " " << fineSolution[(i * 2) + j].getPressure() << endl;
            }
        }
    }
    
    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
    zVelocityFile.close();
    pressureFile.close();
}

void EulerAMRTests::outputSolutionLevel2AMR(tuple<vector<EulerStateVector>, vector<EulerStateVector>, vector<EulerStateVector>, vector<bool>, vector<bool>,
                                            vector<bool>> solution)
{
    vector<EulerStateVector> coarseSolution = get<0>(solution);
    vector<EulerStateVector> intermediateSolution = get<1>(solution);
    vector<EulerStateVector> fineSolution = get<2>(solution);
    
    vector<bool> coarseAMRStructure = get<3>(solution);
    vector<bool> intermediateAMRStructure = get<4>(solution);
    vector<bool> fineAMRStructure = get<5>(solution);
    
    long cellCount = coarseSolution.size();
    double cellSpacing = 1.0 / cellCount;
    
    ofstream densityFile("/Users/jonathangorard/Documents/Cosmos/density.dat");
    ofstream xVelocityFile("/Users/jonathangorard/Documents/Cosmos/xVelocity.dat");
    ofstream yVelocityFile("/Users/jonathangorard/Documents/Cosmos/yVelocity.dat");
    ofstream zVelocityFile("/Users/jonathangorard/Documents/Cosmos/zVelocity.dat");
    ofstream pressureFile("/Users/jonathangorard/Documents/Cosmos/pressure.dat");
    
    for (int i = 0; i < cellCount; i++)
    {
        if (coarseAMRStructure[i])
        {
            densityFile << (cellSpacing * i) << " " << coarseSolution[i].getDensity() << endl;
            xVelocityFile << (cellSpacing * i) << " " << coarseSolution[i].getXVelocity() << endl;
            yVelocityFile << (cellSpacing * i) << " " << coarseSolution[i].getYVelocity() << endl;
            zVelocityFile << (cellSpacing * i) << " " << coarseSolution[i].getZVelocity() << endl;
            pressureFile << (cellSpacing * i) << " " << coarseSolution[i].getPressure() << endl;
        }
        else if (intermediateAMRStructure[(i * 2)])
        {
            densityFile << ((cellSpacing * 0.5) * (i * 2)) << " " << intermediateSolution[(i * 2)].getDensity() << endl;
            xVelocityFile << ((cellSpacing * 0.5) * (i * 2)) << " " << intermediateSolution[(i * 2)].getXVelocity() << endl;
            yVelocityFile << ((cellSpacing * 0.5) * (i * 2)) << " " << intermediateSolution[(i * 2)].getYVelocity() << endl;
            zVelocityFile << ((cellSpacing * 0.5) * (i * 2)) << " " << intermediateSolution[(i * 2)].getZVelocity() << endl;
            pressureFile << ((cellSpacing * 0.5) * (i * 2)) << " " << intermediateSolution[(i * 2)].getPressure() << endl;
            
            if (intermediateAMRStructure[(i * 2) + 1])
            {
                densityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getDensity() << endl;
                xVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getXVelocity() << endl;
                yVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getYVelocity() << endl;
                zVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getZVelocity() << endl;
                pressureFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getPressure() << endl;
            }
            else
            {
                for (int j = 0; j < 2; j++)
                {
                    densityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getDensity() << endl;
                    xVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getXVelocity() << endl;
                    yVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getYVelocity() << endl;
                    zVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getZVelocity() << endl;
                    pressureFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getPressure() << endl;
                }
            }
        }
        else if (fineAMRStructure[(i * 4)] && fineAMRStructure[(i * 4) + 1])
        {
            for (int j = 0; j < 2; j++)
            {
                densityFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << fineSolution[(i * 4) + j].getDensity() << endl;
                xVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << fineSolution[(i * 4) + j].getXVelocity() << endl;
                yVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << fineSolution[(i * 4) + j].getYVelocity() << endl;
                zVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << fineSolution[(i * 4) + j].getZVelocity() << endl;
                pressureFile << ((cellSpacing * 0.25) * ((i * 4) + j)) << " " << fineSolution[(i * 4) + j].getPressure() << endl;
            }
            
            if (fineAMRStructure[(i * 4) + 2] && fineAMRStructure[(i * 4) + 3])
            {
                for (int j = 0; j < 2; j++)
                {
                    densityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getDensity() << endl;
                    xVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getXVelocity() << endl;
                    yVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getYVelocity() << endl;
                    zVelocityFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getZVelocity() << endl;
                    pressureFile << ((cellSpacing * 0.25) * ((i * 4) + 2 + j)) << " " << fineSolution[(i * 4) + 2 + j].getPressure() << endl;
                }
            }
            else
            {
                densityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getDensity() << endl;
                xVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getXVelocity() << endl;
                yVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getYVelocity() << endl;
                zVelocityFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getZVelocity() << endl;
                pressureFile << ((cellSpacing * 0.5) * ((i * 2) + 1)) << " " << intermediateSolution[(i * 2) + 1].getPressure() << endl;
            }
        }
    }
    
    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
    zVelocityFile.close();
    pressureFile.close();
}

void EulerAMRTests::outputSolution2DLevel1AMR(tuple<vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<bool> >,
                                              vector<vector<bool> >> solution)
{
    vector<vector<EulerStateVector> > coarseSolution = get<0>(solution);
    vector<vector<EulerStateVector> > fineSolution = get<1>(solution);
    vector<vector<bool> > coarseAMRStructure = get<2>(solution);
    vector<vector<bool> > fineAMRStructure = get<3>(solution);
    
    long rowCount = coarseSolution.size();
    long columnCount = coarseSolution[0].size();
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
            if (coarseAMRStructure[i][j])
            {
                densityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << rowCellSpacing << " " << columnCellSpacing << " "
                << coarseSolution[i][j].getDensity() << endl;
                xVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << rowCellSpacing << " " << columnCellSpacing << " "
                << coarseSolution[i][j].getXVelocity() << endl;
                yVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << rowCellSpacing << " " << columnCellSpacing << " "
                << coarseSolution[i][j].getYVelocity() << endl;
                zVelocityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << rowCellSpacing << " " << columnCellSpacing << " "
                << coarseSolution[i][j].getZVelocity() << endl;
                pressureFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << rowCellSpacing << " " << columnCellSpacing << " "
                << coarseSolution[i][j].getPressure() << endl;
            }
            else if (fineAMRStructure[(i * 2)][(j * 2)] && fineAMRStructure[(i * 2)][(j * 2) + 1] && fineAMRStructure[(i * 2) + 1][(j * 2)] &&
                     fineAMRStructure[(i * 2) + 1][(j * 2) + 1])
            {
                for (int k = 0; k < 2; k++)
                {
                    for (int l = 0; l < 2; l++)
                    {
                        densityFile << ((rowCellSpacing * 0.5) * ((i * 2) + k)) << " " << ((columnCellSpacing * 0.5) * ((j * 2) + l)) << " " << (rowCellSpacing * 0.5)
                        << " " << (columnCellSpacing * 0.5) << " " << fineSolution[(i * 2) + k][(j * 2) + l].getDensity() << endl;
                        xVelocityFile << ((rowCellSpacing * 0.5) * ((i * 2) + k)) << " " << ((columnCellSpacing * 0.5) * ((j * 2) + l)) << " " << (rowCellSpacing * 0.5)
                        << " " << (columnCellSpacing * 0.5) << " " << fineSolution[(i * 2) + k][(j * 2) + l].getXVelocity() << endl;
                        yVelocityFile << ((rowCellSpacing * 0.5) * ((i * 2) + k)) << " " << ((columnCellSpacing * 0.5) * ((j * 2) + l)) << " " << (rowCellSpacing * 0.5)
                        << " " << (columnCellSpacing * 0.5) << " " << fineSolution[(i * 2) + k][(j * 2) + l].getYVelocity() << endl;
                        zVelocityFile << ((rowCellSpacing * 0.5) * ((i * 2) + k)) << " " << ((columnCellSpacing * 0.5) * ((j * 2) + l)) << " " << (rowCellSpacing * 0.5)
                        << " " << (columnCellSpacing * 0.5) << " " << fineSolution[(i * 2) + k][(j * 2) + l].getZVelocity() << endl;
                        pressureFile << ((rowCellSpacing * 0.5) * ((i * 2) + k)) << " " << ((columnCellSpacing * 0.5) * ((j * 2) + l)) << " " << (rowCellSpacing * 0.5)
                        << " " << (columnCellSpacing * 0.5) << " " << fineSolution[(i * 2) + k][(j * 2) + l].getPressure() << endl;
                    }
                }
            }
        }
    }
    
    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
    zVelocityFile.close();
    pressureFile.close();
}

void EulerAMRTests::outputSolution2DLevel2AMR(tuple<vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >, vector<vector<EulerStateVector> >,
                                              vector<vector<bool> >, vector<vector<bool> >, vector<vector<bool> >> solution)
{
    vector<vector<EulerStateVector> > coarseSolution = get<0>(solution);
    vector<vector<EulerStateVector> > intermediateSolution = get<1>(solution);
    vector<vector<EulerStateVector> > fineSolution = get<2>(solution);
    
    vector<vector<bool> > coarseAMRStructure = get<3>(solution);
    vector<vector<bool> > intermediateAMRStructure = get<4>(solution);
    vector<vector<bool> > fineAMRStructure = get<5>(solution);
    
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
                densityFile << (rowCellSpacing * i) << " " << (columnCellSpacing * j) << " " << rowCellSpacing << " " << columnCellSpacing << " "
                << coarseSolution[i][j].getDensity() << endl;
            }
            else if (intermediateAMRStructure[(i * 2)][(j * 2)])
            {
                densityFile << ((rowCellSpacing * 0.5) * (i * 2)) << " " << ((columnCellSpacing * 0.5) * (j * 2)) << " " << (rowCellSpacing * 0.5)
                << (columnCellSpacing * 0.5) << " " << intermediateSolution[(i * 2)][(j * 2)].getDensity() << endl;
                
                if (intermediateAMRStructure[(i * 2)][(j * 2) + 1])
                {
                    densityFile << ((rowCellSpacing * 0.5) * (i * 2)) << " " << ((columnCellSpacing * 0.5) * ((j * 2) + 1)) << " " << (rowCellSpacing * 0.5) << " "
                    << (columnCellSpacing * 0.5) << " " << intermediateSolution[(i * 2)][(j * 2) + 1].getDensity() << endl;
                }
                else
                {
                    for (int k = 0; k < 2; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            densityFile << ((rowCellSpacing * 0.25) * ((i * 4) + k)) << " " << ((columnCellSpacing * 0.25) * ((j * 4) + 2 + l)) << " "
                            << (rowCellSpacing * 0.25) << " " << (columnCellSpacing * 0.25) << " " << fineSolution[(i * 4) + k][(j * 4) + 2 + l].getDensity() << endl;
                        }
                    }
                }
                
                if (intermediateAMRStructure[(i * 2) + 1][(j * 2)])
                {
                    densityFile << ((rowCellSpacing * 0.5) * ((i * 2) + 1)) << " " << ((columnCellSpacing * 0.5) * (j * 2)) << " " << (rowCellSpacing * 0.5) << " "
                    << (columnCellSpacing * 0.5) << " " << intermediateSolution[(i * 2) + 1][(j * 2)].getDensity() << endl;
                }
                else
                {
                    for (int k = 0; k < 2; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            densityFile << ((rowCellSpacing * 0.25) * ((i * 4) + 2 + k)) << " " << ((columnCellSpacing * 0.25) * ((j * 4) + l)) << " "
                            << (rowCellSpacing * 0.25) << " " << (columnCellSpacing * 0.25) << " " << fineSolution[(i * 4) + 2 + k][(j * 4) + l].getDensity() << endl;
                        }
                    }
                }
                
                if (intermediateAMRStructure[(i * 2) + 1][(j * 2) + 1])
                {
                    densityFile << ((rowCellSpacing * 0.5) * ((i * 2) + 1)) << " " << ((columnCellSpacing * 0.5) * ((j * 2) + 1)) << " " << (rowCellSpacing * 0.5) << " "
                    << (columnCellSpacing * 0.5) << " " << intermediateSolution[(i * 2) + 1][(j * 2) + 1].getDensity() << endl;
                }
                else
                {
                    for (int k = 0; k < 2; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            densityFile << ((rowCellSpacing * 0.25) * ((i * 4) + 2 + k)) << " " << ((columnCellSpacing * 0.25) * ((j * 4) + 2 + l)) << " "
                            << (rowCellSpacing * 0.25) << " " << (columnCellSpacing * 0.25) << " " << fineSolution[(i * 4) + 2 + k][(j * 4) + 2 + l].getDensity() << endl;
                        }
                    }
                }
            }
            else if (fineAMRStructure[(i * 4)][(j * 4)] && fineAMRStructure[(i * 4)][(j * 4) + 1] && fineAMRStructure[(i * 4) + 1][(j * 4)] &&
                     fineAMRStructure[(i * 4) + 1][(j * 4) + 1])
            {
                for (int k = 0; k < 2; k++)
                {
                    for (int l = 0; l < 2; l++)
                    {
                        densityFile << ((rowCellSpacing * 0.25) * ((i * 4) + k)) << " " << ((columnCellSpacing * 0.25) * ((j * 4) + l)) << " "
                        << (rowCellSpacing * 0.25) << " " << (columnCellSpacing * 0.25) << " " << fineSolution[(i * 4) + k][(j * 4) + l].getDensity() << endl;
                    }
                }
                
                if (fineAMRStructure[(i * 4) + 2][(j * 4) + 2] && fineAMRStructure[(i * 4) + 2][(j * 4) + 3] && fineAMRStructure[(i * 4) + 3][(j * 4) + 2] &&
                    fineAMRStructure[(i * 4) + 3][(j * 4) + 3])
                {
                    for (int k = 0; k < 2; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            densityFile << ((rowCellSpacing * 0.25) * ((i * 4) + 2 + k)) << " " << ((columnCellSpacing * 0.25) * ((j * 4) + 2 + l)) << " "
                            << (rowCellSpacing * 0.25) << " " << (columnCellSpacing * 0.25) << " " << fineSolution[(i * 4) + 2 + k][(j * 4) + 2 + l].getDensity() << endl;
                        }
                    }
                }
                else
                {
                    densityFile << ((rowCellSpacing * 0.5) * ((i * 2) + 1)) << " " << ((columnCellSpacing * 0.5) * ((j * 2) + 1)) << " " << (rowCellSpacing * 0.25) << " "
                    << (columnCellSpacing * 0.25) << " " << intermediateSolution[(i * 2) + 1][(j * 2) + 1].getDensity() << endl;
                }
            }
        }
    }
    
    densityFile.close();
}
