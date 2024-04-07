#include "EulerMaterialParameters.hpp"

EulerMaterialParameters::EulerMaterialParameters()
{
    adiabaticIndex = 1.0;
}

EulerMaterialParameters::EulerMaterialParameters(double newAdiabaticIndex)
{
    adiabaticIndex = newAdiabaticIndex;
}

void EulerMaterialParameters::setAdiabaticIndex(double newAdiabaticIndex)
{
    adiabaticIndex = newAdiabaticIndex;
}

double EulerMaterialParameters::getAdiabaticIndex()
{
    return adiabaticIndex;
}
