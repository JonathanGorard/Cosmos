#ifndef EulerMaterialParameters_hpp
#define EulerMaterialParameters_hpp

class EulerMaterialParameters
{
public:
    EulerMaterialParameters();
    EulerMaterialParameters(double newAdiabaticIndex);
    
    void setAdiabaticIndex(double newAdiabaticIndex);
    
    double getAdiabaticIndex();
    
private:
    double adiabaticIndex;
    double stiffeningParameter;
};

#endif
