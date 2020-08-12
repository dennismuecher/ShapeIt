#ifndef SHAPERHO_H
#define SHAPERHO_H
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"

#include <TGraphErrors.h>

class ShapeRho {
    
private:
    ShapeSetting *m_sett;
    ShapeMatrix *m_matrix;
    double e_snail = 9.427;             //neutron separation energy (MeV)
    double rho_snail = 5.89E4;           // level density at e_snail (1/MeV)
    double drho_snail = 1.18E4;          // error level denisty at e_snail (1/NeV)
public:
    
    ShapeRho(ShapeSetting* t_setting, ShapeMatrix* t_matrix);
    TGraphErrors *rhoGraph;
    void Read();
    void Draw();
    void Draw(double alpha, double alpha_l, double alpha_h);
    TGraphErrors* Transform(double A, double alpha);
    double Eval(double ene);
};
#endif
