/* ***********************************************************************
* Copyright (C) 2019-2021, Dennis Muecher.                               *
* All rights reserved.                                                   *
*                                                                        *
* This program is free software: you can redistribute it and/or modify   *
* it under the terms of the GNU General Public License as published by   *
* the Free Software Foundation, either version 3 of the License, or      *
* (at your option) any later version.                                    *
* You should have received a copy of the GNU General Public License      *
* along with this program. If not, see  http://www.gnu.org/licenses/.    *
*************************************************************************/

#ifndef SHAPERHO_H
#define SHAPERHO_H
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"

#include <TGraphErrors.h>

class ShapeRho {
    
private:
    ShapeSetting *m_sett;
    ShapeMatrix *m_matrix;
    //76Ge
    //double e_snail = 9.427;             //neutron separation energy (MeV)
    //double rho_snail = 5.89E4;           // level density at e_snail (1/MeV)
    //double drho_snail = 1.18E4;          // error level denisty at e_snail (1/NeV)
    
    //88Kr
    //double e_snail = 7.05;             //neutron separation energy (MeV)
    //double rho_snail = 3000;           // level density at e_snail (1/MeV)
    //double drho_snail = 300;          // error level denisty at e_snail (1/NeV)
    
    //140Ba
    double e_snail = 6.40;             //normalizationenergy (MeV)
    double rho_snail = 20000;           // level density at e_snail (1/MeV)
    double drho_snail = 2000;          // error level denisty at e_snail (1/NeV)
   
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
