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

#include <iostream>

#include "ShapeSetting.h"
#include "ShapeMatrix.h"

#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

class ShapeRho {
    
private:
    ShapeSetting *m_sett;
    ShapeMatrix *m_matrix;
   
public:
    
    ShapeRho(ShapeSetting* t_setting, ShapeMatrix* t_matrix);
    TGraphErrors *rhoGraph;
    void Read();
    void Draw();
    TGraphAsymmErrors* rhoTrafoGraph(double alpha, double alpha_l, double alpha_h);
    TGraphErrors* Transform(double A, double alpha);
    double Eval(double ene);
};
#endif
