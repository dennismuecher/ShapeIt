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
#include <iomanip>
#include <string>

#include "ShapeSetting.h"

#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

class ShapeRho {
    
private:
    ShapeSetting *m_sett;
   
public:
    
    ShapeRho(ShapeSetting* t_setting);
    TGraphErrors *rhoGraph;
    double rhoScaleTrafo = 1;        //scaling factor for transformed graph
    bool wasConvertedFromMeV = false;  // Flag indicating MeV->keV conversion happened
    double originalMaxEnergy = 0.0;    // Original max energy if converted
    void Read();
    void Draw();
    TGraphAsymmErrors* rhoTrafoGraph(double alpha, double alpha_l, double alpha_h);
    TGraphErrors* Transform(double A, double alpha);
    double Eval(double ene);
};
#endif
