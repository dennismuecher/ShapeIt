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

#ifndef SHAPEGSF_H
#define SHAPEGSF_H

#include <iostream>

#include "ShapeSetting.h"
#include "ShapeMatrix.h"

#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TRandom3.h>
#include <TGraphAsymmErrors.h>

#include <TObjArray.h>
#include <TH1.h>
#include <TRandom3.h>


/* **************************************************************************
*  class ShapeGSF                                                           *
*  base class to store gamma ray strength data for level 1 and level 2;     *
*        also offers the "sewing" of the gSF data                           *
*                                                                           *
*                                                                           *
*************************************************************************/

class ShapeGSF {
    
private:
    
    ShapeSetting*       m_sett;
    ShapeMatrix*        m_matrix;
    TGraphErrors*       levGraph_1;          //TGraph containing the gSF data for level1
    TGraphErrors*       levGraph_2;          //TGraph containing the gSF data for level2
    TGraphErrors*       levGraph;            //TGraph to contain the gSF data for both levels
    double              m_B = 1;             //transfomration gSF scaling
    double              m_alpha = 0;         //trnasformation gSF exponential slope
    
    void                FillgSF();
    void                Merge();
    void                Sewing();
    double              getBgRatio(int bin, int level);
    double              Slope(int i);
    
public:
    
    ShapeGSF(ShapeSetting* t_sett, ShapeMatrix* t_matrix);
    ShapeGSF(ShapeSetting* t_sett);

    void                Print();
    void                Reset();
    void                Scale(double factor);
    void                Transform(double B_t, double alpha_t);
    void                Transform();

    void                SetESize(double ene)    {m_matrix->SetESize(ene);}
    void                SetEne0(double ene)     {m_matrix->SetEne0(ene);}
    void                SetEne1(double ene)     {m_matrix->SetEne1(ene);}


    TGraphErrors*       GetLevGraph();
    TGraphErrors*       GetLevGraph_1();
    TGraphErrors*       GetLevGraph_2();
};
#endif

