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
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>

#include <TObjArray.h>
#include <TH1.h>

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
    
    TGraphErrors*       levGraph_1;                      //TGraph containing the gSF data for level1
    TGraphErrors*       levGraph_2;                      //TGraph containing the gSF data for level2
    TGraphErrors*       levGraph;                        //TGraph to contain the gSF data for both levels
    
    TRandom3            r;

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

    TGraphErrors*       GetLevGraph();
    TGraphErrors*       GetLevGraph_1();
    TGraphErrors*       GetLevGraph_2();
};
#endif

