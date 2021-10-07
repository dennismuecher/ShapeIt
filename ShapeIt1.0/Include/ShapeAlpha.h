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

#ifndef SHAPEALPHA_H
#define SHAPEALPHA_H
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeCollector.h"

#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>

#include <TObjArray.h>
#include <TH1.h>

/* **************************************************************************
*  class ShapeALPHA                                                         *
*  class offering tools to determine the best slope "alpha" in comparison   *
*        to literature Oslo data                                            *
*                                                                           *
*                                                                           *
*************************************************************************/

class ShapeAlpha {
    
private:
    
    ShapeSetting*       m_sett;
    ShapeMatrix*        m_matrix;
    ShapeCollector*     m_collector;
    TGraph*             chi2Graph;      //1-d graph of chi2 values for various alpha values
    
    double const        alphaRange = 0.5;   //alpha will be varied in Chi2Loop() according to this
    int const           nOfPoints = 200;    //number of steps in the chi2 evaluation
    int const           nOfExi = 1;         //number of excitation energies;  not impl.
    double const        exiStep = 1;        //exi step size (in keV) in Chi2Loop(); not impl.
    double              minChi2 = 0;        //minimum chi2 value
    double              minAlpha = 1;       //alpha value for which chi2 is at minimum
    
    void                Chi2Loop();
    void                FindMinimum();      //finds minChi2 and minAlpha
    
public:
    
    //constructor
    ShapeAlpha(ShapeSetting* t_sett, ShapeCollector* t_collector);
    
    TGraph*             getChi2Graph() {return chi2Graph;}
    double              getMinChi2()   {return minChi2;}
    double              getMinAlpha()   {return minAlpha;}

    
};
#endif


