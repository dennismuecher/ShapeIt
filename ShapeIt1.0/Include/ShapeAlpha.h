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
    TGraph*             chi2Graph;               //1-d graph of chi2 values for various alpha values
    
    double const        alphaRange = 1.;        //alpha will be varied in Chi2Loop() according to this
    int const           nOfPoints = 100;         //number of steps in the chi2 evaluation
    int const           nOfExi = 1;              //number of excitation energies;  not impl.
    double const        exiStep = 1;             //exi step size (in keV) in Chi2Loop(); not impl.
    double              minChi2 = 1;             //minimum chi2 value
    double              minAlpha = 1;             //alpha value for which chi2 is at minimum
    double              minAlpha_low;             //lower bound for alpha within chi2 confidence level
    double              minAlpha_high;            //higher bound for alpha within chi2 confidence level
    
    int                 minAlphaPoint = 0;       //index in chi2Graph with minimum chi2 value
    double              conf_level = 0.1;         //confidence level for alpha range; a confidence level of 0.1 means that the chance that a given chi2 value is larger than the given value is 0.1 = 10%
    
    
    void                Chi2Loop();
    void                FindMinimum();      //finds minChi2 and minAlpha
    void                AlphaRange();       //finds lower and upper bound for alpha following chi2 criterium
    
public:
    
    //constructor
    ShapeAlpha(ShapeSetting* t_sett, ShapeCollector* t_collector);
    
    TGraph*             getChi2Graph() {return chi2Graph;}
    double              getMinChi2()   {return minChi2;}
    double              getMinAlpha()   {return minAlpha;}


};
#endif


