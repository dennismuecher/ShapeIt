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

#include "../Include/ShapeAlpha.h"

//constructor using setting file and a matrix
ShapeAlpha::ShapeAlpha(ShapeSetting* t_sett, ShapeCollector* t_collector):m_sett(t_sett), m_collector(t_collector)
{
    Chi2Loop();
}

void ShapeAlpha::FindMinimum() {
    
    double x_min = 0, y_min = 0;
    for (int i = 0; i < chi2Graph->GetN(); i++) {
        double x = chi2Graph->GetX()[i];
        double y = chi2Graph->GetY()[i];
        if ( (y < y_min) || (y_min == 0) ) {
            y_min = y;
            x_min = x;
        }
    }
    minChi2 = y_min;
    minAlpha = x_min;
}

//calculates the chi2 values relative to literature values for a range of alpha values
void ShapeAlpha::Chi2Loop() {
    
    double alphaX[nOfPoints*nOfExi];
    double chi2Y[nOfPoints*nOfExi];
    double enes[nOfPoints*nOfExi];
    int pointC = 0;
    
    //store current values of alpha and exiEne[0]
    double alpha_bak = m_sett->lit_alpha;
    double exi_bak = m_sett->exiEne[0];
    
    for (int i = 0; i < nOfPoints; i++) {
        
        alphaX[pointC] = 2*i*alphaRange / nOfPoints - alphaRange;
        
        //update alpha in settings file
        m_sett->lit_alpha = alphaX[pointC];
        
        //apply transformation
        m_collector->Transform(m_sett->lit_norm, m_sett->lit_alpha);
        chi2Y[pointC] = m_collector->getChi2();
        enes[pointC] = m_sett->exiEne[0];
        pointC++;
        //std::cout <<" i = " << i << " j = " << j << "pointC = " <<pointC <<"alphaX[i]" << alphaX[pointC] << "exi[0] " << sett->exiEne[0] << "lit_chi2 " << chi2Y[pointC-1]<<  std::endl;
    }
    
    //store values in chi2Graph
    delete chi2Graph;
    chi2Graph = new TGraph(nOfPoints, alphaX, chi2Y);
    chi2Graph->SetMarkerStyle(4);
    chi2Graph->SetMarkerColor(kRed);
    
    //restore m_sett values
    m_sett->exiEne[0] = exi_bak;
    m_sett->lit_alpha = alpha_bak;
    
    //find minimum values for alpha and chi2
    FindMinimum();
    
}

