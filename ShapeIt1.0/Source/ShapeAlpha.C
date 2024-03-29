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

}

void ShapeAlpha::FindMinimum() {
    double x_min = 0, y_min = 0;
    int i_min = 0;
    for (int i = 0; i < chi2Graph->GetN(); i++) {
        double x = chi2Graph->GetX()[i];
        double y = chi2Graph->GetY()[i];
     
        if ( (y < y_min) || (y_min == 0) ) {
            y_min = y;
            x_min = x;
            i_min = i;
        }
    }

    if (i_min ==0)
        std::cout <<"Warning: minimum chi2 value found for first chi2Graph point!" <<std::endl;
    else if (i_min == chi2Graph->GetN()-1)
        std::cout <<"Warning: minimum chi2 value found for last chi2Graph point!" <<std::endl;
    
    minChi2 = y_min;
    minAlpha = x_min;
    minAlphaPoint = i_min;
    
}

//provides a TPaveText window with infos on the chi2 analysis results
TPaveText* ShapeAlpha::getPaveTextChi2() {
    
    TPaveText *t=new TPaveText(0.45,0.65,0.9,0.9,"brNDC");
    t->SetTextSize(0.025);
    t->SetTextAlign(13);
    t->SetFillColor(10);
    t->SetTextColor(61);
    
    t->AddText(Form("Minimum #chi^{2} value: %4.2f ",minChi2 ));
    t->AddText(Form("Slope #alpha at minimum: %4.2f ",minAlpha ));
    t->AddText(Form("Chance for larger #chi^{2}: %4.2f %% ",100*minAlphaProb ));
    t->AddText(Form("Chance for smaller #chi^{2}: %4.2f %% ",100*(1-minAlphaProb) ));
    
    if (minAlpha_low < 1E3 && minAlpha_high < 1E3)
        t->AddText(Form("%2.0f %% confidence level #alpha: [%4.2f, %4.2f] ",100*(1-conf_level ),minAlpha_low, minAlpha_high ));
    else
        t->AddText(Form("%2.0f %% confidence level #alpha: NaN! #chi^{2} too large?" ,100*(1-conf_level )));
    return t;
}

//determines the range of alpha values around the minimum which meet the chi2 confidence level conf_level
void ShapeAlpha::AlphaRange() {
    
    minAlpha_high = 1E3;
    minAlpha_low = 1E3;
    
    int n = chi2Graph->GetN();
    int nOfDegFreedom;
    if (m_sett->displayAvg)
        nOfDegFreedom = m_collector->GetNSmooth() -1;
    else
        nOfDegFreedom = m_collector->GetN() -1;

    if (m_sett->verbose) {
        std::cout <<" minChi2 = " << minChi2 <<std::endl;
        std::cout <<" minAlpha = " << minAlpha <<std::endl;
        std::cout <<" nOfDegFreedom = " << nOfDegFreedom <<std::endl;
        std::cout <<"minAlphaPoint = " << minAlphaPoint <<std::endl;
        std::cout <<"confidence level = " << conf_level <<std::endl;

    }

    minAlphaProb = TMath::Prob(minChi2,nOfDegFreedom);
    
    if (m_sett->verbose) {
        if (minAlphaProb <= conf_level) {
            std::cout <<"Cannot determine alpha confidence level: chi2 probability for minimum point is " <<minAlphaProb << " compared to the confidence level " << conf_level <<std::endl;
            return;
        }
        if (minAlphaPoint == 0) {
            std::cout <<"Cannot determine alpha confidence level: chi2 minimum is at first chi2Graph point! " <<std::endl;
            return;
        }
        
        if (minAlphaPoint == n-1) {
            std::cout <<"Cannot determine alpha confidence level: chi2 minimum is at last chi2Graph point! " <<std::endl;
            return;
        }
    }
   
    double alpha,chi2,chi2_prob;

    
    //find the higher limit of alpha for which the chi2 is within the given confidence level
    for ( int i = minAlphaPoint; i < n; i++) {
        alpha  = chi2Graph->GetX()[i];
        chi2 = chi2Graph->GetY()[i];
         
        // probability that an observed Chi-squared exceeds the value chi2 by chance, even for a correct model
        chi2_prob = TMath::Prob(chi2,nOfDegFreedom);
        if ( chi2_prob > conf_level)
            minAlpha_high = alpha;
    }
    if (minAlpha_high == minAlpha) {
        std::cout <<"Cannot determine alpha confidence level: chi2 minimum too close to confidence level!" <<std::endl;
        minAlpha_high=1E3;
    }
    
    //find the lower limit of alpha for which the chi2 is within the given confidence level
    for ( int i = minAlphaPoint; i > 0; i--) {
        alpha  = chi2Graph->GetX()[i];
        chi2 = chi2Graph->GetY()[i];
      
        // probability that an observed Chi-squared exceeds the value chi2 by chance, even for a correct model
        chi2_prob = TMath::Prob(chi2,nOfDegFreedom);
        if ( chi2_prob > conf_level)
            minAlpha_low = alpha;
    }
    if (minAlpha_low == minAlpha) {
        std::cout <<"Cannot determine alpha confidence level: chi2 minimum too close to confidence level!" <<std::endl;
        minAlpha_low =-1E3;
    }
    
    if (m_sett->verbose)
        std::cout <<"Confidence level for alpha = " << minAlpha_low <<" to "<<minAlpha_high <<std::endl;
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
        
        alphaX[pointC] = i * (m_sett->alphaLimit[1]-m_sett->alphaLimit[0] ) / nOfPoints + m_sett->alphaLimit[0];

        //update alpha in settings file
        m_sett->lit_alpha = alphaX[pointC];
        
        //apply transformation
        m_collector->Transform(m_sett->lit_norm, m_sett->lit_alpha);
        chi2Y[pointC] = m_collector->getChi2();
        enes[pointC] = m_sett->exiEne[0];
        
        //std::cout <<" i = " << i << "pointC = " <<pointC <<" alphaX[i] " << alphaX[pointC] << " exi[0] " << m_sett->exiEne[0] << " chi2: " << chi2Y[pointC]<<  std::endl;
        pointC++;
    }
    
    //store values in chi2Graph
    //delete chi2Graph;
    chi2Graph = new TGraph(nOfPoints, alphaX, chi2Y);
    chi2Graph->SetMarkerStyle(22);
    chi2Graph->SetMarkerColor(6);
    chi2Graph->SetTitle("#chi^{2} fit to literature gSF data; slope #alpha; #chi^{2}" );
    //restore m_sett values
    m_sett->exiEne[0] = exi_bak;
    m_sett->lit_alpha = alpha_bak;
    //find minimum values for alpha and chi2
    FindMinimum();

    //find lower and upper bounds for alpha
    if (!m_sett->doMC)
        AlphaRange();
}

