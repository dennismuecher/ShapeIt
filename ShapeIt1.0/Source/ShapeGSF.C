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

#include "../Include/ShapeGSF.h"

//constructor using setting file and a matrix
ShapeGSF::ShapeGSF(ShapeSetting* t_sett, ShapeMatrix* t_matrix):m_sett(t_sett), m_matrix(t_matrix)
{
    
    levGraph_1 = new TGraphErrors();
    levGraph_2 = new TGraphErrors();
    levGraph   = new TGraphErrors();

    levGraph_1->SetMarkerStyle(22);
    levGraph_1->SetMarkerSize(2);
    levGraph_1->SetMarkerColor(6);
    levGraph_2->SetMarkerStyle(22);
    levGraph_2->SetMarkerSize(2);
 
    if (m_sett->colour)
        levGraph_2->SetMarkerColor(7);
    else
        levGraph_2->SetMarkerColor(6);
    
    FillgSF();
    if (m_sett->DoInterpol)
        Sewing();
}

//merges levGraph1 and levGraph2 into one graph and sorts by energy
void ShapeGSF::Merge() {
    
    TObjArray *mArray = new TObjArray();
    mArray->Add(levGraph_1);
    mArray->Add(levGraph_2);
    levGraph->Merge(mArray);
    levGraph->Sort();
}

//fills of levGraph_1 and levGraph_2 with values
void ShapeGSF::FillgSF() {
   
    m_matrix->Integrate();
    m_matrix->IntegrateBg();
    m_matrix->IntegrateSquare();
    m_matrix->IntegrateCube();
   
    //if autofit or background subtraction is set, perform autofit using Gauss
    if (m_sett->mode == 2 || m_sett->doBackground)
        m_matrix->FitIntegral();

    double egamma1, egamma2;
    double gSF1, dgSF1, gSF2, dgSF2;
    
    for (int i = 0; i < m_matrix->integral1Cube.size(); i++ ) {

        //if any of the two peak areas is below minCounts, skip this entire gSF pair
        if ( getBgRatio(i, 1) * m_matrix->integral1[i] < m_sett->minCounts ||
             getBgRatio(i, 2) * m_matrix->integral2[i] < m_sett->minCounts ||
             getBgRatio(i, 1) * m_matrix->integral1[i] <= 0 ||
             getBgRatio(i, 2) * m_matrix->integral2[i] <= 0 )
            continue;

        //calculating the average E_g, weighted by gSF, i.e. egamma_avg = SUM (E_g*gSF) / SUM(gSF). E_g*gSF is contained in the gSFSquare matrix
        egamma1 = m_matrix->integral1Square[i] / m_matrix->integral1Cube[i];
        egamma2 = m_matrix->integral2Square[i] / m_matrix->integral2Cube[i];

        if (egamma1 < m_sett->gammaEne[0] || egamma2 < m_sett->gammaEne[0] )
            continue;
        if (egamma1 > m_sett->gammaEne[1] || egamma2 > m_sett->gammaEne[1] )
            continue;
        
        //calculate gamma ray strength
        gSF1 = getBgRatio(i, 1) * m_matrix->integral1Cube[i] * m_sett->getEffCor(egamma1, 1);
        dgSF1 = TMath::Power(1./m_matrix->integral1[i], 0.5) * gSF1;

        gSF2 = getBgRatio(i, 2) * m_matrix->integral2Cube[i] * m_sett->getEffCor(egamma2, 2);
        dgSF2 = TMath::Power(1./m_matrix->integral2[i], 0.5) * gSF2;

        //recalculate gSF in case autofit is activated
        
        if (m_sett->mode == 2 && m_sett->doBackground) {
            gSF1 = m_sett->getEffCor(egamma1, 1) * m_matrix->fit_integral1Net[i] * m_matrix->integral1Cube[i]  / m_matrix->integral1[i];
            
            gSF2 = m_sett->getEffCor(egamma2, 2) * m_matrix->fit_integral2Net[i] * m_matrix->integral2Cube[i] / m_matrix->integral2[i];
                
            if (m_matrix->fit_integral1Net[i] < 0 || m_matrix->fit_integral2Net[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                continue;
            }
        }
        //autofit without background
        else if( m_sett->mode == 2 && !m_sett->doBackground) {
            gSF1 = m_sett->getEffCor(egamma1, 1) * m_matrix->fit_integral1[i] * m_matrix->integral1Cube[i]  / m_matrix->integral1[i];
            
            gSF2 = m_sett->getEffCor(egamma2, 2) * m_matrix->fit_integral2[i] * m_matrix->integral2Cube[i]  / m_matrix->integral2[i];
            
            if (m_matrix->fit_integral1[i] < 0 || m_matrix->fit_integral2[i] < 0) {
                std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                continue;
            }
        }
        
        //filling results into levGraph vector;
        levGraph_1->SetPoint(levGraph_1->GetN(), egamma1, gSF1);
        levGraph_1->SetPointError(levGraph_1->GetN()-1, 0, dgSF1);
        levGraph_2->SetPoint(levGraph_2->GetN(), egamma2, gSF2);
        levGraph_2->SetPointError(levGraph_2->GetN()-1, 0, dgSF2);
    
        if (m_sett->verbose) {
            std::cout <<"Bin: " <<i+1<<std::endl;
            std::cout <<"energies: " <<egamma1 << " " << egamma2 <<std::endl;
            std::cout <<"gSF1: " <<gSF1 << "+- " << dgSF1 <<std::endl;
            std::cout <<"gSF2: " <<gSF2 << "+- " << dgSF2 <<std::endl;
        }
    }
}

//what exactly does this do?
double ShapeGSF::Slope(int i) {
    return ( levGraph_1->GetPointY(i) - levGraph_2->GetPointY(i) ) /
    ( levGraph_1->GetPointX(i) - levGraph_2->GetPointX(i) );
}

//this is doing the sewing of gSF data
void ShapeGSF::Sewing() {
    if (m_sett->verbose)
        std::cout <<"\nINTERPOLATION OF gSF DATA... " <<endl;
    
    int k = levGraph_1[j]->GetN();
    //if there is only one pair of non-zero gSF values, don't do anything
    if (k < 2) {
        std::cout <<"Nothing to interpolate in this iteration!" <<std::endl;
        return;
    }
    
    //interpolate from first bin to higher energies; there is nothing to do for the first bin
    for (int i = 1; i < k; i++) {
        // do average interpolation; this scaling makes the interpolation independent of the direction
        double scale_avg = ( 2 * levGraph_1->GetPointY(i-1) - ( ( levGraph_1->GetPointX(i-1) - levGraph_2->GetPointX(i)) * Slope(i-1) ) ) / ( 2 * levGraph_2->GetPointY(i) + ( (levGraph_1->GetPointX(i-1) - levGraph_2->GetPointX(i)) * Slope(i) ) );
       
        levGraph_1->GetY()[i] *= scale_avg;
        levGraph_1->GetEY()[i] *= scale_avg;
        levGraph_2->GetY()[i] *= scale_avg;
        levGraph_2->GetEY()[i] *= scale_avg;
    }
}

//calculates the ratio of Peak to background for bin i

double ShapeGSF::getBgRatio(int bin, int level) {
    
    if (!m_sett->doBackground)
        return 1;
    
    if (level == 1 ) {
        
        double peak = m_matrix->integral1[bin];
        double backgr = m_matrix->fit_integral1Bg[bin];
        if (peak > 0)
            return (peak - backgr) / peak;
        else
            return 0;
    }
    
    else if (level == 2 ) {
        
        double peak = m_matrix->integral2[bin];
        double backgr =m_matrix->fit_integral2Bg[bin];
        if (peak > 0)
            return (peak - backgr) / peak;
        else
            return 0;
    }
    return 0;
}

//prints the gSF values, sorted for egamma
void ShapeGSF::Print() {
    std::cout << "\n\nResults for gamma ray strength function: " <<std::endl;
    levGraph_1->Print();
    levGraph_2->Print();
}

void ShapeGSF::Scale(double factor){

    for (int i = 0; i < levGraph_1->GetN(); i++) {
        levGraph_1->GetY()[i] *= factor;
        levGraph_2->GetY()[i] *= factor;
        levGraph_1->GetEY()[i] *= factor;
        levGraph_2->GetEY()[i] *= factor;
    }
}

