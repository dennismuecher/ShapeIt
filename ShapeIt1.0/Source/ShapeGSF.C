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
    if (m_sett->doInterpol)
        Sewing();
}

//constructor using the literature values for gSF
ShapeGSF::ShapeGSF(ShapeSetting* t_sett):m_sett(t_sett)
{
    
    levGraph_1 = new TGraphErrors();
    levGraph_2 = new TGraphErrors();
    levGraph   = new TGraphErrors();
    
    TRandom3 *r = new TRandom3(0);

    //read file data into gSF
    if (m_sett->osloFileName == "") {
        std::cout << "No Literature File loaded!"<<std::endl;
        return;
    }
    
    if (m_sett->verbose)
        std::cout <<"\nReading OSLO DATA... " <<std::endl;
    std::ifstream inp;
    inp.open(m_sett->osloFileName.c_str());
    if (inp.is_open() ) {
        
        double e_gamma;
        double oslo_gSF_high, oslo_gSF_low, oslo_gSF, oslo_dgSF;
        
        while ( !inp.eof() ) {
            inp >> e_gamma >> oslo_gSF_high >>oslo_gSF_low;
            
            oslo_gSF = ( oslo_gSF_high + oslo_gSF_low ) / 2;
            oslo_dgSF =( oslo_gSF_high - oslo_gSF_low ) / 2;
            
            if (m_sett->doMC) {
                 levGraph->SetPoint(levGraph->GetN(), e_gamma,r->Gaus(oslo_gSF,oslo_dgSF));
                 levGraph->SetPointError(levGraph->GetN() - 1, 0, 0);
            }
            else {
                levGraph->SetPoint(levGraph->GetN(), e_gamma, oslo_gSF);
                levGraph->SetPointError(levGraph->GetN() - 1, 0, oslo_dgSF);
            }
        }
        
        //apply transformation according to the settings file
        //Transform();
        
        levGraph->SetFillColor(4);
        levGraph->SetFillStyle(3010);
        if (m_sett->verbose)
            levGraph->Print();
    }
}

TGraphErrors* ShapeGSF::GetLevGraph() {
    
    //merge all data in levGraph_1 and levGraph_2 into levGraph
    Merge();
    return levGraph;
}

TGraphErrors* ShapeGSF::GetLevGraph_1() {
    
    return levGraph_1;
}

TGraphErrors* ShapeGSF::GetLevGraph_2() {
    
    return levGraph_2;
}

//resets all stored gSF data
void ShapeGSF::Reset() {
    levGraph->Set(0);
    levGraph_1->Set(0);
    levGraph_2->Set(0);
}

//merges levGraph1 and levGraph2 into one graph and sorts by energy
void ShapeGSF::Merge() {
    
    // in case of literature data, don't do anything here
    if (levGraph_1->GetN() == 0)
        return;
    
    levGraph->Set(0);
    
    //merge
    TObjArray *mArray = new TObjArray();
    mArray->Add(levGraph_1);
    mArray->Add(levGraph_2);
    levGraph->Merge(mArray);
    //sort
    levGraph->Sort();
}

//transforms literature values using setting file
void ShapeGSF::Transform() {

}

//transforms literature gSF values via B*exp(alpha E_gamma)
void ShapeGSF::Transform(double B_t, double alpha_t) {
    
    //gSF was previously transformed via B and Alpha, so only transform according to the change in B_t and Alpha_t
    
    for (int i = 0; i < levGraph->GetN() ; i++ ) {
        levGraph->GetY()[i] *=  B_t / m_B * TMath::Exp(( alpha_t - m_alpha) * levGraph->GetX()[i] / 1000.);
        levGraph->GetEY()[i] *=  B_t / m_B * TMath::Exp(( alpha_t - m_alpha) * levGraph->GetX()[i] / 1000.);
    }
    
    m_B = B_t;
    m_alpha = alpha_t;
}

//fills of levGraph_1 and levGraph_2 with values
void ShapeGSF::FillgSF() {
   
    //clean up matrix
    m_matrix->Reset();
    
    m_matrix->Diag();
    m_matrix->Integrate();
    m_matrix->IntegrateBg();
    m_matrix->IntegrateSquare();
    m_matrix->IntegrateCube();
   
    //if autofit or background subtraction is set, perform autofit using Gauss
    if (m_sett->mode == 2 || m_sett->doBackground)
        m_matrix->FitIntegral();

    double egamma1, egamma2;
    double gSF1, dgSF1, gSF2, dgSF2;

    TRandom3 *r = new TRandom3(0);

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
        
        //in case of MC mode, shuffle gSF according to the gauss distribution
        if (m_sett->doMC) {
            //gaus smear according to fit eerror bars
            gSF1 = r->Gaus(gSF1,dgSF1);
            gSF2 = r->Gaus(gSF2,dgSF2);
            //now set error bars to zero, they are not used for anything anymore
            dgSF1 = 0;
            dgSF2 = 0;
        }
        
        //add a constant fractional error do dgSF as otherwise points with high statistics dominate all fitting results for the slope; in a way, this tries to account for hte fluctuations in the shape method
        
        //if (dgSF1 < 0.25 * gSF1) dgSF1 = 0.25 * gSF1;
        //if (dgSF2 < 0.25 * gSF2) dgSF2 = 0.25 * gSF2;

        
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
    
    //Scale data according to settings file; this also updates levGraph
    
    Scale(m_sett->gSF_norm);
}

double ShapeGSF::Slope(int i) {
    return ( levGraph_1->GetPointY(i) - levGraph_2->GetPointY(i) ) /
    ( levGraph_1->GetPointX(i) - levGraph_2->GetPointX(i) );
}

//sewing of gSF data using the "average" interpolation approach
void ShapeGSF::Sewing() {
    if (m_sett->verbose)
        std::cout <<"\nINTERPOLATION OF gSF DATA... " <<std::endl;
    
    int k = levGraph_1->GetN();
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
    //update levGraph
    Merge();
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
    for (int i = 0; i < levGraph->GetN(); i++) {
        std::cout << levGraph->GetX()[i] <<" "<<levGraph->GetY()[i] <<" "<<levGraph->GetEY()[i]<<std::endl;
    }
    //levGraph->Print();
}

void ShapeGSF::Scale(double factor){
    for (int i = 0; i < levGraph_1->GetN(); i++) {
        levGraph_1->GetY()[i] *= factor;
        levGraph_2->GetY()[i] *= factor;
        levGraph_1->GetEY()[i] *= factor;
        levGraph_2->GetEY()[i] *= factor;
    }
    //update levGraph
    Merge();
}

