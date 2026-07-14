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

#include "../Include/ShapeRho.h"

ShapeRho::ShapeRho(ShapeSetting* t_setting) {
    m_sett = t_setting;
    Read();
}

void ShapeRho::Read() {
    if (m_sett->rhoFileName !="") {
        rhoGraph = new TGraphErrors(m_sett->rhoFileName.c_str(),"%lg %lg %lg");
        
        // Auto-detect units: check if input is in MeV or keV
        // Find maximum energy in the dataset
        double maxEnergy = 0;
        for (int i=0; i < rhoGraph->GetN(); i++) {
            double energy = rhoGraph->GetX()[i];
            if (energy > maxEnergy) maxEnergy = energy;
        }
        
        // Heuristic: if max energy < 50, assume MeV (since excitations > 50 MeV are unphysical)
        // Otherwise assume keV (default expected format)
        bool isInMeV = (maxEnergy < 50.0);
        
        if (isInMeV) {
            // Store conversion info for WebSocket notification
            wasConvertedFromMeV = true;
            originalMaxEnergy = maxEnergy;
            
            // Print warning to terminal (visible when running ROOT directly)
            std::cout << "\n" << std::string(70, '=') << std::endl;
            std::cout << "WARNING: NLD file appears to be in MeV units" << std::endl;
            std::cout << "         (max energy = " << std::fixed << std::setprecision(2) << maxEnergy << " MeV)" << std::endl;
            std::cout << "         Auto-converting to keV for internal consistency..." << std::endl;
            std::cout << std::string(70, '=') << "\n" << std::endl;
        } else {
            wasConvertedFromMeV = false;
            originalMaxEnergy = 0.0;
        }
        
        // Apply scaling factor and convert units if needed
        for (int i=0; i < rhoGraph->GetN(); i++) {
            double energy = rhoGraph->GetX()[i];
            double rho = rhoGraph->GetY()[i];
            double rho_err = rhoGraph->GetEY()[i];
            
            if (isInMeV) {
                // Convert from MeV to keV: E[MeV] * 1000 = E[keV], rho[MeV^-1] / 1000 = rho[keV^-1]
                energy = energy * 1000.0;
                rho = (rho / 1000.0) * m_sett->rhoScale;
                rho_err = (rho_err / 1000.0) * m_sett->rhoScale;
            } else {
                // Already in keV, just apply scale factor
                rho = rho * m_sett->rhoScale;
                rho_err = rho_err * m_sett->rhoScale;
            }
            
            rhoGraph->SetPoint(i, energy, rho);
            rhoGraph->SetPointError(i, 0, rho_err);
        }
            
        if (m_sett->verbose)
            std::cout <<"Read " << rhoGraph->GetN() <<" data points from level density file " <<m_sett->rhoFileName.c_str() <<std::endl;
    }
    else {
        std::cout <<"No level density file given!"<<std::endl;
        rhoGraph = 0;
        wasConvertedFromMeV = false;
        originalMaxEnergy = 0.0;
    }
}

void ShapeRho::Draw() {
    if (!rhoGraph)
        return;
    std::cout <<"rho graphs contains: " <<rhoGraph->GetN()<<std::endl;
    rhoGraph->SetMarkerStyle(4);
    rhoGraph->SetMarkerColor(kBlue);
    rhoGraph->SetTitle("exp. level density - this work; energy (keV); level density (1/keV)");
    rhoGraph->SetFillColorAlpha(4,0.5);
    rhoGraph->SetFillStyle(3010);

    //rhoGraph->Draw("AP3*");
    rhoGraph->Draw("same");
}

TGraphAsymmErrors* ShapeRho::rhoTrafoGraph(double alpha, double alpha_l, double alpha_h) {

    if (!rhoGraph)
        return NULL;
    
TGraphErrors* graph_t_mid = Transform(1,alpha);
TGraphErrors* graph_t_low =Transform(1,alpha_l);
TGraphErrors* graph_t_high = Transform(1,alpha_h);

TGraphAsymmErrors* graph_t = new TGraphAsymmErrors();

    for (int i=0; i < rhoGraph->GetN(); i++) {
        graph_t->SetPoint( graph_t->GetN(), rhoGraph->GetX()[i], rhoScaleTrafo*graph_t_mid->GetY()[i]);

        double EY_l =TMath::Abs(graph_t_low->GetY()[i] -graph_t_mid->GetY()[i]);
        double EY_h =TMath::Abs(graph_t_high->GetY()[i] -graph_t_mid->GetY()[i]);
     
        EY_l = TMath::Sqrt(TMath::Power(EY_l,2) + TMath::Power(graph_t_mid->GetEY()[i],2));
        EY_h = TMath::Sqrt(TMath::Power(EY_h,2) + TMath::Power(graph_t_mid->GetEY()[i],2));
        
        graph_t->SetPointError(graph_t->GetN()-1,0,0,rhoScaleTrafo*EY_l,rhoScaleTrafo*EY_h);
    }
    graph_t->SetMarkerStyle(22);
    graph_t->SetMarkerColor(kBlue);
    graph_t->SetTitle("present work; energy (keV); level density (1/keV)");
    graph_t->SetFillColorAlpha(kRed,0.2);
    graph_t->SetFillStyle(3010);
    //printing results to terminal
    if (m_sett->verbose) {
        std::cout <<"Results for transformed level density: \n" <<std::endl;
        for (int i=0; i < graph_t->GetN(); i++) {
            std::cout << graph_t->GetX()[i] <<" " <<graph_t->GetY()[i] <<" " <<graph_t->GetEYhigh()[i] << " " << graph_t->GetEYlow()[i] <<std::endl;
        }
    }
    return graph_t;
}

TGraphErrors* ShapeRho::Transform(double A, double alpha) {
    TGraphErrors* graph_t = new TGraphErrors();
    // alpha is in units of [1/MeV], but energy is now in keV
    // Convert alpha to [1/keV]: alpha[1/MeV] / 1000 = alpha[1/keV]
    double alpha_keV = alpha / 1000.0;
    
    //normalize at 500 keV (was 0.5 MeV)
    double scale = 1 / TMath::Exp(alpha_keV * 500.0);
    for (int i=0; i < rhoGraph->GetN(); i++) {
        double Y = scale * TMath::Exp(alpha_keV * rhoGraph->GetX()[i]) * rhoGraph->GetY()[i];
        double EY = scale * TMath::Exp(alpha_keV * rhoGraph->GetX()[i]) * rhoGraph->GetEY()[i];

        graph_t->SetPoint( graph_t->GetN(), rhoGraph->GetX()[i],Y);
        graph_t->SetPointError( graph_t->GetN()-1, 0, EY);

    }

    return graph_t;
}

double ShapeRho::Eval(double ene) {
    return rhoGraph->Eval(ene);
}
