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
        
        //manually multiplying experimental level density!!!!!!!
      for (int i=0; i < rhoGraph->GetN(); i++) {
        rhoGraph->SetPoint( i, rhoGraph->GetX()[i], m_sett->rhoScale*rhoGraph->GetY()[i]);
    }
            
        if (m_sett->verbose)
            std::cout <<"Read " << rhoGraph->GetN() <<" data points from level density file " <<m_sett->rhoFileName.c_str() <<std::endl;
    }
    else {
        std::cout <<"No level density file given!"<<std::endl;
        rhoGraph = 0;
    }
}

void ShapeRho::Draw() {
    if (!rhoGraph)
        return;
    std::cout <<"rho graphs contains: " <<rhoGraph->GetN()<<std::endl;
    rhoGraph->SetMarkerStyle(4);
    rhoGraph->SetMarkerColor(kBlue);
    rhoGraph->SetTitle("exp. level density - this work; energy (MeV); level density (1/MeV)");
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
    graph_t->SetTitle("present work; energy (MeV); level density (1/MeV)");
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
    //normalize at 0.5MeV
    double scale = 1 / TMath::Exp(alpha * 0.5);
    for (int i=0; i < rhoGraph->GetN(); i++) {
        double Y = scale * TMath::Exp(alpha * rhoGraph->GetX()[i]) * rhoGraph->GetY()[i];
        double EY = scale * TMath::Exp(alpha * rhoGraph->GetX()[i]) * rhoGraph->GetEY()[i];

        graph_t->SetPoint( graph_t->GetN(), rhoGraph->GetX()[i],Y);
        graph_t->SetPointError( graph_t->GetN()-1, 0, EY);

    }

    return graph_t;
}

double ShapeRho::Eval(double ene) {
    return rhoGraph->Eval(ene);
}
