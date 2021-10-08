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

ShapeRho::ShapeRho(ShapeSetting* t_setting, ShapeMatrix* t_matrix) {
    m_sett = t_setting;
    m_matrix = t_matrix;
    Read();
}

void ShapeRho::Read() {
    if (m_sett->rhoFileName !="") {
        rhoGraph = new TGraphErrors(m_sett->rhoFileName.c_str(),"%lg %lg %lg");
        if (m_sett->addSnail) {
            rhoGraph->SetPoint(rhoGraph->GetN(), e_snail, rho_snail);
            rhoGraph->SetPointError(rhoGraph->GetN()-1, 0, drho_snail);
        }
        if (m_sett->verbose)
            std::cout <<"Read " << rhoGraph->GetN() <<" data points from level density file " <<m_sett->rhoFileName.c_str() <<std::endl;
    }
    else {
        std::cout <<"No level density file given!"<<std::endl;
        rhoGraph = 0;
    }
}

void ShapeRho::Read(const char* filename) {
    
    rhoGraph = new TGraphErrors(filename,"%lg %lg %lg");
    if (m_sett->addSnail) {
        rhoGraph->SetPoint(rhoGraph->GetN(), e_snail, rho_snail);
        rhoGraph->SetPointError(rhoGraph->GetN()-1, 0, drho_snail);
        
        if (m_sett->verbose)
            std::cout <<"Read " << rhoGraph->GetN() <<" data points from level density file " <<m_sett->rhoFileName.c_str() <<std::endl;
    }
}


void ShapeRho::Draw() {
    if (!rhoGraph)
        return;
    rhoGraph->SetMarkerStyle(4);
    rhoGraph->SetMarkerColor(kBlue);
    rhoGraph->SetTitle("level density; energy (MeV); level density (1/MeV)");
    rhoGraph->SetFillColorAlpha(4,0.5);
    rhoGraph->SetFillStyle(3010);

    rhoGraph->Draw("AP3*");
    //rhoGraph->Draw("same");
}

void ShapeRho::Draw(double alpha, double alpha_l, double alpha_h, int m_color) {


    if (!rhoGraph)
        return;
    
TGraphErrors* graph_t_mid = Transform(1,alpha);
TGraphErrors* graph_t_low =Transform(1,alpha_l);
TGraphErrors* graph_t_high = Transform(1,alpha_h);

TGraphAsymmErrors* graph_t = new TGraphAsymmErrors();

    for (int i=0; i < rhoGraph->GetN(); i++) {
        graph_t->SetPoint( graph_t->GetN(), rhoGraph->GetX()[i], graph_t_mid->GetY()[i]);

        double EY_l =TMath::Abs(graph_t_low->GetY()[i] -graph_t_mid->GetY()[i]);
        double EY_h =TMath::Abs(graph_t_high->GetY()[i] -graph_t_mid->GetY()[i]);
        if (graph_t_mid->GetX()[i] < e_snail) {
            EY_l = TMath::Sqrt(TMath::Power(EY_l,2) + TMath::Power(graph_t_mid->GetEY()[i],2));
            EY_h = TMath::Sqrt(TMath::Power(EY_h,2) + TMath::Power(graph_t_mid->GetEY()[i],2));
        }
        graph_t->SetPointError(graph_t->GetN()-1,0,0,EY_l,EY_h);
    }
    graph_t->SetMarkerStyle(4);
    if (m_color == 1)
        graph_t->SetMarkerColor(kRed);
    else if (m_color == 2)
           graph_t->SetMarkerColor(kBlue);
    else if (m_color == 3)
        graph_t->SetMarkerColor(kPink);
    
    graph_t->SetMarkerColor(m_color);
    graph_t->SetTitle("level density; energy (MeV); level density (1/MeV)");
    graph_t->SetFillColorAlpha(kRed,0.5);
    graph_t->SetFillStyle(3010);
    //graph_t->Draw("aP3*");
    graph_t->Draw("same");

    //printing results to terminal
    std::cout <<"Results for transformed level density: \n" <<std::endl;
    for (int i=0; i < graph_t->GetN(); i++) {
        std::cout << graph_t->GetX()[i] <<" " <<graph_t->GetY()[i] <<" " <<graph_t->GetEYhigh()[i] << " " << graph_t->GetEYlow()[i] <<std::endl;
    }
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


//temp attempt to get level density from first generation matrix
//previously used in ShapeGSF class

/*void ShapeGSF::Rho(){
    TH2* fgMatrix;
    double MeV=1;
    rhoDiagram = new TGraphErrors();
    TFile* fgFile =new TFile("../FirstGenerationFiles/Ge76_fg.root","read");
    fgFile->GetObject("matrix1", fgMatrix);
 
    TAxis *xaxis = fgMatrix->GetXaxis();
    TAxis *yaxis = fgMatrix->GetYaxis();
    
    Int_t XNum  = fgMatrix->GetNbinsX();
    Int_t YNum  = fgMatrix->GetNbinsY();
    
    double y_width = yaxis->GetBinWidth(YNum) * MeV;
    
    double eMax_x  =  ( fgMatrix->GetXaxis()->GetBinCenter(XNum) * MeV) + (fgMatrix->GetXaxis()->GetBinWidth(XNum) * MeV /2 ) ;
    double eMax_y  =  ( fgMatrix->GetYaxis()->GetBinCenter(YNum) * MeV) + (fgMatrix->GetYaxis()->GetBinWidth(YNum) * MeV /2 ) ;
    double eMin_x = fgMatrix->GetXaxis()->GetXmin();
    double eMin_y = fgMatrix->GetYaxis()->GetXmin();
    
    double ene0 = sett->exiEne[0];
    double ene1 = sett->exiEne[1];
    double esize = sett->exi_size[0];
   
    double xbins = (int) (eMax_x - eMin_x) / fgMatrix->GetXaxis()->GetBinWidth(1);
    int ybins = (int) ( ene1 - ene0 ) / esize;
    std::cout <<"ybins: " <<ybins <<std::endl;

    TH1F* rhoGraph = new TH1F("rho","rho",100,0,ene1);

    
    matrixEx= new TH2F("matrix_ex","Gamma ray energy vs excitation energy",xbins,eMin_x, eMax_x, ybins, ene0, ene1);
    
    //fill the diag matrices and histos by looping over fgMatrix
    int counter = 0;
    for (int j = 1; j < fgMatrix->GetNbinsY() + 1; j++) {
        double Y = yaxis->GetBinCenter(j);
       
        if (Y > ene0 && Y < ene1) {
            //std::cout <<"Y: " <<Y <<std::endl;
            int y_t = (int) (Y - ene0) / esize;
            //std::cout <<"y_t: " <<y_t <<std::endl;
            
            
            double y_low = (ene0 + (y_t*esize))/y_width;
            double y_high = (ene0 + ((y_t+1)*esize))/y_width;
            
            double integ = fgMatrix->Integral(1, XNum,y_low, y_high);

            //std::cout <<"integ: " <<integ <<std::endl;
            //std::cout <<"bin_low " << y_low  <<std::endl;
            //std::cout <<"bin_high " << y_high <<std::endl;


            for (int i = 1; i < fgMatrix->GetNbinsX() + 1; i++) {
                double X = xaxis->GetBinCenter(i);
                double eRho = Y - X;
                if (eRho > 200 && eRho < ene1) {
                double g = InterpolValueSort(X)[0];
                if (g !=0 && integ !=0 && X!=0) {
                    //double rho = fgMatrix->GetBinContent(i,j) /g /integ * eRho /TMath::Power(X,3);
                    double rho = fgMatrix->GetBinContent(i,j) /g /integ  /TMath::Power(X,3);
                    if (rho > 0) {
                        //std::cout <<"Bin content: " <<fgMatrix->GetBinContent(i,j)<<std::endl;
                        counter++;
                        rhoGraph->Fill(eRho,rho);
                        //rhoDiagram->SetPoint(rhoDiagram->GetN(),eRho,rho);
                        //std::cout <<"E_x: " <<Y <<" "<<"E_g: " <<X<<" "<<"gSF: " <<g<< " "<<"rho: " <<rho<<std::endl;
                    }
                }
                }
            }
        }
    }
    
   
    for(int i = 1; i <= rhoGraph->GetNbinsX(); i++ ) rhoDiagram->SetPoint(i-1, rhoGraph->GetBinCenter(i), rhoGraph->GetBinContent(i));
    
    
    //rhoDiagram->Draw("A*");
   // std::cout <<"Counter: " <<counter <<std::endl;
    //rhoGraph->Draw("HIST");
}*/

/*TGraphErrors* ShapeGSF::GetRhoDiagram() {
    return rhoDiagram;
}*/





