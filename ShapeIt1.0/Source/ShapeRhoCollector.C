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

#include "../Include/ShapeRhoCollector.h"
#include "../Include/ShapeMultiGraph.h"
#include "../Include/ShapeTalys.h"


//consstructor

ShapeRhoCollector::ShapeRhoCollector(ShapeSetting* p_sett) {
    sett = p_sett;
    rho = new ShapeRho(sett);
    collGraph = new TMultiGraph();
   
    //initialize ldmodel objects and set display switch
    for (int i = 0; i < 6; i++) {
        ldmodel[i] = NULL;
        ldmodelDisplay[i] = true;
        ldmodelChi2Display[i] = false;
    }
    ldmodelChi2Display[4] = true;

    discreteLevel = NULL;
    discBand = NULL;
    
    //read the discrete level density information from file
    ReadDiscrete();
    
    //supplies the original and transformed experimental level density using slope from settings file
    FillTrafoGraph();
    
    //run Talys calculations, as requested
    FillTalys();
    
    //determine exp normalization and re-read level densities
    //ExpNorm();
    
    //run MC to find experimental errors
    ExpMC();
    
    //determine band of discrete level fits
    DiscreteBand();
    
    //add literature Graph
    //AddLitGraph("/Users/dmuecher/iCloud/ShapeIt/ShapeIt1.0/Analysis/76Ge/NLD/ld_70_new_av.dat");
}

//fills the graph off transformed level densities using the literture values
void ShapeRhoCollector::FillTrafoGraph() {
   
    rhoTrafo = rho->rhoTrafoGraph(sett->lit_alpha,sett->lit_alpha - sett->lit_alpha_error[0], sett->lit_alpha + sett->lit_alpha_error[1] );

    rhoTrafo->SetMarkerColor(kBlack);
    rhoTrafo->SetLineColor(kBlack);
    rhoTrafo->SetFillColor(kWhite);

}


void ShapeRhoCollector::Draw() {

    //add band of discrete level fits
    if (discBand) {
        collGraph->Add(discBand,"3");
        collGraph->Add(discBandUpper,"L");
        collGraph->Add(discBandLower,"L");
        
    }
    //experimental error band from MC
    for (int i = 0 ; i < 6; i++) {
        if (ldmodel[i] && ldmodelChi2Display[i])
           // collGraph->Add((TGraph*)ldmodel[i]->expBand->Clone(),"3");
            collGraph->Add(ldmodel[i]->expBand,"3");
    }

  
    //include talys ldmodel graphs, if requested
    for (int i = 0 ; i < 6; i++) {
        if (ldmodel[i] && ldmodelDisplay[i]) {
            collGraph->Add(ldmodelGraph[i],"L");
        }
    }
   
    
    //add untransfomred level density
    if (rhoOrigDisplay) {
        TGraphErrors* origGraph =rho->rhoGraph;
        origGraph->SetFillColor(kAzure-9);
        origGraph->SetLineColor(kAzure-9);
        origGraph->SetTitle("^{76}Ge Spyrou 2014");
        collGraph->Add(rho->rhoGraph,"3");
    }
    
    //add literature Graphs
    if (litGraphDisplay) {
        for (int i = 0; i < litGraphs.size(); i++) {
            litGraphs[i]->SetMarkerColor(kRed);
            litGraphs[i]->SetMarkerStyle(24);
            litGraphs[i]->SetMarkerSize(1);
            litGraphs[i]->SetLineColor(kRed);
            litGraphs[i]->SetLineWidth(2);
            

            
            //collGraph->Add(litGraphs[i],"3");
            collGraph->Add(litGraphs[i],"P");
        }
    }
    
    //add rho data to Multigraph
    rhoTrafo->SetMarkerColor(kAzure+3);
    rhoTrafo->SetMarkerSize(1.6);
    rhoTrafo->SetLineColor(kAzure+3);
    rhoTrafo->SetLineWidth(2);

    collGraph->Add(rhoTrafo,"APE");
    
  
    
    //draw
    collGraph->GetYaxis()->SetRangeUser(0.5,1E4);
    collGraph->Draw("apl");
    collGraph->GetYaxis()->SetTitleOffset(1.4);
    collGraph->GetYaxis()->SetTitle("Level density #rho (E) (MeV^{-1})");
    collGraph->GetXaxis()->SetTitle("Energy (MeV)");
    collGraph->SetTitle("Level Density");

    //add discrete levels
    if (discreteLevel) {
        discreteLevel->SetFillColorAlpha(kAzure-9,0.4);
        discreteLevel->SetFillStyle(3002);
        discreteLevel->SetLineColorAlpha(kBlack,0.6);

        discreteLevel->Draw("same hist");
    }

}

//runs the Monte Carlo simulations for each defined ldmodel
void ShapeRhoCollector::ExpMC() {
    
    for (int i = 0 ; i < 6; i++) {
        if (ldmodel[i] && ldmodelChi2Display[i])
            ldmodel[i]->Chi2PartialLoopMC(2.2, 4);
    }
}

//calucates the difference of the discrete levels to the experimental levels for the set energy range
double ShapeRhoCollector::Chi2DiscreteExp(double eMin, double eMax, double norm) {
    
    //loop over discreteLevel histogram
    TAxis *xaxis = discreteLevel->GetXaxis();

    int XNum  = discreteLevel->GetNbinsX();
    double ene, disRho, chi2;
    chi2 = 0;
    for (int i = 1; i < XNum + 1; i++) {
         ene = xaxis->GetBinCenter(i);
        disRho = discreteLevel->GetBinContent(i);
        
        //energy cutoff
        if (ene >= eMin && ene <= eMax &&disRho > 0) {
            chi2 += TMath::Power((norm*rhoTrafo->Eval(ene)-disRho)/disRho,1);        }
    }
    return TMath::Abs(chi2);
}

//determines normalization of experimental level denisties
void ShapeRhoCollector::ExpNorm() {
    
    double chi2, chi2_min, norm_min;
    chi2_min = 1E5;
    for (double m_norm = 0.6; m_norm < 1.6; m_norm+=0.01 ) {
        chi2 = Chi2DiscreteExp(0.0, 3.4, m_norm);
        std::cout <<"chi2 = " <<chi2 <<" for norm factor of: "<<m_norm <<std::endl;

        if (chi2 < chi2_min) {
            chi2_min = chi2;
            norm_min = m_norm;
        }
    }
    std::cout <<"best fit gives chi2 = " <<chi2_min <<" for norm factor of: "<<norm_min <<std::endl;
    norm = norm_min;
    
    //set normalization for level denisty object rho and re-transform
    rho->rhoScaleTrafo = norm;
    FillTrafoGraph();
}

//adds a new literature graph to the plot
void ShapeRhoCollector::AddLitGraph(string litName) {
    litGraphs.push_back(new TGraphErrors(litName.c_str(),"%lg %lg %lg"));
    litGraphs[litGraphs.size()-1]->SetTitle("^{76}Ge Voinov 2019");
}

//fills the uncertainty band for fitting the discrete levels within the set ctable and ptable extremes
void ShapeRhoCollector::DiscreteBand() {
    
    double cTableExtreme = -0.7;
    int ldmodelNr = 5;
    double ctableBak = sett->cTable[ldmodelNr-1];
    ShapeTalys* discFit[15];
    ShapeMultiGraph* disc_sgraph = new ShapeMultiGraph();
    
    for (int i =0; i <15; i++) {
        cTableExtreme = -0.7 +( 0.1 * i );
        sett->cTable[ldmodelNr-1] = cTableExtreme;
        discFit[i] = new ShapeTalys(sett, rhoTrafo, ldmodelNr);
        discFitGraph[i] = new TGraph();
        
        double m_ptable = discFit[i]->PTableFromCTableDiscrete(cTableExtreme,2,3);
        std::cout <<"Best result for ctable = " <<cTableExtreme <<" : m_ptable: " <<m_ptable <<std::endl;
        discFit[i]->SetPCTable(m_ptable, cTableExtreme);
        discFitGraph[i] = discFit[i]->getDenPartialGraphTrans();
        discFitGraph[i]->SetLineWidth(4);
        discFitGraph[i]->SetLineColor(kBlack);
        disc_sgraph->Add(discFitGraph[i],"L");
        //disc_sgraph->Add(discFitGraph,"L");
        
        double m_ptable2 = discFit[i]->PTableFromCTableDiscrete(cTableExtreme,2,2.8);
        std::cout <<"Best result for ctable = " <<cTableExtreme <<" : m_ptable: " <<m_ptable2 <<std::endl;
        discFit[i]->SetPCTable(m_ptable2, cTableExtreme);
        discFitGraph[i] = discFit[i]->getDenPartialGraphTrans();
        discFitGraph[i]->SetLineWidth(4);
        discFitGraph[i]->SetLineColor(kBlack);
        //disc_sgraph->Add(discFitGraph[i],"L");

    }

    disc_sgraph->doFill(1,15);

    //discBand = (TGraphErrors*)disc_sgraph->fillGraph->Clone();
    discBand = new TGraphErrors();
    *discBand = *disc_sgraph->fillGraph;
    discBandUpper = new TGraph();
    *discBandUpper = *disc_sgraph->fillGraphUpper;
    discBandLower = new TGraph();
    *discBandLower = *disc_sgraph->fillGraphLower;
    
    discBand->SetTitle("fit to discrete exp. levels");
   
    discBand->SetFillColorAlpha(color->GetColor(244,207,223),1);
}

//loops through the settings file and fills the ldmodel graphs
void ShapeRhoCollector::FillTalys() {
    
    ldmodelCounter = 0;
    for (int i = 0 ; i < 6; i++) {
        if (sett->ldFileName[i]=="")
            continue;
        ldmodelCounter++;
        
        ldmodel[i] = new ShapeTalys(sett, rhoTrafo, i+1);
        
        ldmodelGraph[i] = ldmodel[i]->getDenPartialGraphTrans();
        
        ldmodelGraph[i]->SetTitle((talysNames[i]).c_str());
        ldmodelGraph[i]->SetLineColor(kBlack);
        ldmodelGraph[i]->SetLineWidth(3);
        switch(i) {
          case 3:
            ldmodelGraph[i]->SetLineStyle(2);
            break;
          case 4:
            ldmodelGraph[i]->SetLineStyle(1);
            break;
          case 5:
            ldmodelGraph[i]->SetLineStyle(3);
            break;
          
        }
        
        
    }
}

void ShapeRhoCollector::ReadDiscrete() {
    ifstream file(sett->discreteLevelFile);
    if (!file)
    {
        cerr << "cannot read the discrete level file!"
        << strerror(errno) << endl;
        return;
    }
    
    double e, disc;
    delete gROOT->FindObject("discreteLevel");
    std::vector <double> ene;
    std::vector <double> discLev;
    
    //read the discrete level file data
    while (file >> e >> disc) {
        ene.push_back(e);
        discLev.push_back(disc);
    }
    
    //find the maximum energy
    double discreteMax = *max_element(ene.begin(), ene.end());
    
    std::cout <<"discreteMax: " <<discreteMax<<std::endl;
    //define histogram
    discreteLevel = new TH1F("discreteLevel","discrete levels",((1000.*(discreteMax)/sett->discreteBins)),0,discreteMax);
    
    //fill histogram
    for (int i = 0; i < ene.size(); i++)
        discreteLevel->Fill(ene[i],discLev[i]);
    
    discreteLevel->SetLineColor(kBlack);
    
}

         
          
          
