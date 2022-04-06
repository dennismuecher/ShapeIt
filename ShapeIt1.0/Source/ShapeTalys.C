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

#include "../Include/ShapeTalys.h"
#include "../Include/ShapeMultiGraph.h"

//consstructor

ShapeTalys::ShapeTalys(ShapeSetting* p_sett, TGraphAsymmErrors* p_rhoGraph, int p_levelmodelNr) {
    sett = p_sett;
    levelmodelNr = p_levelmodelNr;
    //if levelmodelNr == 0, only the discrete level File is requested
    if (levelmodelNr >0) {
        parityFlag = sett->parityFlag[levelmodelNr-1];
        rhoGraph = p_rhoGraph;
        format = sett->formatFlag[levelmodelNr-1];
        talysOutFile = sett->ldFileName[levelmodelNr-1];
        ptable = sett->pTable[levelmodelNr-1];
        ctable =sett->cTable[levelmodelNr-1];
        NewReadTree();
        SetPCTable();
    }
    discreteLevelFile = sett->discreteLevelFile;
    ReadDiscrete();
}

//calculates a new MC representation of the level density graph
TGraph* ShapeTalys::MCRhoGraph() {
    
    double e, r, dr;
    TGraph* mcGraph = new TGraph();
    TRandom3 *ran = new TRandom3(0);

    for (int i =0; i < rhoGraph->GetN(); i++) {
        e = rhoGraph->GetPointX(i);
        r = rhoGraph->GetPointY(i);
        dr = (rhoGraph->GetEYlow()[i] + rhoGraph->GetEYhigh()[i]) / 2.;
        mcGraph->SetPoint(mcGraph->GetN(), e, ran->Gaus(r,dr));
    }
    rhoGraphMC = mcGraph;
    return mcGraph;
}

//calculates Chi2 between actual theoretical level density for the partial graph and the experimental MC values
double ShapeTalys::GetChi2PartialMC(double lower_ene, double higher_ene) {
    
    double chi2 = 0;
    double e, rho, rho_theo;
    
    for (int i =0; i < rhoGraphMC->GetN(); i++) {
        e = rhoGraphMC->GetPointX(i);
        rho = rhoGraphMC->GetPointY(i);
        if (e > lower_ene && e < higher_ene ) {
            rho_theo = denPartialGraphTrans->Eval(e);
            //chi2 += TMath::Power((rho-rho_theo),2)/rho_theo;
            chi2 += TMath::Power((rho-rho_theo)/rho_theo,2);
        }
    }
    return chi2;
}

//calculates Chi2 between actual theoretical level density for the partial graph and the experimental values stored in rho
double ShapeTalys::GetChi2Partial(double lower_ene, double higher_ene) {
    
    double chi2 = 0;
    double e, rho, rho_theo;
    int n =0;
    for (int i =0; i < rhoGraph->GetN(); i++) {
        e = rhoGraph->GetPointX(i);
        rho = rhoGraph->GetPointY(i);
        double drhoh = rhoGraph->GetErrorYhigh(i);
        double drhol = rhoGraph->GetErrorYlow(i);
        double drho = (drhoh + drhol)/2;
        
        if (e > lower_ene && e < higher_ene ) {
            rho_theo = denPartialGraphTrans->Eval(e);
            if (drho >0) {
                chi2 += TMath::Power((rho-rho_theo)/drho,2);
                //chi2 += TMath::Power((rho-rho_theo),2)/rho_theo;
                n++;
            }
        }
        
    }
    nOfDegFreedom = n;
   //return ( TMath::Prob(chi2,n) );
    if (n > 1)
        return chi2;
    else
        return (-1);
}

//calculates Chi2 between actual theoretical level density for the partial graph and the discrete levels
double ShapeTalys::GetChi2Discrete(double lower_ene, double higher_ene) {
    
    double chi2 = 0;
    double e, rho, rho_theo;
    
    //number of degrees of freedom
    int n =0;

    for (int i = 1; i < discreteHist->GetNbinsX() + 1; i++) {
        e = discreteHist->GetXaxis()->GetBinCenter(i);
        rho = discreteHist->GetBinContent(i);
        
        if (e > lower_ene && e < higher_ene ) {
            rho_theo = denPartialGraphTrans->Eval(e);
            chi2 += TMath::Power(rho-rho_theo,2)/rho_theo;
            //chi2 += TMath::Power((rho-rho_theo)/rho_theo,2);
            n++;
        }
    }
    return chi2/n;
}

//determines a best fit to the discrete levels for ptable using a fixed ctable
double ShapeTalys::PTableFromCTableDiscrete(double m_ctable, double lower_ene, double higher_ene) {
    
    double min_ptable = -4;
    double max_ptable = +4;
    double m_ptable;
    double best_ptable;
    double red_chi2;
    double min_chi2 = 1E5;
    for (int i = 0; i < 100; i++) {
        m_ptable = i* (max_ptable - min_ptable)/100 + min_ptable;
        SetPCTable(m_ptable, m_ctable);
        red_chi2 = GetChi2Discrete(lower_ene,higher_ene);

        if (red_chi2 < min_chi2) {
            min_chi2 = red_chi2;
            best_ptable = m_ptable;
        }
            
    }
    return best_ptable;
}

void ShapeTalys::Chi2PartialLoopMC(double lower_ene, double higher_ene) {

    bool showPlot = false;
    TMultiGraph* m_graph = new TMultiGraph();
    //TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,600,480);
    gPad->SetLogy();

    if (showPlot) {
        //m_graph->Add(rhoGraph, "APE");
        rhoGraph->Draw("APE");
    }
    double min_ptable = 0.8;
    double max_ptable = +2;
    
    double min_ctable = 0.5;
    double max_ctable = 1.4;
    
    double ptable_mean = 1.36;
    double ptable_sigma = 0.15;
    
    double ctable_mean = 0.85;
    double ctable_sigma = 0.13;
    
    ShapeMultiGraph* s_graph = new ShapeMultiGraph();
    
    delete gROOT->FindObject("bestFitMC");
    bestFitMC = new TH2D("bestFitMC","ptable vs ctable from MC simulation",40,min_ptable,max_ptable,40,min_ctable,max_ctable);
    
    int nOfIter = 51;
    int nOfGraphs = 0;
    
    TGraph* m_rhoGraph = new TGraph();
    double ptableAccepted[100];
    double ctableAccepted[100];
    for (int mc_run = 0; mc_run < nOfIter; mc_run++) {
        //get new representation of level densities
        m_rhoGraph = MCRhoGraph();
        
        if (showPlot) {
            rhoGraph->Draw("APE");
            m_rhoGraph->Draw(" * same");
            m_graph->Add(rhoGraphMC,"AP");
        }
        //the best chi2 result
        chi2_min = 1E5;
        double chi2;
        
       
        
        for (double p = min_ptable; p < max_ptable; p +=0.025) {
            
            for (double c = min_ctable; c < max_ctable; c +=0.025) {

                SetPCTable(p,c);
                chi2 = GetChi2PartialMC(lower_ene, higher_ene);
                if (chi2 < chi2_min) {
                    chi2_min = chi2;
                    ptablePartialMC = p;
                    ctablePartialMC = c;
                }
            }
           
        }
        if (mc_run % 10 ==0) {
            std::cout <<"MC Iteration: " <<mc_run <<std::endl;
            std::cout <<"accepted fits in 1-sigma: " <<nOfGraphs <<std::endl;
        }
        
        bestFitMC->Fill(ptablePartialMC,ctablePartialMC);
        
        //cut on 1-sigma accetance for ptable and ctable; the values for the mean and sigma of ptable and ctable can be determined through the MC process
        if (TMath::Abs(ptablePartialMC - ptable_mean) > ptable_sigma)
           continue;
        if (TMath::Abs(ctablePartialMC - ctable_mean) > ctable_sigma)
            continue;
        nOfGraphs++;
        if (nOfGraphs < 100) {
            ptableAccepted[nOfGraphs] =ptablePartialMC;
            ctableAccepted[nOfGraphs] =ctablePartialMC;

        }
        std::cout << std::setprecision(3)<<ptablePartialMC << " " << ctablePartialMC <<std::endl;
        
        SetPCTable(ptablePartialMC,ctablePartialMC);
        //denPartialGraphTrans->Draw("same");
        m_graph->Add(denPartialGraphTrans,"L");
        denPartialGraphTrans->Draw("L same");
        s_graph->Add((TGraph*)denPartialGraphTrans->Clone(),"L");

        // make sure to update display
        //m_graph->Draw("apl");
        //gPad->Modified();
        //gPad->Update();
        //m_graph->GetHistogram()->GetYaxis()->SetRangeUser(1E-2,1E3);
        //m_graph->GetHistogram()->GetXaxis()->SetRangeUser(1,7);

        //gSystem->ProcessEvents(); gSystem->ProcessEvents();
        //gSystem->Sleep(10);
        
    }
    TColor* color;
    std::cout <<"ptable = {";
    for (int i = 0; i < nOfGraphs; i++) {
        std::cout << std::setprecision(3) << ptableAccepted[i] <<", ";
    }
    std::cout <<"};"<<std::endl;
    std::cout <<"ctable = {";

    for (int i = 0; i < nOfGraphs; i++) {
        std::cout << std::setprecision(3) << ctableAccepted[i] <<", ";
    }
    s_graph->doFill(1,nOfGraphs);
    expBand = (TGraphErrors*)s_graph->fillGraph->Clone();
    //expBand->SetFillColorAlpha(kMagenta, 0.8);
    expBand->SetFillColorAlpha(color->GetColor(87,132,186),0.7);
    //expBand->SetFillStyle(3010);
    string s = talysNames[levelmodelNr-1]+" (1#sigma fit to data)";
    expBand->SetTitle(s.c_str());

}

void ShapeTalys::Chi2PartialLoop(double lower_ene, double higher_ene) {
    
    TCanvas *canv;
    bool plot_fit = false;

    if (plot_fit) {
        canv = new TCanvas("c", "c", 300, 300);
        rhoGraph->Draw("APC");
    }
    double min_ptable = -1.5;
    double max_ptable = +1.5;
    
    double min_ctable = -1.5;
    double max_ctable = 1.5;
    delete gROOT->FindObject("chi2Fit");
    TH2D* t = new TH2D("chi2Fit","Chi2 value for partial level density fit",200,min_ptable,max_ptable,200,min_ctable,max_ctable);
    
    //the best chi2 result
    chi2_min = 1E5;
    //degrees of freedom for best fit
    int n;
    for (double p = min_ptable; p < max_ptable; p +=0.05) {
       
        double chi2;
        double c_min;

        for (double c = min_ctable; c < max_ctable; c +=0.05) {
            SetPCTable(p,c);
            chi2 = GetChi2Partial(lower_ene, higher_ene);
            if (chi2 < chi2_min) {
                chi2_min = chi2;
                c_min = c;
                ptablePartial = p;
                ctablePartial = c;
                n = nOfDegFreedom;
            }
            
            //fill graph of chi2 values if probability for larger chi2 above critical value
            if (TMath::Prob(chi2,nOfDegFreedom) >= 0.10) {
               t->Fill(p,c,chi2);
                //denPartialGraphTrans->Clone()->Draw("same");
            }
        }
        
        //display best fit for this loop
        if (plot_fit) {
            SetPCTable(p,c_min);
            //rhoGraph->Draw("APC");
            denPartialGraphTrans->Draw("same");
            gPad->SetLogy();
            gSystem->ProcessEvents(); gSystem->ProcessEvents();
            gSystem->Sleep(10);
            canv->Modified();canv->Update();
        }
        //std::cout<< "Minimum: " << chi2_min <<" "<< p <<" "<<" " << c_min<<std::endl;
    }
   
    if (plot_fit)
        t->Draw("colz");
    std::cout <<"Best Fit values for ldmodel "<< levelmodelNr << ": ptable =" <<ptablePartial <<"  ctable ="<<ctablePartial<<" chi2: " <<chi2_min << std::endl;
}

//Sets ptable and ctable to best fest results of model to partial level density; requires Chi2PartialLoop to run before
void ShapeTalys::BestFitPartial()
{
    if (ptablePartial == 0 && ctablePartial == 0)
        std::cout <<"Cannot deliver best fit to partial level density; run chi2partialloop, first!"<<std::endl;
    else
        SetPCTable(ptablePartial, ctablePartial);
}

//sets ptable to new value and transforms all level density graphs
void ShapeTalys::SetPCTable()
{
    SetPCTable(ptable, ctable);
}


//sets ptable to new value and transforms all level density graphs
void ShapeTalys::SetPCTable(double m_ptable, double m_ctable)
{
    double e,f;
    //reset all transformed graphs
    denTotGraphTrans->Set(0);
    denPartialGraphTrans->Set(0);
    for (int j = 0; j <nOfSpins; j++)
        denSpinGraphTrans[j]->Set(0);
    
    //fill transformed graphs
    for (int i = 0; i < denTotGraph->GetN(); i++) {
        
        //shift of energy
        e = denTotGraph->GetPointX(i) + m_ptable;
        
        //energy can be smaller zero at this point
        if ( (e >= 0.5) && (e <=10) )
            f = TMath::Exp((m_ctable)*TMath::Sqrt(denTotGraph->GetPointX(i) ));
        else
            continue;
        
        int k = denTotGraphTrans->GetN();
        denTotGraphTrans->SetPoint(k, e, f * denTotGraph->GetPointY(i));
        denPartialGraphTrans->SetPoint(k, e, f * denPartialGraph->GetPointY(i));

         for (int j = 0; j <nOfSpins; j++) {
             denSpinGraphTrans[j]->SetPoint(k, e, f * denSpinGraph[j]->GetPointY(i));
         }
    }
    ptable = m_ptable;
    ctable = m_ctable;
}


int ShapeTalys::NewReadTree()
{
    ifstream file(talysOutFile.c_str());
    if (!file)
    {
        cerr << "cannot read the file"
        << strerror(errno) << endl;
        return -1;
    }
    
    float e, denTot;
    float den[nOfSpins];
    float a, sigma;
    if (format == 0) {
        while (file >> e >> denTot >>den[0] >>den[1] >>den[2] >>den[3]>>den[4] >>den[5]>>den[6] >>den[7] >>den[8] ) {
            p_energy.push_back(e);
            p_densityTot.push_back(denTot);
            for (int i =0; i < nOfSpins; i++)
                p_densitySpin[i].push_back(den[i]);
        }
    }
    else {
        while (file >> e >> a >> sigma >> denTot >>den[0] >>den[1] >>den[2] >>den[3]>>den[4] >>den[5]>>den[6] >>den[7] >>den[8] ) {
            p_energy.push_back(e);
            p_densityTot.push_back(denTot);
            for (int i =0; i < nOfSpins; i++)
                p_densitySpin[i].push_back(den[i]);
        }
    }
        
    // in case of files with two parities, add level densities, accordingly; otherwise multiply with 2
    if (parityFlag) {
        int n =(int) p_energy.size()/2;
        for (int i =0; i < n; i++) {
            densityTot.push_back(p_densityTot[i]+p_densityTot[i+n]);
            energy.push_back(p_energy[i]);
            double partial = 0;
            for (int j =0; j < nOfSpins; j++) {
                densitySpin[j].push_back(p_densitySpin[j][i]+p_densitySpin[j][i+n]);
                if (j >=sett->spinLow && j <=sett->spinHigh)
                    partial +=densitySpin[j][i];
            }
            densityPartial.push_back(partial);
        }
    }
    else {
        int n =(int) p_energy.size();
        for (int i =0; i < n; i++) {
            densityTot.push_back(2*p_densityTot[i]);
            energy.push_back(p_energy[i]);
            double partial = 0;
            for (int j =0; j < nOfSpins; j++) {
                densitySpin[j].push_back(2*p_densitySpin[j][i]);
                if (j >=sett->spinLow && j <=sett->spinHigh)
                    partial +=densitySpin[j][i];
            }
            densityPartial.push_back(partial);
        }
    }
    
    //create TGraphs
    int n =densityTot.size();
    denTotGraph = new TGraph(n, &energy[0],&densityTot[0]);
    denTotGraphTrans = new TGraph(n, &energy[0],&densityTot[0]);

    denPartialGraph = new TGraph(n, &energy[0],&densityPartial[0]);
    denPartialGraphTrans = new TGraph(n, &energy[0],&densityPartial[0]);

    for (int i = 0; i <nOfSpins; i++) {
        denSpinGraph[i] = new TGraph(n, &energy[0],&densitySpin[i][0]);
        denSpinGraphTrans[i] = new TGraph(n, &energy[0],&densitySpin[i][0]);

        denSpinGraph[i]->SetLineColor(i);
    }
    
    return 0;
}

void ShapeTalys::ReadDiscrete() {
    ifstream file(discreteLevelFile);
    if (!file)
    {
        cerr << "cannot read the discrete level file!"
        << strerror(errno) << endl;
    }
    
    float e, disc;
    delete gROOT->FindObject("disLevel");
    std::vector <float> ene;
    std::vector <float> discLev;
    
    //read the discrete level file data
    while (file >> e >> disc) {
        ene.push_back(e);
        discLev.push_back(disc);
    }
    
    //find the maximum energy
    discreteMax = *max_element(ene.begin(), ene.end());
    
    //define histogram
    discreteHist = new TH1F("disLevel","discrete levels",((1000.*(discreteMax)/sett->discreteBins)),0,discreteMax);
    
    //fill histogram
    for (int i = 0; i < ene.size(); i++)
        discreteHist->Fill(ene[i],discLev[i]);
    
    discreteHist->SetLineColor(kBlack);
}

TGraph* ShapeTalys::getDenPartialGraphTrans() {
    TGraph* test = new TGraph();
    *test = *denPartialGraphTrans;
    return test;
    
}

