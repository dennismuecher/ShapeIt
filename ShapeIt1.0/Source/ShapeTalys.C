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

//consstructor
ShapeTalys::ShapeTalys(ShapeSetting* p_sett, std::string p_talysOutFile, TGraphAsymmErrors* p_rhoGraph, bool p_parityFlag, bool p_format, double p_ptable, double p_ctable) {
    sett = p_sett;
    parityFlag = p_parityFlag;
    rhoGraph = p_rhoGraph;
    format = p_format;
    talysOutFile = p_talysOutFile;
    discreteLevelFile = sett->discreteLevelFile;
    ptable = p_ptable;
    ctable = p_ctable;
    NewReadTree();
    SetPCTable();
    ReadDiscrete();
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

/*//calculates Chi2 between actual theoretical level density for the partial graph and the discrete levels
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
            n++;
        }
    }
    return chi2/n;
}

*/
void ShapeTalys::Chi2PartialLoop(double lower_ene, double higher_ene) {
    
    TCanvas *canv;
    bool plot_fit = false;
    //rhoGraph->Draw("APC");

    if (plot_fit) {
        canv = new TCanvas("c", "c", 300, 300);
        rhoGraph->Draw("APC");
    }
    double min_ptable = -1.5;
    double max_ptable = +1.5;
    
    double min_ctable = -1.5;
    double max_ctable = +1.5;
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
    //rhoGraph->SetLineColor(kBlack);
     // rhoGraph->SetMarkerColor(kBlack);
    //rhoGraph->Draw("same");
    if (plot_fit)
        t->Draw("colz");
    std::cout <<"Best Fit values: " <<ptablePartial <<" "<<ctablePartial<<" " <<chi2_min<<" "<<n<< std::endl;
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
        if (e >= 0)
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
                if (j >=spinLow && j <=spinHigh)
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
                if (j >=spinLow && j <=spinHigh)
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
    discreteHist = new TH1F("disLevel","discrete Levels",int(8000/sett->discreteBins),-1,7);
        while (file >> e >> disc) {
            discreteHist->Fill(e,disc);
        }
    discreteHist->SetLineColor(kBlack);
}


