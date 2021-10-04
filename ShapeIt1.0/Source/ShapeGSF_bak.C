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
    litGraph = new TGraphErrors();
    mc_litGraph = new TGraph();
    mc_litGraphAll = new TMultiGraph();
    levGraphAll = new TGraphErrors();
    mc_levGraphAll = new TGraph();
   

    levGraphSmooth = new TGraphAsymmErrors();
    if (m_sett->doOslo) {
        ReadLit();
    }
}

//reads the literature values for gSF
void ShapeGSF::ReadLit()
{
    //read file data into gSF
    if (m_sett->osloFileName == "") {
        std::cout << "No Literature File loaded!"<<std::endl;
        return NULL;
    }
    
    if (m_sett->verbose)
        std::cout <<"\nReading OSLO DATA... " <<endl;
    ifstream inp;
    inp.open(m_sett->osloFileName.c_str());
    if (inp.is_open() ) {
        
        double e_gamma;
        double oslo_gSF_high, oslo_gSF_low;
        litGraph->Set(0);
        while ( !inp.eof() ) {
            inp >> e_gamma >> oslo_gSF_high >>oslo_gSF_low;
            litGraph->SetPoint(litGraph->GetN(), e_gamma,
                                    ( oslo_gSF_high + oslo_gSF_low ) / 2 );
            
            litGraph->SetPointError(litGraph->GetN() - 1, 0,
                                    ( oslo_gSF_high - oslo_gSF_low ) / 2 );
        }
        
        litGraph->SetFillColor(4);
        litGraph->SetFillStyle(3010);
        if (m_sett->verbose)
            litGraph->Print();
    }
}

// filles mc_litGraph with new values from litGraph using MC smearing
void ShapeGSF::doMCGraph()
{
    TRandom3 s;
    mc_litGraph->Set(0);
    for (int i = 0; i < litGraph->GetN(); i++)
    {
        //important: reset transformation as mc_litGraph is regenerated using litGraph, which has not been transformed
        m_B = 1;
        m_alpha = 0;
        mc_litGraph->SetPoint(i, litGraph->GetX()[i], s.Gaus(litGraph->GetY()[i], litGraph->GetEY()[i] ) );
        //mc_litGraph->SetPoint(i, litGraph->GetX()[i], litGraph->GetY()[i] );
    }
    //litGraph->Print();
    //mc_litGraph->Print();
    mc_litGraphAll->Add(mc_litGraph);
}


//calcualtions a normalization factor to match T1 and T2
double ShapeGSF::Norm(TGraphErrors* T1, TGraphErrors* T2 ) {
    double c1 = 0, c2 = 0;
    std::vector<double> a;
    std::vector<double> da;
    std::vector<double> b;
    for (int i = 0; i < T1->GetN(); i++) {
        a.push_back( T1->GetY()[i] );
        da.push_back( T1->GetEY()[i] );
        b.push_back( T2->Eval( T1->GetX()[i] ) );
    }
    for (int i = 0; i < T1->GetN(); i++) {
        //in case of MC mode, there are no error bars
        if (m_sett->doMC) {
            c1 += ( a[i] * b[i] );
            c2 += ( a[i] * a[i] );
        }
        else if ( da[i] !=0 ) {
          c1 += ( a[i] * b[i] / da[i] );
          c2 += ( a[i] * a[i] / da[i] );
        }
        else {
            std::cout <<"Error is Norm(): zero error value detected" <<std::endl;
            return (0);
        }
    }
    if (c2 !=0)
        return (c1 / c2);
    else {
        std::cout <<"Error is Norm(): zero c2 value detected" <<std::endl;
        return (0);
    }
}

//returns the chi2 for the MC mode, comparing the smooth gSF to the Oslo literature data

double ShapeGSF::mc_getChi2(){
    
    double chi2 = 0;
    for (int i = 0; i < levGraphAll->GetN(); i++) {
        double x = levGraphAll->GetX()[i];
        double y = levGraphAll->GetY()[i];
        double y_oslo = mc_litGraph->Eval(x);
        if (y > 0) {
            chi2 += TMath::Power( (y_oslo - y) /y, 2);
        }
        else {
            std::cout <<"Error in getChi2: zero or negative value for y_oslo error for point " <<i<<" ! Skipping this data point" <<std::endl;
            continue;
        }
    }
    return (chi2);
}


//returns the chi2 for the levGraphAll graph, comparing the smooth gSF to the Oslo literature data

double ShapeGSF::getChi2All(){
    
    double chi2 = 0;
    for (int i = 0; i < levGraphAll->GetN(); i++) {
        double x = levGraphAll->GetX()[i];
        double y = levGraphAll->GetY()[i];
        double dyh = levGraphAll->GetErrorYhigh(i);
        double dyl = levGraphAll->GetErrorYlow(i);
        double dy = (dyh + dyl)/2;
        double y_oslo = litGraph->Eval(x);
        
        if (dy == 0) {
            std::cout <<"Error in getChi2: zero value for gSF error for point " <<i<<" ! Skipping this data point" <<std::endl;
            continue;
        }
        else
            chi2 += TMath::Power( ( y_oslo - y ) / dy, 2);
    }
    return (chi2);
}

//returns the chi2 for the levGraphSmooth graph, comparing the smooth gSF to the Oslo literature data

double ShapeGSF::getChi2Smooth(){
    
    double chi2 = 0;
    int n = levGraphSmooth->GetN();
    for (int i = 0; i < n; i++) {
        double x = levGraphSmooth->GetX()[i];
        double y = levGraphSmooth->GetY()[i];
        double dyh = levGraphSmooth->GetErrorYhigh(i);
        double dyl = levGraphSmooth->GetErrorYlow(i);
        double dy = (dyh + dyl)/2;
        double y_oslo = litGraph->Eval(x);
        
        if (dy == 0) {
            std::cout <<"Error in getChi2: zero value for gSF error for point " <<i<<" ! Skipping this data point" <<std::endl;
            continue;
        }
        else
            chi2 += TMath::Power( ( y_oslo - y ) / dy, 2);
    }
    //calculating the normal-distributed chi2 value, following this website:
    //https://www.phys.hawaii.edu/~varner/PHYS305-Spr12/DataFitting.html
    n = n-1; //number of degrees of freedom; assumes that levGraphSmooth points are statistically independent
    double chi2_distr = TMath::Sqrt(2*chi2) - TMath::Sqrt(2*n -1);
    std::cout <<"Normal distributed chi2 value: " <<chi2_distr<< " for chi2 " <<chi2 <<" for alpha: "<< m_sett->lit_alpha <<std::endl;
    // probability that an observed Chi-squared exceeds the value chi2 by chance, even for a correct model
    double chi2_prob = TMath::Prob(chi2,n);
    std::cout <<"chi2 probability: " <<chi2_prob << std::endl;
    return (chi2);
}

double ShapeGSF::getChi2(){
    if (m_sett->doMC)
        return mc_getChi2();
    else if (m_sett->displayAvg)
        return getChi2Smooth();
    else
        return getChi2All();
}
 

//fills the actual vector of levGraph_1 and levGraph_2 with values; also fills levGraph
void ShapeGSF::FillgSF() {
   
    m_matrix->Integrate();
    m_matrix->IntegrateBg();
    m_matrix->IntegrateSquare();
    m_matrix->IntegrateCube();
   
    //if autofit or background subtraction is set, perform autofit using Gauss
    if (m_sett->mode == 2 || m_sett->doBackground)
        m_matrix->FitIntegral();

    double elevel1 = ( m_sett->levEne[0] + m_sett->levEne[1] ) / 2;
    double elevel2 = ( m_sett->levEne[2] + m_sett->levEne[3] ) / 2;
    double egamma1, egamma2;
    double gSF1, dgSF1, gSF2, dgSF2;
    
    //create new vector element; j is the the index of the first empty container, i.e. "Update" will fill a new vector container
    
    if (m_sett->doMC) {
        levGraph.clear();
        levGraph_1.clear();
        levGraph_2.clear();
        levGraphAll->Set(0);
        levGraphSmooth->Set(0);
        //get a new Monte Carlo version of the literature values
        doMCGraph();
    }
    
    levGraph.push_back(new TGraphErrors());
    levGraph_1.push_back(new TGraphErrors());
    levGraph_2.push_back(new TGraphErrors());
    
    int j = levGraph_1.size() - 1;
    
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
        
        //in case of MC mode, shuffle gSF according to the gauss distribution
        if (m_sett->doMC) {
            //gaus smear according to fit eerror bars
            gSF1 = r.Gaus(gSF1,dgSF1);
            gSF2 = r.Gaus(gSF2,dgSF2);
            //now set error bars to zero, they are not used for anything anymore
            dgSF1 = 0;
            dgSF2 = 0;
        }

        
        //filling results into levGraph vector;
        levGraph_1[j]->SetPoint(levGraph_1[j]->GetN(), egamma1, gSF1);
        levGraph_1[j]->SetPointError(levGraph_1[j]->GetN()-1, 0, dgSF1);
        levGraph_2[j]->SetPoint(levGraph_2[j]->GetN(), egamma2, gSF2);
        levGraph_2[j]->SetPointError(levGraph_2[j]->GetN()-1, 0, dgSF2);
    
        if (m_sett->verbose) {
            std::cout <<"Bin: " <<i+1<<std::endl;
            std::cout <<"energies: " <<egamma1 << " " << egamma2 <<std::endl;
            std::cout <<"gSF1: " <<gSF1 << "+- " << dgSF1 <<std::endl;
            std::cout <<"gSF2: " <<gSF2 << "+- " << dgSF2 <<std::endl;
        }
    }
    
    //scale gSF values according to user input
    Scale(j, m_sett->gSF_norm);
    scale_bak = m_sett->gSF_norm;
    //run interpolation
    if (m_sett->doInterpol)
        DoInterpol();
}

//merges levGraph1 and levGraph2 into one graph and sorts by energy; used for normalization of all graphs
void ShapeGSF::Merge() {

    for (int i = 0; i < levGraph_1.size(); i++) {
        TObjArray *mArray = new TObjArray();
        mArray->Add(levGraph_1[i]);
        mArray->Add(levGraph_2[i]);
        levGraph[i]->Set(0);
        levGraph[i]->Merge(mArray);
        levGraph[i]->Sort();
    }
}


//merge all lelvGraphs into single graph "levGraphAll"
void ShapeGSF::MergeAll() {
    TObjArray *mArray = new TObjArray();
    for (int i = 0; i < levGraph.size(); i++) mArray->Add(levGraph[i]);
    levGraphAll->Set(0);
    levGraphAll->Merge(mArray);
    levGraphAll->Sort();
}

/*
//merge all mc_lelvGraphs into single graph "mc_levGraphAll"
void ShapeGSF::mc_MergeAll() {
    TObjArray *mArray = new TObjArray();
    for (int i = 0; i < mc_levGraph.size(); i++) mArray->Add(mc_levGraph[i]);
    mc_levGraphAll->Set(0);
    mc_levGraphAll->Merge(mArray);
    mc_levGraphAll->Sort();
}*/

//returns a graph of gSF based on levGraphAll, but with averaged gSF values according to the bin size res

void ShapeGSF::Smooth(int res) {
    levGraphSmooth->Set(0);
    // if res == 0 use exi_size[0] (user input) for the bin size
    if (res ==0)
        res = m_sett->exi_size[0];
    auto nPoints = levGraphAll->GetN();
    double m_e = levGraphAll->GetX()[0] + res;
    std::vector <int> i_bin;          //the point number limits of levGraphAll belonging to the bins with size res
    i_bin.push_back(0);             //first bin always starts at zero
    for(int i =1; i < nPoints; i++) {
        if (levGraphAll->GetX()[i] > m_e) {
            i_bin.push_back(i);
            m_e += res;
            i--;
        }
    }
    
    //now calculate the average gSF values for each bin
    
    for (int i = 0 ; i < i_bin.size() -1; i++) {
        double x = 0, y = 0;
        double dy = 0;
        for (int j = i_bin[i]; j < i_bin[i+1]; j++) {
            x += levGraphAll->GetX()[j];
            y += levGraphAll->GetY()[j];
            dy += TMath::Power(levGraphAll->GetEY()[j] , 2);
            
        }
        int nOfP = i_bin[i+1] - i_bin[i];     //number of points in this bin
        if (nOfP ==0) continue;
        x = x / nOfP;
        y = y / nOfP;
        //this error is the error of the mean value, based on the propagation of the indiviual errors, only; this does not take into account the fluctuations in the data points (i.e. the standard deviation)
        dy = TMath::Sqrt(dy) / nOfP;
        levGraphSmooth->SetPoint(levGraphSmooth->GetN(),x,y);
        
       /*
        double st_dev = 0;
        for (int j = i_bin[i]; j < i_bin[i+1]; j++) {
            st_dev += TMath::Power(levGraphAll->GetY()[j] - y, 2);
        }
        st_dev = TMath::Sqrt(st_dev / (nOfP -1) );
        */
        //calculate max errors including errors of individual data points

        double dyl=0,dyh=0;
        for (int j = i_bin[i]; j < i_bin[i+1]; j++) {
                //difference of average to upper error bar
                double m_dyh = levGraphAll->GetY()[j] + levGraphAll->GetEY()[j] - y;
                if (m_dyh > dyh)
                    dyh = m_dyh;
                //difference of average to lower error bar
                double m_dyl = y - levGraphAll->GetY()[j] + levGraphAll->GetEY()[j];
                if (m_dyl > dyl)
                dyl = m_dyl;
        }
        levGraphSmooth->SetPointError(levGraphSmooth->GetN()-1, 0, 0, dyl, dyh);
        //some old stuff
        
       /* int upper = 0, lower = 0;
        double dyl=0,dyh=0;
        for (int j = i_bin[i]; j < i_bin[i+1]; j++) {
            //if ( levGraphAll->GetY()[j] > y) {
                dyh += TMath::Power(levGraphAll->GetY()[j]
                                +levGraphAll->GetEY()[j]
                                -y, 2);
                upper++;
            //}
            //else {
                dyl += TMath::Power(levGraphAll->GetY()[j]
                                -levGraphAll->GetEY()[j]
                                -y, 2);
                lower++;
            //}
            
        }
        
        dyl = TMath::Sqrt(dyl /lower );
        
        if (dyh == 0) {
            dyh = dyl;
        }
        else {
            dyh = TMath::Sqrt(dyh / upper );
        }
         
        levGraphSmooth->SetPointError(levGraphSmooth->GetN()-1,
                                      0, 0, dyl, dyh);*/
    }
        
    levGraphSmooth->SetMarkerStyle(22);
    levGraphSmooth->SetMarkerSize(2);
    levGraphSmooth->SetMarkerColor(1);
}


//creates the graph of all gSF data; for this, first interpolate the individual runs, if requested by user; then add them to a multigraph
TMultiGraph* ShapeGSF::getMultGraph() {
    
    //define Multigraph
    TMultiGraph* multGraph = new TMultiGraph();
    multGraph->SetTitle("Gamma Ray Strength Function from Shape Method; E_{#gamma} (keV); f(E_{#gamma} (MeV^{-3})" );
    
    //add literature values to graph; also applies transformation
    //if (m_sett->doOslo && !m_sett->doMC) {
        Transform(m_sett->lit_norm, m_sett->lit_alpha);
        multGraph->Add(litGraph,"3A");
   // }
   /* else if (m_sett->doOslo && m_sett->doMC) {
        Transform(m_sett->lit_norm, m_sett->lit_alpha);
        multGraph->Add(mc_litGraph);
    }*/
    //merge gSF data into levGraph
    Merge();
    
    //apply autoscale, if requested
    
    if (m_sett->doAutoScale) {
        m_sett->gSF_norm *= Norm(levGraph[0], litGraph );
    }
    
    //apply user (auto) scaling factor
    ScaleAll(m_sett->gSF_norm);
    
    //merge again to refresh gSF data
    Merge();
    
    for (int i = 0; i < levGraph_1.size(); i++) {
   
        levGraph_1[i]->SetMarkerStyle(22);
        levGraph_1[i]->SetMarkerSize(2);
        levGraph_1[i]->SetMarkerColor(6);
        levGraph_2[i]->SetMarkerStyle(22);
        levGraph_2[i]->SetMarkerSize(2);
     
        if (m_sett->colour)
            levGraph_2[i]->SetMarkerColor(7);
        else
            levGraph_2[i]->SetMarkerColor(6);
        
        //do the normalization to the first graph
        if (i != 0) {
            Scale(i, Norm(levGraph[i], levGraph[0] ) );
        }
        if (m_sett->displaySingle) {
            multGraph->Add(levGraph_1[i],"P");
            multGraph->Add(levGraph_2[i],"P");
        }
    }
    Merge();
    MergeAll();
    //add smoothed graph, if requested
    if (m_sett->displayAvg) {
        Smooth(0);
        multGraph->Add(getSmoothGraph(),"P");
    }
    return ( multGraph );
}

double ShapeGSF::Slope(int i) {
    int j = levGraph_1.size()-1;
    return ( levGraph_1[j]->GetPointY(i) - levGraph_2[j]->GetPointY(i) ) /
    ( levGraph_1[j]->GetPointX(i) - levGraph_2[j]->GetPointX(i) );
}

//this is doing the sewing of gSF data
void ShapeGSF::DoInterpol() {
    if (m_sett->verbose)
        std::cout <<"\nINTERPOLATION OF gSF DATA... " <<endl;
    
    int j = levGraph_1.size() -1 ;
    int k = levGraph_1[j]->GetN();
    //if there is only one pair of non-zero gSF values, don't do anything
    if (k < 2) {
        std::cout <<"Nothing to interpolate in this iteration!" <<std::endl;
        return;
    }
    
    //interpolate from first bin to higher energies; there is nothing to do for the first bin
    for (int i = 1; i < k; i++) {
        // do average interpolation; this scaling makes the interpolation independent of the direction
        double scale_avg = ( 2 * levGraph_1[j]->GetPointY(i-1) - ( ( levGraph_1[j]->GetPointX(i-1) - levGraph_2[j]->GetPointX(i)) * Slope(i-1) ) ) / ( 2 * levGraph_2[j]->GetPointY(i) + ( (levGraph_1[j]->GetPointX(i-1) - levGraph_2[j]->GetPointX(i)) * Slope(i) ) );
       
        levGraph_1[j]->GetY()[i] *= scale_avg;
        levGraph_1[j]->GetEY()[i] *= scale_avg;
        levGraph_2[j]->GetY()[i] *= scale_avg;
        levGraph_2[j]->GetEY()[i] *= scale_avg;
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
    int j = levGraph_1.size()-1;
    levGraph_1[j]->Print();
    levGraph_2[j]->Print();
}

//transforms literature gSF values via B*exp(alpha E_gamma)
void ShapeGSF::Transform(double B_t, double alpha_t) {
    
    //gSF was previously transformed via B and Alpha, so only transform according to the change in B_t and Alpha_t
    if (m_sett->doMC) {
        for (int i = 0; i < mc_litGraph->GetN() ; i++ ) {
            mc_litGraph->GetY()[i] *=  B_t / m_B * TMath::Exp(( alpha_t - m_alpha) * mc_litGraph->GetX()[i] / 1000.);
        }
    }
    else {
        for (int i = 0; i < litGraph->GetN() ; i++ ) {
            litGraph->GetY()[i] *=  B_t / m_B * TMath::Exp(( alpha_t - m_alpha) * litGraph->GetX()[i] / 1000.);
            litGraph->GetEY()[i] *=  B_t / m_B * TMath::Exp(( alpha_t - m_alpha) * litGraph->GetX()[i] / 1000.);
        }
    }
    m_B = B_t;
    m_alpha = alpha_t;
}

void ShapeGSF::Scale(int j, double factor){

    for (int i = 0; i < levGraph_1[j]->GetN(); i++) {
        levGraph_1[j]->GetY()[i] *= factor;
        levGraph_2[j]->GetY()[i] *= factor;
        levGraph_1[j]->GetEY()[i] *= factor;
        levGraph_2[j]->GetEY()[i] *= factor;
    }
}

//scales all graphs with factor with respect to the raw data; the previous scaling was scale_bak

void ShapeGSF::ScaleAll(double factor){

    for (int j = 0; j < levGraph_1.size(); j++) {
        Scale(j, factor / scale_bak);
    }
    scale_bak = factor;
}

