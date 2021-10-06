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

//constructor using setting file and a matrix
ShapeCollector::ShapeCollector(ShapeSetting* t_sett, ShapeMatrix* t_matrix):m_sett(t_sett), m_matrix(t_matrix)
{
    //get literature gSF data
    litCollector = new ShapeGSF(m_sett);
    gSFGraph = new TGraphErrors();
    gSFGraphSmooth = new TGraphAsymmErrors();

}

// runs all different iterations of gSF following the user inout stored in the settings file m_sett
void ShapeCollector::Collect() {
    
    m_matrix->SetEne0( m_sett->exiEne[0] );
    m_matrix->SetEne1( m_sett->exiEne[1] );
    m_matrix->SetESize( m_sett->exi_size[0] );
    
    //reset gSF collector vector
    gSFCollector.clear();
    
    do {
        //sliding window
        for (int i = 1; i < kmax; i++) {
            
            // the sliding window moves from the inital position in kmax steps to the end position, which is one bin to the "left"
            // the high energy is kept at the initital value, at all times
            m_matrix->SetEne0( m_sett->exiEne[0] - ( (double) (i-1) * m_matrix->GetESize() / (kmax-1)));
            
           //get gSF data
            gSFCollector.push_back( new ShapeGSF(m_sett, m_matrix));
            
            //if no sliding window is requested, stop here
            if (!m_sett->doSlidingWindow)
                break;
        }
        
        //check if bin size should be varried
        if (!m_sett->doBinVariation) break;
        
        //change bin size
        m_matrix->SetESize( m_matrix->GetESize() + 50 );
        //update matrix and recalculate gSF in case doSlidingWindow is not active; otherwise this is done in the sliding window loop
        
    } while (m_matrix->GetESize() <= m_sett->exi_size[1]);
    
    //normalize all data to each other
    NormCollect();
}

void ShapeCollector::Draw() {
    
    getMultGraph()->Draw("AP");
}

void ShapeCollector::Print() {
    for (auto coll : gSFCollector)
        coll->Print();
    if (m_sett->doOslo)
        litCollector->Print();
}

//normalizes all gSF iterations to each other
void ShapeCollector::NormCollect() {
    
    double norm_fac = 1;
    
    for (int i = 0; i < gSFCollector.size(); i++) {
        if (m_sett->doOslo && m_sett->doAutoScale)
            norm_fac = Norm(gSFCollector[i], litCollector);
        else
            norm_fac = Norm(gSFCollector[i], gSFCollector[0]);
    
        gSFCollector[i]->Scale(norm_fac);
    }
    Merge();
}
    
TMultiGraph* ShapeCollector::getMultGraph() {
    
    //define Multigraph
    TMultiGraph* multGraph = new TMultiGraph();
    multGraph->SetTitle("Gamma Ray Strength Function from Shape Method; E_{#gamma} (keV); f(E_{#gamma} (MeV^{-3})" );
    
    //add literature values
    if (m_sett->doOslo)
        multGraph->Add(litCollector->GetLevGraph(),"3");
    
    //add all the gSF data stored in the collector to the MultiGraph
    for (int i = 0; i < gSFCollector.size(); i++) {
        multGraph->Add(gSFCollector[i]->GetLevGraph_1());
        multGraph->Add(gSFCollector[i]->GetLevGraph_2());
    }
    
    //add smoothed graph, if requested
    if (m_sett->displayAvg) {
        Smooth(0);
        multGraph->Add(gSFGraphSmooth,"P");
    }
    return ( multGraph );
}

//merges all gSF data in the gSFCollector vector into a single graph and sorts via gamma energy
void ShapeCollector::Merge() {
    
    TObjArray *mArray = new TObjArray();
    gSFGraph->Set(0);
    
    for (int i = 0; i < gSFCollector.size(); i++)             mArray->Add(gSFCollector[i]->GetLevGraph());
    
    gSFGraph->Merge(mArray);
    gSFGraph->Sort();
}

//calcualtion of a normalization factor to match T1 to T2
double ShapeCollector::Norm(ShapeGSF* T1, ShapeGSF* T2 ) {
    
    double c1 = 0, c2 = 0;
    std::vector<double> a;
    std::vector<double> da;
    std::vector<double> b;
    
    for (int i = 0; i < T1->GetLevGraph()->GetN(); i++) {
        a.push_back( T1->GetLevGraph()->GetY()[i] );
        da.push_back( T1->GetLevGraph()->GetEY()[i] );
        b.push_back( T2->GetLevGraph()->Eval( T1->GetLevGraph()->GetX()[i] ) );
    }
    
    for (int i = 0; i < T1->GetLevGraph()->GetN(); i++) {
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
            std::cout <<"Error in Norm(): zero error value detected" <<std::endl;
            return (0);
        }
    }
    if (c2 !=0)
        return (c1 / c2);
    else {
        std::cout <<"Error in Norm(): zero c2 value detected" <<std::endl;
        return (0);
    }
}

void ShapeCollector::Smooth(int res) {
    
    gSFGraphSmooth->Set(0);
    // if res == 0 use exi_size[0] (user input) for the bin size
    if (res ==0)
        res = m_sett->exi_size[0];
    auto nPoints = gSFGraph->GetN();
    double m_e = gSFGraph->GetX()[0] + res;
    std::vector <int> i_bin;          //the point number limits of levGraphAll belonging to the bins with size res
    i_bin.push_back(0);             //first bin always starts at zero
    for(int i =1; i < nPoints; i++) {
        if (gSFGraph->GetX()[i] > m_e) {
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
            x += gSFGraph->GetX()[j];
            y += gSFGraph->GetY()[j];
            dy += TMath::Power(gSFGraph->GetEY()[j] , 2);
            
        }
        int nOfP = i_bin[i+1] - i_bin[i];     //number of points in this bin
        if (nOfP ==0) continue;
        x = x / nOfP;
        y = y / nOfP;
        //this error is the error of the mean value, based on the propagation of the indiviual errors, only; this does not take into account the fluctuations in the data points (i.e. the standard deviation)
        dy = TMath::Sqrt(dy) / nOfP;
        gSFGraphSmooth->SetPoint(gSFGraphSmooth->GetN(),x,y);
        
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
                double m_dyh = gSFGraph->GetY()[j] + gSFGraph->GetEY()[j] - y;
                if (m_dyh > dyh)
                    dyh = m_dyh;
                //difference of average to lower error bar
                double m_dyl = y - gSFGraph->GetY()[j] + gSFGraph->GetEY()[j];
                if (m_dyl > dyl)
                dyl = m_dyl;
        }
        gSFGraphSmooth->SetPointError(gSFGraphSmooth->GetN()-1, 0, 0, dyl, dyh);
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
        
    gSFGraphSmooth->SetMarkerStyle(22);
    gSFGraphSmooth->SetMarkerSize(2);
    gSFGraphSmooth->SetMarkerColor(1);
}
