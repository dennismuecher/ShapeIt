#include "../Include/ShapeGSF_new.h"
//#include <algorithm>

//constructor using setting file and a matrix
ShapeGSF_new::ShapeGSF_new(ShapeSetting* t_sett, ShapeMatrix* t_matrix):m_sett(t_sett), m_matrix(t_matrix)
{
    //levGraph_1[0] = new TGraphErrors();
    //levGraph_2[0] = new TGraphErrors();
    multGraph = new TMultiGraph();
    multGraph->SetTitle("Gamma Ray Strength Function from Shape Method; E_{#gamma} (keV); f(E_{#gamma} (MeV^{-3})" );
    litGraph = new TGraphErrors();
}

//reads the literature values for gSF
void ShapeGSF_new::ReadLit()
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
        
        while ( !inp.eof() ) {
            inp >> e_gamma >> oslo_gSF_high >>oslo_gSF_low;
            litGraph->SetPoint(litGraph->GetN(), e_gamma,
                                    ( oslo_gSF_high + oslo_gSF_low ) / 2 );
            
            litGraph->SetPointError(litGraph->GetN() - 1, 0,
                                    ( oslo_gSF_high - oslo_gSF_low ) / 2 );
        }
        if (m_sett->verbose)
            litGraph->Print();
    }
}

//merges levGraph_1 and levGraph_2 into levGraph and sorts by energies
/*void ShapeGSF_new::Merge() {
    
    TObjArray *mergeGraph = new TObjArray(2);
    mergeGraph->Add(levGraph_1[j]);
    mergeGraph->Add(levGraph_2[j]);
    std::cout <<"Printing results for levGraph 1: " <<std::endl;
    levGraph_1[j]->Print();
    std::cout <<"Printing results for levGraph 2: " <<std::endl;
    levGraph_2[j]->Print();
    levGraph[j]->Merge(mergeGraph);
    std::cout <<"Printing results for levGraph merge: " <<std::endl;
    levGraph[j]->Print();
    
    //sort
    levGraph[j]->Sort();
}*/

//fills the actual vector of levGraph_1 and levGraph_2 with values; also fills levGraph
void ShapeGSF_new::FillgSF() {
   
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
    
    levGraph.push_back(new TGraphErrors());
    levGraph_1.push_back(new TGraphErrors());
    levGraph_2.push_back(new TGraphErrors());
    int j = levGraph_1.size() - 1;
    
    std::cout <<"levGraph_1.size(): " << levGraph_1.size() << std::endl;
    
    for (int i = 0; i < m_matrix->integral1Cube.size(); i++ ) {

        //if any of the two peak areas is below minCounts, skip this entire gSF pair
        if ( getBgRatio(i, 1) * m_matrix->integral1[i] < m_sett->minCounts ||
             getBgRatio(i, 2) * m_matrix->integral2[i] < m_sett->minCounts)
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
        
        //filling results into levGraph vector;
        //this adds a 10% systemnatic uncertainty; should not be hard coded!
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
    if (m_sett->doInterpol)
        DoInterpol();
}

//creates the graph of all gSF data; for this, first interpolate the individual runs, if requested by user; then add them to a multigraph
TMultiGraph* ShapeGSF_new::getMultGraph() {
    for (int i = 0; i < levGraph.size(); i++) {
       
        levGraph_1[i]->SetMarkerStyle(22);
        levGraph_1[i]->SetMarkerSize(2);
        levGraph_1[i]->SetMarkerColor(6);
        levGraph_2[i]->SetMarkerStyle(22);
        levGraph_2[i]->SetMarkerSize(2);
        if (m_sett->colour)
            levGraph_2[i]->SetMarkerColor(7);
        else
            levGraph_2[i]->SetMarkerColor(6);
         
        multGraph->Add(levGraph_1[i],"P");
        multGraph->Add(levGraph_2[i],"P");

    }
    //add literature values to graph
    if (m_sett->doOslo) {
        ReadLit();
        multGraph->Add(litGraph,"3A");
    }
    
    return ( multGraph );
}


double ShapeGSF_new::Slope(int i) {
    int j = levGraph_1.size()-1;
    return ( levGraph_1[j]->GetPointY(i) - levGraph_2[j]->GetPointY(i) ) /
    ( levGraph_1[j]->GetPointX(i) - levGraph_2[j]->GetPointX(i) );
}

void ShapeGSF_new::DoInterpol() {
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
        
        //this is the interpolated value for gSF
        double gSF2 = ( levGraph_2[j]->GetPointX(i) - levGraph_2[j]->GetPointX(i-1)) * Slope(i-1)  + levGraph_2[j]->GetPointY(i-1);
        double gSF1 = levGraph_1[j]->GetPointY(i) * gSF2 / levGraph_2[j]->GetPointY(i);
        double dgSF1 = levGraph_1[j]->GetEY()[i] * gSF2 / levGraph_2[j]->GetPointY(i);
        double dgSF2 = levGraph_2[j]->GetEY()[i] * gSF2 / levGraph_2[j]->GetPointY(i);
    
       if (m_sett->verbose > 1) {
            std::cout <<"Interpolation of bin " << i << std::endl;
            std::cout <<"\nValues of gSF after interpolation, peak 1: "<< levGraph_1[j]->GetPointX(i) <<" "<< gSF1 << endl;
            std::cout <<"\nValues of gSF after interpolation, peak 2: "<< levGraph_2[j]->GetPointX(i) <<" "<< gSF2 << endl;
        }
       
        
        // do average interpolation; this scaling makes the interpolation independent of the direction
        
        double scale_avg = ( 2 * levGraph_1[j]->GetPointY(i-1) - ( ( levGraph_1[j]->GetPointX(i-1) - levGraph_2[j]->GetPointX(i)) * Slope(i-1) ) ) / ( 2 * levGraph_2[j]->GetPointY(i) + ( (levGraph_1[j]->GetPointX(i-1) - levGraph_2[j]->GetPointX(i)) * Slope(i) ) );
        
        //levGraph_1[j]->SetPoint(i, levGraph_1[j]->GetPointX(i), levGraph_1[j]->GetPointY(i) * scale_avg );
        //levGraph_2[j]->SetPoint(i, levGraph_2[j]->GetPointX(i), levGraph_2[j]->GetPointY(i) * scale_avg );
        levGraph_1[j]->SetPointY(i, levGraph_1[j]->GetPointY(i) * scale_avg );
        levGraph_1[j]->SetPointError(i, 0, levGraph_1[j]->GetEY()[i] * scale_avg );
        levGraph_2[j]->SetPointY(i, levGraph_2[j]->GetPointY(i) * scale_avg );
        levGraph_2[j]->SetPointError(i, 0, levGraph_2[j]->GetEY()[i] * scale_avg );

        
        /*if (sett->verbose > 1) {
           std::cout <<"Calculating average interpolation " << std::endl;
          std::cout <<"\nValues of gSF after average interpolation, peak 1: "<< gSF[i].egamma1 <<" "<< gSF[i].value1 << endl;
          std::cout <<"\nValues of gSF after average interpolation, peak 2: "<< gSF[i].egamma2 <<" "<< gSF[i].value2 << endl;
        }*/
    }
}

//calculates the ratio of Peak to background for bin i

double ShapeGSF_new::getBgRatio(int bin, int level) {
    
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
void ShapeGSF_new::Print() {
    std::cout << "\n\nResults for gamma ray strength function: " <<std::endl;
    int j = levGraph_1.size()-1;
    levGraph_1[j]->Print();
    levGraph_2[j]->Print();
}

//transforms all gSF values via B*exp(alpha E_gamma)
void ShapeGSF_new::Transform(double B_t, double alpha_t) {
    
    //gSF was previously transformed via B and Alpha, so only transform according to the change in B_t and Alpha_t
   /* for (int i = 0; i < gSF_sort.size() ; i++ ) {
        gSF_sort[i].value =  B_t / B * TMath::Exp(( alpha_t - alpha) * gSF_sort[i].egamma / 1000.) *gSF_sort[i].value;
        gSF_sort[i].dvalue = B_t / B * TMath::Exp(( alpha_t - alpha) * gSF_sort[i].egamma / 1000.) *gSF_sort[i].dvalue;
    }
    B = B_t;
    alpha = alpha_t;*/
}

void ShapeGSF_new::Scale(double factor){
    int j = levGraph_1.size()-1;
    
    for (int i = 0; i < levGraph_1[j]->GetN(); i++) {
        levGraph_1[j]->GetY()[i] *= factor;
        levGraph_2[j]->GetY()[i] *= factor;
        levGraph_1[j]->GetEY()[i] *= factor;
        levGraph_2[j]->GetEY()[i] *= factor;
    }
}











