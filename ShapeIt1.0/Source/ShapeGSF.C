#include "../Include/ShapeGSF.h"
#include <algorithm>

ShapeGSF::ShapeGSF(ShapeSetting* setting_t, ShapeMatrix* matrix_t)
{
    sett = setting_t;
    gSF_matrix = matrix_t;
    gSF.clear();
    binRange[0] = 1;
    binRange[1] = sett->nOfBins;
}

//returns an interpolated value of gSF for the gamma ray energy ene
std::vector <double> ShapeGSF::InterpolValue(double ene) {
    //first, find the closest gamma ray energy for which there is data on gSF
    double ene1 = 1E5;
    double val1 = 1;
    double dval1 = 1;
    
    double ene2 = 1E5;
    double val2 = 1;
    double dval2 = 1;
    
    double t_ene;
    double t_val;
    double t_dval;
    
    int nOfActiveBins = (binRange[1] - binRange[0] +1);
    int index_l = -1;
    int index_r = -1;
    
    if (sett->verbose == 2) {
        for (int i = 0; i < nOfActiveBins; i++) {
            std::cout <<"Current gSF energies gamma1: "<< gSF[binRange[0]+i-1].egamma1 <<std::endl;
             std::cout <<"Current gSF energies gamma2: "<< gSF[binRange[0]+i-1].egamma2 <<std::endl;
        }
    }
   
    if (sett->verbose == 2)
        std::cout <<"Energy for Interpolation: " <<ene<<std::endl;
    //find clodest gamma ray lower than ene
    for (int i = 0; i < nOfActiveBins; i++) {
        t_ene =gSF[binRange[0]+i-1].egamma1 - ene;
        t_val =gSF[binRange[0]+i-1].value1;
        t_dval =gSF[binRange[0]+i-1].dvalue1;
        
        if (TMath::Abs(t_ene) < TMath::Abs(ene1) && (t_val > 0 ) && (t_ene <= 0)) {
            ene1 = t_ene;
            val1 = t_val;
            dval1 = t_dval;
            index_l = i;
            if (sett->verbose == 2)
                std::cout <<"Found closer energy smaller than ene: " <<ene1<< " with value: " <<val1 <<   std::endl;
        }
        t_ene =gSF[binRange[0]+i-1].egamma2 - ene;
        t_val =gSF[binRange[0]+i-1].value2;
        t_dval =gSF[binRange[0]+i-1].dvalue2;

        if (TMath::Abs(t_ene) < TMath::Abs(ene1) && (t_val > 0 ) && (t_ene <= 0) ) {
            ene1 = t_ene;
            val1 = t_val;
            dval1 =t_dval;
            index_l = i;
            if (sett->verbose == 2)
                std::cout <<"Found closer energy smaller than ene: " <<ene1<< " with value: " <<val1 <<   std::endl;
        }
    }
    //find closest gamma ray higher than ene
    for (int j = 0; j < nOfActiveBins; j++) {
        t_ene =gSF[binRange[0]+j-1].egamma1 - ene;
        t_val =gSF[binRange[0]+j-1].value1;
        t_dval =gSF[binRange[0]+j-1].dvalue1;
        if (TMath::Abs(t_ene) < TMath::Abs(ene2) && (t_val > 0 ) && (t_ene > 0)) {
            ene2 = t_ene;
            val2 = t_val;
            dval2 = t_dval;
            index_r = j;
            if (sett->verbose == 2)
                std::cout <<"Found closer energy higher than ene: " <<ene2<< " with value: " <<val2 <<   std::endl;
        }
        t_ene =gSF[binRange[0]+j-1].egamma2 - ene;
        t_val =gSF[binRange[0]+j-1].value2;
        t_dval =gSF[binRange[0]+j-1].dvalue2;
        if (TMath::Abs(t_ene) < TMath::Abs(ene2) && (t_val > 0 ) && (t_ene > 0) ) {
            ene2 = t_ene;
            val2 = t_val;
            dval2 = t_dval;
            index_r = j;
            if (sett->verbose == 2)
                std::cout <<"Found closer energy higher than ene: " <<ene2<< " with value: " <<val2 <<   std::endl;
        }
    }
    //check if ene1 and ene2 are defined
    if (index_l == -1) {
        ene1 = gSF[binRange[0]+index_r].egamma2 - ene;
        val1 = gSF[binRange[0]+index_r].value2;
        dval1 = gSF[binRange[0]+index_r].dvalue2;
        if (sett->verbose == 2)
            std::cout <<"ene1 wasn't defined, now setting to: " <<ene1<< " with value: " <<val1 <<   std::endl;
    }
    else if (index_r == -1) {
        ene2 = gSF[binRange[0]+index_l-2].egamma1 - ene;
        val2 = gSF[binRange[0]+index_l-2].value1;
        dval2 = gSF[binRange[0]+index_l-2].dvalue1;
        if (sett->verbose == 2)
            std::cout <<"ene2 wasn't defined, now setting to: " <<ene2<< " with value: " <<val2 <<   std::endl;

    }
    
    //transform to gamma ray energy
    ene1 = ene1 + ene;
    ene2 = ene2 + ene;
    if (sett->verbose == 2)
        std::cout <<"Final values ene1 and ene2 " <<ene1<< " and " <<ene2 <<   std::endl;
    //now interpolate gSF
    double m = ( val1 - val2) / ( ene1 - ene2);
    if (sett->verbose == 2)
        std::cout <<"slope calculated as : " << m <<   std::endl;
    std::vector <double> result;
    result.push_back ( ( ene - ene1) * m + val1);
    m = ( dval1 - dval2) / ( ene1 - ene2);
    result.push_back ( ( ene - ene1) * m + dval1);
    return result;
}

void ShapeGSF::DoInterpol() {
    if (sett->verbose)
        std::cout <<"\nINTERPOLATION OF gSF DATA... " <<endl;
    
    //slope of the last two points
    double slope = 0;
    
    //interpolate from interPoint to lower energies
    for (int i = (sett->interPoint -1 ); i > (-1) ; i--) {
        
        //do nothing for the first loop
        if ( slope != 0 ) {
            
            //this is the interpolated value for gSF
            double y_norm = ( gSF[i].egamma1 -gSF[i+1].egamma1 ) * slope  + gSF[i+1].value1;
            
            if (gSF[i].value1 <= 0) {
                //in case the count rates are zero (or neagtive) the gSF value could still be defined through interpolation, but its error bar would not be defined anymore. One also cannot interpolate gSF[i].value2 anymore, so the entire bin, and all following bins, are ignored
                //erase all elements in gSF from current bin onwards
                binRange[0] = i+2; //will not include the current bin
                if (sett->verbose == 2)
                    std::cout <<"Setting lowest active bin to "<< i+2 <<std::endl;
                break;
            }
            else if (gSF[i].value2 <= 0) {
                //in this case, gSF[i].value1 could still be useful, but interpolation beyond this bin will not be useful, so all earlier bins will be ignored
                gSF[i].value2 = 0; gSF[i].egamma2 = 0;
                gSF[i].dvalue1 = gSF[i].dvalue1 * y_norm / gSF[i].value1;
                gSF[i].value1 = y_norm;
                //binRange[0] = i+1;//will include the current bin and terminate the interpolation
                binRange[0] = i+2;//will NOT include the current bin and terminate the interpolation
                if (sett->verbose == 2)
                    std::cout <<"Setting lowest active bin to "<< i+1 <<std::endl;
                break;
            }
            
            gSF[i].value2 = gSF[i].value2 * y_norm / gSF[i].value1;
            gSF[i].dvalue2 = gSF[i].dvalue2 * y_norm / gSF[i].value1;
            
            gSF[i].dvalue1 = gSF[i].dvalue1 * y_norm / gSF[i].value1;
            gSF[i].value1 = y_norm;
            
            if (sett->verbose == 2)
                std::cout <<"\nValues of gSF after interpolation: "<< gSF[i].value1 <<" "<< gSF[i].value2 << endl;
        }
        
        //calculate slope for next iteration
      
        slope = ( gSF[i].value1-gSF[i].value2 ) / ( gSF[i].egamma1 - gSF[i].egamma2 );
        
        if (sett->verbose == 2)
            std::cout << "\nSlope for interpolation:" <<slope<<endl;
    }
    
    //interpolate from interPoint to higher energies
    slope = 0;
    for (int i = (sett->interPoint -1 ); i < sett->nOfBins ; i++) {
       
        //do nothing for the first loop
        if ( slope != 0 ) {
            
            //this is the interpolated value for gSF
            double y_norm = ( gSF[i].egamma2 -gSF[i-1].egamma2 ) * slope  + gSF[i-1].value2;
            
            if (gSF[i].value2 <= 0) {
                //in case the count rates are zero (or neagtive) the gSF value could still be defined through interpolation, but its error bar would not be defined anymore. One also cannot interpolate gSF[i].value1 anymore, so the entire bin, and all following bins, are ignored
                binRange[1] = i; //will not include the current bin
                if (sett->verbose == 2)
                    std::cout <<"Setting number of active bins to "<< i <<std::endl;
                break;
            }
            else if (gSF[i].value1 <= 0) {
                //in this case, gSF[i].value2 could still be useful, but interpolation beyond this bin will not be useful, so all following bins will be ignored
                gSF[i].value1 = 0; gSF[i].egamma1 = 0;
                gSF[i].dvalue2 = gSF[i].dvalue2 * y_norm / gSF[i].value2;
                gSF[i].value2 = y_norm;
                //binRange[1] = i+1; //will include the current bin
                binRange[1] = i; //will NOT include the current bin
                if (sett->verbose == 2)
                    std::cout <<"Setting number of active bins to "<< i+1 <<std::endl;
                break;
            }
            else {
                //both gSF values are > 0
                gSF[i].value1 = gSF[i].value1 * y_norm / gSF[i].value2;
                gSF[i].dvalue1 = gSF[i].dvalue1 * y_norm / gSF[i].value2;
                gSF[i].dvalue2 = gSF[i].dvalue2 * y_norm / gSF[i].value2;
                gSF[i].value2 = y_norm;
                
            }
            
            if (sett->verbose == 2)
                std::cout << "\nSlope for interpolation:" <<slope<<endl;
            if (sett->verbose == 2)
                std::cout <<"\nValues of gSF after interpolation: "<< gSF[i].value1 <<" "<< gSF[i].value2 << endl;
        }
        //calculate slope for next iteration
        slope = ( gSF[i].value1-gSF[i].value2 ) / ( gSF[i].egamma1 - gSF[i].egamma2 );

    }
}

//calculates the ratio of Peak to background for bin i

double ShapeGSF::getBgRatio(int bin, int level) {
    
    if (!sett->doBackground)
        return 1;
    
    if (level == 1 ) {
        
        double peak = gSF_matrix->integral1[bin];
        double backgr = gSF_matrix->fit_integral1Bg[bin];
        if (peak > 0)
            return (peak - backgr) / peak;
        else
            return 0;
    }
    
    if (level == 2 ) {
        
        double peak = gSF_matrix->integral2[bin];
        double backgr =gSF_matrix->fit_integral2Bg[bin];
        if (peak > 0)
            return (peak - backgr) / peak;
        else
            return 0;
    }
    return 0;
}


void ShapeGSF::FillgSF() {
    gSF_str gSF_t;
    gSF.clear();
    
    gSF_matrix->Integrate();
    gSF_matrix->IntegrateBg();
    gSF_matrix->IntegrateSquare();
    gSF_matrix->IntegrateCube();
    
    //if autofit or background usbtraction is set, perform autofit using Gauss
    if (sett->mode == 2 || sett->doBackground)
        gSF_matrix->FitIntegral();
    
    double elevel1 = ( sett->levEne[0] + sett->levEne[1] ) / 2;
    double elevel2 = ( sett->levEne[2] + sett->levEne[3] ) / 2;
    double egamma1, egamma2;
    
    for (int i = 0; i < gSF_matrix->integral1Cube.size(); i++ ) {
        
        //calculate gamma ray strength
        gSF_t.value1 = getBgRatio(i, 1) * gSF_matrix->integral1Cube[i] * sett->gSF_norm;
        gSF_t.value2 = getBgRatio(i, 2) * gSF_matrix->integral2Cube[i] * sett->eff_corr * sett->gSF_norm;
        
        if (getBgRatio(i, 1) * gSF_matrix->integral1[i] < sett->minCounts || gSF_matrix->integral1[i] == 0 ) {
            gSF_t.value1 = 0;
            gSF_t.dvalue1 = 0;
            gSF_t.egamma1 = ( gSF_matrix->ybins[i] + gSF_matrix->ybins[i+1] ) / 2 - elevel1;
        }
        else {
            gSF_t.dvalue1 = TMath::Power(1./gSF_matrix->integral1[i], 0.5)*gSF_t.value1;
            //this is the better way of calculating E_gamma by calculating the average E_g, weighted by gSF, i.e. egamma_avg = SUM (E_g*gSF) / SUM(gSF). E_g*gSF is contained in the gSFSquare matrix
            gSF_t.egamma1 = gSF_matrix->integral1Square[i] / gSF_matrix->integral1Cube[i];
            //this calculates the average E_gamma just using the middle of each excitation bin; not the preferred way as it doesn't take into account of feeding as a function iof excitation energy
            //gSF_t.egamma1 = ( gSF_matrix->ybins[i] + gSF_matrix->ybins[i+1] ) / 2 - elevel1;
            
            //recalculate gSF in case autofit is activated
            if (sett->mode == 2 && sett->doBackground) {
                gSF_t.value1 = gSF_matrix->fit_integral1Net[i] * gSF_matrix->integral1Cube[i] * sett->gSF_norm / gSF_matrix->integral1[i];
                if (gSF_matrix->fit_integral1Net[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    gSF_t.value1 = 0;
                }
            }
            else if( sett->mode == 2 && !sett->doBackground) {
                gSF_t.value1 = gSF_matrix->fit_integral1[i] * gSF_matrix->integral1Cube[i] * sett->gSF_norm / gSF_matrix->integral1[i];
                if (gSF_matrix->fit_integral1[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    gSF_t.value1 = 0;
                }
            }
        }
        
        
        if (getBgRatio(i, 2) * gSF_matrix->integral2[i] < sett->minCounts || gSF_matrix->integral2[i] == 0 ) {
            gSF_t.value2 = 0;
            gSF_t.dvalue2 = 0;
            gSF_t.egamma2 = ( gSF_matrix->ybins[i] + gSF_matrix->ybins[i+1] ) / 2 - elevel2;
            
        }
        else {
            gSF_t.dvalue2 = TMath::Power(1./gSF_matrix->integral2[i], 0.5)*gSF_t.value2;
            //this is the better way of calculating E_gamma by calculating the average E_g, weighted by gSF, i.e. egamma_avg = SUM (E_g*gSF) / SUM(gSF). E_g*gSF is contained in the gSFSquare matrix
            gSF_t.egamma2 = gSF_matrix->integral2Square[i] / gSF_matrix->integral2Cube[i];
            //this calculates the average E_gamma just using the middle of each excitation bin; not the preferred way as it doesn't take into account of feeding as a function iof excitation energy
            //gSF_t.egamma2 = ( gSF_matrix->ybins[i] + gSF_matrix->ybins[i+1] ) / 2 - elevel2;
            
            if (sett->mode == 2 && sett->doBackground) {
                gSF_t.value2 = sett->eff_corr * gSF_matrix->fit_integral2Net[i] * gSF_matrix->integral2Cube[i] * sett->gSF_norm / gSF_matrix->integral2[i];
                if (gSF_matrix->fit_integral2Net[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    gSF_t.value2 = 0;
                }
            }
            else if (sett->mode == 2 && !sett->doBackground) {
                gSF_t.value2 = sett->eff_corr * gSF_matrix->fit_integral2[i] * gSF_matrix->integral2Cube[i] * sett->gSF_norm / gSF_matrix->integral2[i];
                if (gSF_matrix->fit_integral2[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    gSF_t.value2 = 0;
                }
            }
        }
        
        gSF.push_back ( gSF_t);
        if (sett->verbose) {
            std::cout <<"Bin: " <<i+1<<std::endl;
            std::cout <<"energies: " <<gSF[i].egamma1 << " " << gSF[i].egamma2 <<std::endl;
            std::cout <<"gSF1: " <<gSF[i].value1 << "+- " << gSF[i].dvalue1 <<std::endl;
            std::cout <<"gSF2: " <<gSF[i].value2 << "+- " << gSF[i].dvalue2 <<std::endl;   }
    }
    
}

void ShapeGSF::Update() {
    //calcualte gSF
    FillgSF();

    if (sett->doInterpol)
        DoInterpol();
    if (sett->verbose)
        std::cout << "Active bins: " <<binRange[0] << " - " <<binRange[1] <<std::endl;
}


//collects the values of gSF and egamma in a vector
void ShapeGSF::gSF_Collect() {
    
    int nOfActiveBins = (binRange[1] - binRange[0] +1);
    gSF_sor s;
    for (int i = 0; i < nOfActiveBins ; i++ ) {
        
        s.egamma = gSF[binRange[0]+i-1].egamma1;
        s.value = gSF[binRange[0]+i-1].value1;
        s.dvalue = gSF[binRange[0]+i-1].dvalue1;
        gSF_sort.push_back(s);
        
        s.egamma = gSF[binRange[0]+i-1].egamma2;
        s.value = gSF[binRange[0]+i-1].value2;
        s.dvalue = gSF[binRange[0]+i-1].dvalue2;
        gSF_sort.push_back(s);
    }
}

bool compare(gSF_sor &a, gSF_sor &b) { return a.egamma < b.egamma; }


//prints the gSF values, sorted for egamma
void ShapeGSF::gSF_Print() {
    std::cout << "\n\nResults for gamma ray strength function: " <<std::endl;
    std::sort(gSF_sort.begin(), gSF_sort.end(), compare);
    for (int i = 0; i < gSF_sort.size(); i++ )
        std::cout << gSF_sort[i].egamma <<"     "<< gSF_sort[i].value<<"     " << gSF_sort[i].dvalue <<std::endl;
}

//creats TGraphError using the sorted data and drawing option for plotitng a band
TGraphErrors* ShapeGSF::gSF_SortHisto() {
    std::sort(gSF_sort.begin(), gSF_sort.end(), compare);
    int nOfPoints = gSF_sort.size();
    double x[nOfPoints], y[nOfPoints];
    double dx[nOfPoints], dy[nOfPoints];
    for (int i = 0; i < nOfPoints ; i++ ) {
        x[i] = gSF_sort[i].egamma;
        y[i] = gSF_sort[i].value;
        dx[i] = 1;
        dy[i] = gSF_sort[i].dvalue;
    }
    TGraphErrors *sortPlot = new TGraphErrors(nOfPoints,x,y,dx,dy);
    return sortPlot;
}



TGraphErrors* ShapeGSF::gSF_Histo() {
    
    int nOfActiveBins = (binRange[1] - binRange[0] +1);
    double x[2*nOfActiveBins], y[2*nOfActiveBins];
    double dx[2*nOfActiveBins], dy[2*nOfActiveBins];
    for (int i = 0; i < nOfActiveBins ; i++ ) {
        x[2*i] = gSF[binRange[0]+i-1].egamma1; dx[2*i] = 1;
        y[2*i] = gSF[binRange[0]+i-1].value1;  dy[2*i] = gSF[binRange[0]+i-1].dvalue1;
        
        x[2*i+1] = gSF[binRange[0]+i-1].egamma2; dx[2*i+1] = 1;
        y[2*i+1] = gSF[binRange[0]+i-1].value2;  dy[2*i+1] = gSF[binRange[0]+i-1].dvalue2;
    }
    TGraphErrors *gSFPlot = new TGraphErrors(2*nOfActiveBins,x,y,dx,dy);
    return gSFPlot;
}

void ShapeGSF::Scale(double factor){
    
    for (int i = 0; i < sett->nOfBins; i++) {
        gSF[i].value1 = gSF[i].value1 * factor;
        gSF[i].value2 = gSF[i].value2 * factor;
        gSF[i].dvalue1 = gSF[i].dvalue1 * factor;
        gSF[i].dvalue2 = gSF[i].dvalue2 * factor;
    }
}
TGraphErrors* ShapeGSF::plotLit() {
    if (sett->osloFileName == "") {
        std::cout << "No Literature File loaded!"<<std::endl;
        return NULL;
    }
        
    if (sett->verbose)
        std::cout <<"\nPLOTTING OSLO DATA... " <<endl;
    ifstream inp;
    inp.open(sett->osloFileName.c_str());
    if (inp.is_open() ) {
        double oslo_gamma[ 100 ];
        double oslo_gSF_high[ 100 ];
        double oslo_gSF_low[ 100 ];
        double oslo_gSF[ 100 ];
        
        double oslo_dgamma[ 100 ];
        double oslo_dgSF[ 100 ];
        int i = 0;
        while ( !inp.eof() ) {
            inp >> oslo_gamma[i] >> oslo_gSF_high[i] >>oslo_gSF_low[i];
            if (sett->verbose)
                std::cout<<oslo_gamma[i] << " "<< oslo_gSF_high[i] <<" "<<oslo_gSF_low[i] <<endl;
            oslo_gSF[i] = ( oslo_gSF_high[i] + oslo_gSF_low[i] ) /2;
            oslo_dgSF[i] = TMath::Abs( oslo_gSF_high[i] - oslo_gSF[i] );
            oslo_dgamma[i] = 1; //not sure what this should be
            //if (oslo_gamma[i] > ( gSF_matrix->eMax_y ) ) break;
            i = i + 1;
        }
        TGraphErrors *oslo_data = new TGraphErrors(i-1, oslo_gamma, oslo_gSF, oslo_dgamma, oslo_dgSF);
        
        oslo_data->SetFillColor(4);
        oslo_data->SetFillStyle(3010);
        return oslo_data;
    }
    else
        throw "Cannot open Oslo Data File! Provide correct file or set 'plot_Oslo' to 0 in input file! Exciting...";
    return NULL;
}
