#include "../Include/ShapeGSF.h"
#include <algorithm>

//constructor using setting file and a matrix
ShapeGSF::ShapeGSF(ShapeSetting* setting_t, ShapeMatrix* matrix_t)
{
    sett = setting_t;
    gSF_matrix = matrix_t;
    gSF.clear();
    binRange[0] = 1;
    binRange[1] = gSF_matrix->GetYBins();
}

//constructor using the literature input file specified in the settings file, containing gamma energies and gSF data
ShapeGSF::ShapeGSF(ShapeSetting* setting_t)
{
     sett = setting_t;
    //reset gSF_sort vector
    gSF_sort.clear();
    
    //read file data into gSF
    if (sett->osloFileName == "") {
        std::cout << "No Literature File loaded!"<<std::endl;
        return NULL;
    }
    
    //temporary vector to collect the file gSF data
    gSF_sor s;

    if (sett->verbose)
        std::cout <<"\nReading OSLO DATA... " <<endl;
    ifstream inp;
    inp.open(sett->osloFileName.c_str());
    if (inp.is_open() ) {
        
        double oslo_gSF_high;
        double oslo_gSF_low;
        
        while ( !inp.eof() ) {
            inp >> s.egamma >> oslo_gSF_high >>oslo_gSF_low;
            s.value = ( oslo_gSF_high + oslo_gSF_low ) /2;
            s.dvalue = oslo_gSF_high - s.value;
            s.peakID = 1;               //peak ID has no meaning here...so just setting this to 1
            gSF_sort.push_back(s);
            
            if (sett->verbose)
                std::cout<< s.egamma << " "<< s.value <<" "<<s.dvalue <<endl;
        }
    }
}

//returns an interpolated value of gSF for the gamma ray energy ene based on the gSF_sort vector
std::vector <double> ShapeGSF::InterpolValueSort(double ene) {
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
    
    int nOfActiveBins = gSF_sort.size();
    int index_l = -1;
    int index_r = -1;
    
    
    if (sett->verbose > 1)
        std::cout <<"Energy for Interpolation: " <<ene<<std::endl;
    
    //find clodest gamma ray lower than ene
    for (int i = 0; i < nOfActiveBins; i++) {
        t_ene =gSF_sort[i].egamma - ene;
        t_val =gSF_sort[i].value;
        t_dval =gSF_sort[i].dvalue;
        
        if (TMath::Abs(t_ene) < TMath::Abs(ene1) && (t_val > 0 ) && (t_ene <= 0)) {
            ene1 = t_ene;
            val1 = t_val;
            dval1 = t_dval;
            index_l = i;
            if (sett->verbose > 1)
                std::cout <<"Found closer energy smaller than ene: " <<ene1<< " with value: " <<val1 <<   std::endl;
        }
        
    }
    //find closest gamma ray higher than ene
    for (int j = 0; j < nOfActiveBins; j++) {
        t_ene =gSF_sort[j].egamma - ene;
        t_val =gSF_sort[j].value;
        t_dval =gSF_sort[j].dvalue;
        if (TMath::Abs(t_ene) < TMath::Abs(ene2) && (t_val > 0 ) && (t_ene > 0)) {
            ene2 = t_ene;
            val2 = t_val;
            dval2 = t_dval;
            index_r = j;
            if (sett->verbose > 1)
                std::cout <<"Found closer energy higher than ene: " <<ene2<< " with value: " <<val2 <<   std::endl;
        }
    }
    //check if ene1 and ene2 are defined
    if (index_l == -1) {
        ene1 = gSF_sort[0].egamma - ene;
        val1 = gSF_sort[0].value;
        dval1 = gSF_sort[0].dvalue;
        ene2 = gSF_sort[1].egamma - ene;
        val2 = gSF_sort[1].value;
        dval2 = gSF_sort[1].dvalue;
    }
    else if (index_r == -1) {
        ene1 = gSF_sort[nOfActiveBins-2].egamma - ene;
        val1 = gSF_sort[nOfActiveBins-2].value;
        dval1 = gSF_sort[nOfActiveBins-2].dvalue;
        ene2 = gSF_sort[nOfActiveBins-1].egamma - ene;
        val2 = gSF_sort[nOfActiveBins-1].value;
        dval2 = gSF_sort[nOfActiveBins-1].dvalue;
    }
    
    //transform to gamma ray energy
    ene1 = ene1 + ene;
    ene2 = ene2 + ene;
    if (sett->verbose > 1)
        std::cout <<"Final values ene1 and ene2 " <<ene1<< " and " <<ene2 <<   std::endl;
    //now interpolate gSF
    double m = ( val1 - val2) / ( ene1 - ene2);
    if (sett->verbose > 1)
        std::cout <<"slope calculated as : " << m <<   std::endl;
    std::vector <double> result;
    result.push_back ( ( ene - ene1) * m + val1);
    m = ( dval1 - dval2) / ( ene1 - ene2);
    result.push_back ( ( ene - ene1) * m + dval1);
    return result;
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
    
    if (sett->verbose > 1) {
        for (int i = 0; i < nOfActiveBins; i++) {
            std::cout <<"Current gSF energies gamma1: "<< gSF[binRange[0]+i-1].egamma1 <<std::endl;
             std::cout <<"Current gSF energies gamma2: "<< gSF[binRange[0]+i-1].egamma2 <<std::endl;
        }
    }
   
    if (sett->verbose > 1)
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
            if (sett->verbose > 1)
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
            if (sett->verbose > 1)
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
            if (sett->verbose > 1)
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
            if (sett->verbose > 1)
                std::cout <<"Found closer energy higher than ene: " <<ene2<< " with value: " <<val2 <<   std::endl;
        }
    }
    //check if ene1 and ene2 are defined
    if (index_l == -1) {
        ene1 = gSF[binRange[0]+index_r].egamma2 - ene;
        val1 = gSF[binRange[0]+index_r].value2;
        dval1 = gSF[binRange[0]+index_r].dvalue2;
        if (sett->verbose > 1)
            std::cout <<"ene1 wasn't defined, now setting to: " <<ene1<< " with value: " <<val1 <<   std::endl;
    }
    else if (index_r == -1) {
        ene2 = gSF[binRange[0]+index_l-2].egamma1 - ene;
        val2 = gSF[binRange[0]+index_l-2].value1;
        dval2 = gSF[binRange[0]+index_l-2].dvalue1;
        if (sett->verbose > 1)
            std::cout <<"ene2 wasn't defined, now setting to: " <<ene2<< " with value: " <<val2 <<   std::endl;

    }
    
    //transform to gamma ray energy
    ene1 = ene1 + ene;
    ene2 = ene2 + ene;
    if (sett->verbose > 1)
        std::cout <<"Final values ene1 and ene2 " <<ene1<< " and " <<ene2 <<   std::endl;
    //now interpolate gSF
    double m = ( val1 - val2) / ( ene1 - ene2);
    if (sett->verbose > 1)
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
    
    //interpolate from interEne to lower energies
      for (int i = gSF_matrix->energyToBinY(sett->sewingEne); i > (-1) ; i--) {
        //do nothing for the first loop
        if ( slope != 0 ) {
            
            //this is the interpolated value for gSF
            double y_norm = ( gSF[i].egamma1 -gSF[i+1].egamma1 ) * slope  + gSF[i+1].value1;
            if (gSF[i].value1 <= 0) {
                //in case the count rates are zero (or neagtive) the gSF value could still be defined through interpolation, but its error bar would not be defined anymore. One also cannot interpolate gSF[i].value2 anymore, so the entire bin, and all following bins, are ignored
                //erase all elements in gSF from current bin onwards
                binRange[0] = i+2; //will not include the current bin
                if (sett->verbose > 1)
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
                if (sett->verbose > 1)
                    std::cout <<"Setting lowest active bin to "<< i+1 <<std::endl;
                break;
            }
            
            gSF[i].value2 = gSF[i].value2 * y_norm / gSF[i].value1;
            gSF[i].dvalue2 = gSF[i].dvalue2 * y_norm / gSF[i].value1;
            
            gSF[i].dvalue1 = gSF[i].dvalue1 * y_norm / gSF[i].value1;
            gSF[i].value1 = y_norm;
            
            if (sett->verbose > 1)
                std::cout <<"\nValues of gSF after interpolation: "<< gSF[i].value1 <<" "<< gSF[i].value2 << endl;
        }
        
        //calculate slope for next iteration
      
        slope = ( gSF[i].value1-gSF[i].value2 ) / ( gSF[i].egamma1 - gSF[i].egamma2 );
        
        if (sett->verbose > 1)
            std::cout << "\nSlope for interpolation:" <<slope<<endl;
    }
    
    //interpolate from interEne to higher energies
    slope = 0;
      for (int i = gSF_matrix->energyToBinY(sett->sewingEne); i < gSF_matrix->GetYBins(); i++) {
        //do nothing for the first loop
        if ( slope != 0 ) {
            
            //this is the interpolated value for gSF
            double y_norm = ( gSF[i].egamma2 -gSF[i-1].egamma2 ) * slope  + gSF[i-1].value2;
            
            if (gSF[i].value2 <= 0) {
                //in case the count rates are zero (or neagtive) the gSF value could still be defined through interpolation, but its error bar would not be defined anymore. One also cannot interpolate gSF[i].value1 anymore, so the entire bin, and all following bins, are ignored
                binRange[1] = i; //will not include the current bin
                if (sett->verbose > 1)
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
                if (sett->verbose > 1)
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
            
            if (sett->verbose > 1)
                std::cout << "\nSlope for interpolation:" <<slope<<endl;
            if (sett->verbose > 1)
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
    
    else if (level == 2 ) {
        
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
    p_ratio.clear();
    gSF_matrix->Integrate();
    gSF_matrix->IntegrateBg();
    gSF_matrix->IntegrateSquare();
    gSF_matrix->IntegrateCube();
    
    //if autofit or background subtraction is set, perform autofit using Gauss
    if (sett->mode == 2 || sett->doBackground)
        gSF_matrix->FitIntegral();
    
    double elevel1 = ( sett->levEne[0] + sett->levEne[1] ) / 2;
    double elevel2 = ( sett->levEne[2] + sett->levEne[3] ) / 2;
    double egamma1, egamma2;
    
    for (int i = 0; i < gSF_matrix->integral1Cube.size(); i++ ) {
        
		if (getBgRatio(i, 1) * gSF_matrix->integral1[i] < sett->minCounts || gSF_matrix->integral1[i] == 0 ) {
            gSF_t.value1 = 0;
            gSF_t.dvalue1 = 0;
            //gSF_t.egamma1 = ( gSF_matrix->ybins[i] + gSF_matrix->ybins[i+1] ) / 2 - elevel1;
			gSF_t.egamma1 = gSF_matrix->GetEne0() + ( i + 0.5) * gSF_matrix->GetESize()  - elevel1;
        }
        else {
			
			//calculate E_gamma
            //this is the better way of calculating E_gamma by calculating the average E_g, weighted by gSF, i.e. egamma_avg = SUM (E_g*gSF) / SUM(gSF). E_g*gSF is contained in the gSFSquare matrix
            gSF_t.egamma1 = gSF_matrix->integral1Square[i] / gSF_matrix->integral1Cube[i];
            //this calculates the average E_gamma just using the middle of each excitation bin; not the preferred way as it doesn't take into account of feeding as a function iof excitation energy
            //gSF_t.egamma1 = ( gSF_matrix->ybins[i] + gSF_matrix->ybins[i+1] ) / 2 - elevel1;
            
	        //calculate gamma ray strength
			double M1 =  getBgRatio(i, 1) * gSF_matrix->integral1Cube[i] * sett->getEffCor(gSF_t.egamma1, 1);
            gSF_t.value1 = M1;
            gSF_t.dvalue1 = TMath::Power(1./gSF_matrix->integral1[i], 0.5)*gSF_t.value1;
            
            //recalculate gSF in case autofit is activated
            if (sett->mode == 2 && sett->doBackground) {
                gSF_t.value1 = sett->getEffCor(gSF_t.egamma1, 1) * gSF_matrix->fit_integral1Net[i] * gSF_matrix->integral1Cube[i]  / gSF_matrix->integral1[i];
                
                if (gSF_matrix->fit_integral1Net[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    gSF_t.value1 = 0;
                }
            }
            else if( sett->mode == 2 && !sett->doBackground) {
                gSF_t.value1 = sett->getEffCor(gSF_t.egamma1, 1) * gSF_matrix->fit_integral1[i] * gSF_matrix->integral1Cube[i]  / gSF_matrix->integral1[i];
                if (gSF_matrix->fit_integral1[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    gSF_t.value1 = 0;
                }
            }
        }
        
        
        if (getBgRatio(i, 2) * gSF_matrix->integral2[i] < sett->minCounts || gSF_matrix->integral2[i] == 0 ) {
            gSF_t.value2 = 0;
            gSF_t.dvalue2 = 0;
            gSF_t.egamma2 = gSF_t.egamma1 = gSF_matrix->GetEne0() + ( i + 0.5) * gSF_matrix->GetESize() - elevel2;
        }
        else {
			//calculate E_gamma
            //this is the better way of calculating E_gamma by calculating the average E_g, weighted by gSF, i.e. egamma_avg = SUM (E_g*gSF) / SUM(gSF). E_g*gSF is contained in the gSFSquare matrix
            gSF_t.egamma2 = gSF_matrix->integral2Square[i] / gSF_matrix->integral2Cube[i];
            //this calculates the average E_gamma just using the middle of each excitation bin; not the preferred way as it doesn't take into account of feeding as a function iof excitation energy
            //gSF_t.egamma2 = ( gSF_matrix->ybins[i] + gSF_matrix->ybins[i+1] ) / 2 - elevel2;
			
			//calculate gamma ray strength
			//gSF_t.value2 = getBgRatio(i, 2) * gSF_matrix->integral2Cube[i] * sett->eff_corr ;
			gSF_t.value2 = getBgRatio(i, 2) * gSF_matrix->integral2Cube[i] * sett->getEffCor(gSF_t.egamma2, 2);
            gSF_t.dvalue2 = TMath::Power(1./gSF_matrix->integral2[i], 0.5)*gSF_t.value2;
            
			//recalculate gSF in case autofit is activated
            if (sett->mode == 2 && sett->doBackground) {
                
				//gSF_t.value2 = sett->eff_corr * gSF_matrix->fit_integral2Net[i] * gSF_matrix->integral2Cube[i]  / gSF_matrix->integral2[i];
         		gSF_t.value2 = sett->getEffCor(gSF_t.egamma2, 2) * gSF_matrix->fit_integral2Net[i] * gSF_matrix->integral2Cube[i] / gSF_matrix->integral2[i];
                if (gSF_matrix->fit_integral2Net[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    gSF_t.value2 = 0;
                }
            }
            else if (sett->mode == 2 && !sett->doBackground) {
                gSF_t.value2 = sett->getEffCor(gSF_t.egamma2, 2) * gSF_matrix->fit_integral2[i] * gSF_matrix->integral2Cube[i]  / gSF_matrix->integral2[i];
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
            std::cout <<"gSF2: " <<gSF[i].value2 << "+- " << gSF[i].dvalue2 <<std::endl;
        }
        
        //calculate peak area ratios; this is useful for debugging and understanding detector efficiencies
        double rat;
        
        if (sett->mode == 1 && !sett->doBackground) {
            if (gSF_matrix->integral2[i] > 0)
                rat = gSF_matrix->integral1[i] / gSF_matrix->integral2[i];
            else
                rat = 0;
        }
        
        else if (sett->mode == 1 && sett->doBackground) {
            if (gSF_matrix->integral2[i] > 0)
                rat = getBgRatio(i, 1) * gSF_matrix->integral1[i] / gSF_matrix->integral2[i] / getBgRatio(i, 2);
            else
                rat = 0;
        }
        
        else if (sett->mode == 2 && sett->doBackground) {
            if (gSF_matrix->fit_integral2Net[i] > 0)
                rat = gSF_matrix->fit_integral1Net[i] / gSF_matrix->fit_integral2Net[i];
            else
                rat = 0;
        }
        
        else if (sett->mode == 2 && !sett->doBackground) {
            if (gSF_matrix->fit_integral2[i] > 0)
                rat = gSF_matrix->fit_integral1[i] / gSF_matrix->fit_integral2[i];
            else
                rat = 0;
        }
        //fill vector with peak ratio value
        p_ratio.push_back(rat);
        if (sett->verbose)
            std::cout <<"Ratio of Level 1 and Level 2 for bin " <<i <<" is: " <<rat<<std::endl;
        
    }
    
}

//returns a graph displaying the peak area ratios for level 1 and level 2 for each bin, taking into account the actual settings (interation autofit, background subtraction)

TGraph* ShapeGSF::getRatioGraph() {
    double x[p_ratio.size()];
    double y[p_ratio.size()];
   
    for (int i =0; i < p_ratio.size(); i++) {
        x[i] = gSF_matrix->GetEne0() + ( i + 0.5) * gSF_matrix->GetESize();
        y[i] = p_ratio[i];
    }
    TGraph *T = new TGraph(p_ratio.size(), x, y);
    return T;
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
        //adding a 10% systematic uncertainty due to fluctuations
        s.dvalue = gSF[binRange[0]+i-1].dvalue1 + ( 0.1 * s.value );
		s.peakID = 1;
        gSF_sort.push_back(s);
        
        s.egamma = gSF[binRange[0]+i-1].egamma2;
        s.value = gSF[binRange[0]+i-1].value2;
        //adding a 10% systematic uncertainty due to fluctuations
        s.dvalue = gSF[binRange[0]+i-1].dvalue2 + ( 0.1 * s.value );
		s.peakID = 2;
        gSF_sort.push_back(s);
    }
}



//prints the gSF values, sorted for egamma
void ShapeGSF::gSF_Print() {
    std::cout << "\n\nResults for gamma ray strength function: " <<std::endl;
    for (int i = 0; i < gSF_sort.size(); i++ )
        std::cout << gSF_sort[i].egamma <<"     "<< sett->gSF_norm * ( gSF_sort[i].value + gSF_sort[i].dvalue ) <<"     " << sett->gSF_norm *  ( gSF_sort[i].value - gSF_sort[i].dvalue ) << std::endl;
}

//transforms all gSF values via B*exp(alpha E_gamma)
void ShapeGSF::Transform(double B_t, double alpha_t) {
    
    //gSF was previously transformed via B and Alpha, so only transform according to the change in B_t and Alpha_t
    for (int i = 0; i < gSF_sort.size() ; i++ ) {
        gSF_sort[i].value =  B_t / B * TMath::Exp(( alpha_t - alpha) * gSF_sort[i].egamma / 1000.) *gSF_sort[i].value;
        gSF_sort[i].dvalue = B_t / B * TMath::Exp(( alpha_t - alpha) * gSF_sort[i].egamma / 1000.) *gSF_sort[i].dvalue;
    }
    B = B_t;
    alpha = alpha_t;
}

//provides a plot of gSF read from a file
TGraphErrors* ShapeGSF::plotLit() {
    
    int nOfPoints = gSF_sort.size();
    
    double x[nOfPoints], y[nOfPoints];
    double dx[nOfPoints], dy[nOfPoints];

    for (int i = 0; i < nOfPoints ; i++ ) {
        x[i] = gSF_sort[i].egamma;
        y[i] = gSF_sort[i].value;
        dx[i] = 1;
        dy[i] = gSF_sort[i].dvalue;
    }
    TGraphErrors *lit_data = new TGraphErrors(nOfPoints,x,y,dx,dy);
    lit_data->SetFillColor(4);
    lit_data->SetFillStyle(3010);
    return lit_data;
}

//sorting of gSF_sort vector

bool compare(gSF_sor &a, gSF_sor &b) { return a.egamma < b.egamma; }

void ShapeGSF::gSF_Sort() {
    std::sort(gSF_sort.begin(), gSF_sort.end(), compare);
}

//creats TGraphError using the sorted data and drawing option for plotitng a band
TMultiGraph* ShapeGSF::gSF_SortHisto(bool colour) {
    int nOfPoints = gSF_sort.size();
    
	//nOfPoints must always be an even number, or something is terribly wrong
	
	if (nOfPoints % 2 != 0) {
		std::cout <<"This is a bug! The number of gSF data points is not an even number! " <<std::endl;
		exit(0);
	}
	
	TMultiGraph *mg = new TMultiGraph();
	
	nOfPoints = nOfPoints / 2;
		
	double x1[nOfPoints], y1[nOfPoints];
	double dx1[nOfPoints], dy1[nOfPoints];
    	
	double x2[nOfPoints], y2[nOfPoints];
	double dx2[nOfPoints], dy2[nOfPoints];
    int id_1 = 0;
	int id_2 = 0;	
	for (int i = 0; i < 2*nOfPoints ; i++ ) {
        	
		if ( gSF_sort[i].peakID == 1 ) {
			x1[id_1] = gSF_sort[i].egamma;
			y1[id_1] = gSF_sort[i].value * sett->gSF_norm;
			dx1[id_1] = 1;
			dy1[id_1] = gSF_sort[i].dvalue * sett->gSF_norm;
			id_1++;
		}
		else {
			x2[id_2] = gSF_sort[i].egamma;
			y2[id_2] = gSF_sort[i].value * sett->gSF_norm;
			dx2[id_2] = 1;
			dy2[id_2] = gSF_sort[i].dvalue * sett->gSF_norm;
			id_2++;
		}
	}	
	TGraphErrors *gSFPlot_1 = new TGraphErrors(nOfPoints,x1,y1,dx1,dy1);
    gSFPlot_1->SetMarkerStyle(22);
    gSFPlot_1->SetMarkerSize(2);
    gSFPlot_1->SetMarkerColor(6);
	
	mg->Add(gSFPlot_1,"P");
	 
  	TGraphErrors *gSFPlot_2 = new TGraphErrors(nOfPoints,x2,y2,dx2,dy2);
    gSFPlot_2->SetMarkerStyle(22);
    gSFPlot_2->SetMarkerSize(2);
    if (colour) 
		gSFPlot_2->SetMarkerColor(7);
	else
		gSFPlot_2->SetMarkerColor(6);
    mg->Add(gSFPlot_2,"P");
	
	return mg;
}

//scales the sorted vector of gSF
void ShapeGSF::ScaleSort(double factor){
    
    for (int i = 0; i < gSF_sort.size(); i++) {
        gSF_sort[i].value = gSF_sort[i].value * factor;
        gSF_sort[i].dvalue = gSF_sort[i].dvalue * factor;
    }
}

void ShapeGSF::Scale(double factor){
    
    for (int i = 0; i < gSF_matrix->GetYBins(); i++) {
        gSF[i].value1 = gSF[i].value1 * factor;
        gSF[i].value2 = gSF[i].value2 * factor;
        gSF[i].dvalue1 = gSF[i].dvalue1 * factor;
        gSF[i].dvalue2 = gSF[i].dvalue2 * factor;
    }
}


//temp attempt to get level density from first generation matrix

void ShapeGSF::Rho(){
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
}

TGraphErrors* ShapeGSF::GetRhoDiagram() {
    return rhoDiagram;
}
