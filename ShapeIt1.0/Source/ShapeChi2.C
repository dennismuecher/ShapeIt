

#include "../Include/ShapeChi2.h"


//this constructor copies the gSF values during initialization and later compres the ref values to the up-to-date values of the gSF_t pointer
ShapeChi2::ShapeChi2(ShapeGSF *gSF_t, ShapeSetting *sett_t) {
    gSF = gSF_t;
    sett = sett_t;
    //fill vector with initital energies and gSF values
    ref_ene.clear();
    ref_val.clear();
    for (int i = gSF->binRange[0]-1; i < gSF->binRange[1]; i++) {
        ref_ene.push_back(gSF->gSF[i].egamma1);
        ref_val.push_back(gSF->gSF[i].value1);
        ref_ene.push_back(gSF->gSF[i].egamma2);
        ref_val.push_back(gSF->gSF[i].value2);
    }
    mode = 0;
}

//this constructor compares values in gSF_ref to those in gSF_t
ShapeChi2::ShapeChi2(ShapeGSF *gSF_t, ShapeGSF *gSF_ref_t, ShapeSetting *sett_t) {
    gSF = gSF_t;
    gSF_ref = gSF_ref_t;
    sett = sett_t;
    //fill vector with initital energies and gSF values
    ref_ene.clear();
    ref_val.clear();
    ref_dval.clear();
    for (int i = 0; i < gSF_ref->gSF_sort.size(); i++) {
        ref_ene.push_back(gSF_ref->gSF_sort[i].egamma);
        ref_val.push_back(gSF_ref->gSF_sort[i].value);
        ref_dval.push_back(gSF_ref->gSF_sort[i].dvalue);
        
    }
    mode = 1;
}

//prints a list of energies and scaling factors which aligns gSF with the reference values
void ShapeChi2::printScalingSort() {
    double e;
    double factor;
    //std::cout <<"Scaling factors to match data gSF to literature gSF:"<<std::endl;
    std::vector <double> d;
     for (int i = 0; i < gSF->gSF_sort.size() ;  i++) {
         e = gSF->gSF_sort[i].egamma;
         d = gSF_ref->InterpolValueSort(e);
         factor = d[0] / ( sett->gSF_norm * gSF->gSF_sort[i].value );
        // std::cout <<e <<" " <<factor<<std::endl;
        // std::cout << "gSF value from data: " << sett->gSF_norm * gSF->gSF_sort[i].value<<std::endl;
         //std::cout << "Interpolated value from Oslo data: " << d[0] <<std::endl;
         //std::cout << "Peak ID: " << gSF->gSF_sort[i].peakID <<std::endl;
     }
}

//chi2 minimization for sorted vector of gSF
double ShapeChi2::minChi2Sort(double scale) {
    double chi2m = 0;
    std::vector <double> d;
    for (int i = 0; i < gSF->gSF_sort.size()-1 ;  i++) {
        d = gSF_ref->InterpolValueSort(gSF->gSF_sort[i].egamma);
        
        //calculate chi2
        if (gSF->gSF_sort[i].dvalue != 0 )
            chi2m += TMath::Power(scale * d[0] - (gSF->gSF_sort[i].value ),2) / TMath::Power(gSF->gSF_sort[i].dvalue,2);
        
    }
    //return ( chi2m / gSF->gSF_sort.size() );
    return chi2m;

}

double ShapeChi2::minChi2(double scale) {
    double chi2m = 0;
    std::vector <double> d;
    
	//this should be fixed....a proper way of defining what to take into account during chi2 min! 
			
	for (int i = 0; i < ref_ene.size() - 1; i++) {
        //std::cout << "in loop: " <<i <<std::endl;
        d = gSF->InterpolValue(ref_ene[i]);
        
        if (d[0] < 0 ) {
            std::cout <<"negative gSF value found in Chi2 routine "<<std::endl;
            continue;
        }
        if (d[0] == 0 ) {
            std::cout <<"zero gSF value found in Chi2 routine "<<std::endl;
            continue;
        }
        if (d[1] == 0 ) {
            std::cout <<"zero gSF error value found in Chi2 routine "<<std::endl;
            continue;
        }
		if (sett->verbose > 1)
            std::cout <<"Chi2 in loop " <<i <<"is: " << chi2m << std::endl;
        if (sett->verbose > 1 )
            std::cout <<"scale in loop " <<i <<"is: " << scale << std::endl;
        if (sett->verbose > 1)
            std::cout <<"d values in loop " <<i <<"are: " << d[0] <<" " <<d[1] << std::endl;
        if (sett->verbose > 1)
            std::cout <<"ref value in loop " <<i <<"is: " << ref_val[i] << std::endl;
        
        //calculate chi2
        if (d[1] != 0 )
            chi2m += TMath::Power( d[0] - ( scale * ref_val[i] ) , 2 ) / TMath::Power(d[1],2);
       // std::cout << d[1] <<std::endl;
    }
    
    return chi2m;
}


double ShapeChi2::GetScale() {
    chi2 = 0;
    double chi2_t = 1E9;
    double scale = 1;
    double scale_tot = 1;
    double dscale = 0.5;
    int i = 0;
    bool walk = true; // defines if scale is increasing or decreasing
    do {
        scale_tot = scale_tot * scale;
        if (mode == 0)
            chi2 = minChi2(scale_tot);
        else
            chi2 = minChi2Sort(scale_tot);
        
        if (sett->verbose > 1)
            std::cout <<"New chi2 value: " << chi2 << "with factor " << scale_tot << endl;
        if (chi2 > chi2_t) { //new chi2 is larger -->change direction and decrese step size
            walk = !walk;
            dscale = 0.5 * dscale;
        }
        if (walk)
            scale = 1 + dscale;
        else
            scale = 1 - dscale;
        i++;
        if (i > 100) {
            std::cout <<"No convergence reached in chi2 fitting...giving up after 100 steps" <<std::endl;
            chi2 = -1;
            break;
        }
        chi2_t = chi2;
        //std::cout <<"Chi2 is loop " <<i<< " is " <<chi2 <<std::endl;
    } while (dscale > 0.005);
    
    if (sett->verbose )
        std::cout <<"Scaling factor: " << scale_tot<<endl;
    return scale_tot;
}

