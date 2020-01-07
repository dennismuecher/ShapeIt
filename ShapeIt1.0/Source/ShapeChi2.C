

#include "../Include/ShapeChi2.h"
//the gSF values of gSF_fit are going to be scaled to those of gSF_ref
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
}

double ShapeChi2::getChi2(double scale) {
    double chi2 = 0;
    std::vector <double> d;
    for (int i = 0; i < ref_ene.size() - 2; i++) {
        d = gSF->InterpolValue(ref_ene[i]);
        if (d[0] > 0.0001) {
            std::cout <<"Unusual gSF detected in Chi2 in loop " <<i <<std::endl;
         sett->verbose = 1;
        }
        //this is a chi2, weighted by the relative error for each gSF value;
        
        if (d[0] < 0 ) {
            std::cout <<"negative gSF value zero found in Chi2 routine "<<std::endl;
            sett->verbose = 1;
            continue;
        }
        if (d[0] == 0 ) {
            std::cout <<"zero gSF value zero found in Chi2 routine "<<std::endl;
            continue;
        }
        if (d[1] == 0 ) {
            std::cout <<"zero gSF error value zero found in Chi2 routine "<<std::endl;
            continue;
        }
        if (sett->verbose > 1)
            std::cout <<"Chi2 in loop " <<i <<"is: " << chi2 << std::endl;
        if (sett->verbose > 1 )
            std::cout <<"scale in loop " <<i <<"is: " << scale << std::endl;
        if (sett->verbose > 1)
            std::cout <<"d values in loop " <<i <<"are: " << d[0] <<" " <<d[1] << std::endl;
        if (sett->verbose > 1)
            std::cout <<"ref value in loop " <<i <<"is: " << ref_val[i] << std::endl;
        //chi2 += d[0]/d[1] * (d[0] - ( scale * ref_val[i] ) );
        chi2 += (d[0] - ( scale * ref_val[i] ) );
    }
    
    return chi2;
}


double ShapeChi2::GetScale() {
    double chi2 = 0;
    double scale = 1;
    double scale_tot = 1;
    double dscale = 0.5;
    int i = 0;
    int walk =1; //1: positive direction; -1: negative direction
    do {
        scale_tot = scale_tot * scale;
        chi2 = getChi2(scale_tot);
        if (sett->verbose > 1)
            std::cout <<"New chi2 value: " << chi2 << "with factor " << scale_tot << endl;
        if (chi2 > 0) { // ref values larger than fit values
            if (walk == -1) {
                dscale = dscale * 0.5;
                walk = 1;
            }
            scale = 1 + dscale;
        }
        else {
            if (walk == 1) {
                dscale = dscale * 0.5;
                walk = -1;
            }
            scale = 1 - dscale;
        }
        i++;
        if (i > 50) {
            std::cout <<"No convergence reached in chi2 fitting...giving up after 50 steps" <<std::endl;
            break;
        }
    } while (dscale > 0.005);
     if (sett->verbose )
         std::cout <<"Scaling factor: " << scale_tot<<endl;
    return scale_tot;
}
