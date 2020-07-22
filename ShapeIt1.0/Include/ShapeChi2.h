#ifndef SHAPECHI2_H
#define SHAPECHI2_H
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"
#include "../Include/ShapeGSF.h"



class ShapeChi2 {
    
private:
    ShapeGSF *gSF;
    ShapeGSF *gSF_ref;
    
    double scale = 1;            //the scaling factor
    ShapeSetting *sett;
    std::vector <double> ref_ene;
    std::vector <double> ref_val;
    std::vector <double> ref_dval;
    double chi2 = -1;                //result of chi2 minimization; -1 means that no minimum has been found, yet
    int mode = 0;                   //0: uing gSF for chi2; 1: using gSF_sort instead
public:
    ShapeChi2(ShapeGSF *gSF_t, ShapeSetting *sett_t);
    ShapeChi2(ShapeGSF *gSF_t, ShapeGSF *gSF_ref, ShapeSetting *sett_t);
    double GetScale();
    double minChi2(double scale);
    double minChi2Sort(double scale);
    double getChi2() {return chi2;}
};
#endif
