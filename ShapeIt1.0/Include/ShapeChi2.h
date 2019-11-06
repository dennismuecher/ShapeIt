#ifndef SHAPECHI2_H
#define SHAPECHI2_H
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"
#include "../Include/ShapeGSF.h"



class ShapeChi2 {
    
private:
    ShapeGSF *gSF;
    double scale = 1;            //the scaling factor
    ShapeSetting *sett;
    std::vector <double> ref_ene;
    std::vector <double> ref_val;

public:
    ShapeChi2(ShapeGSF *gSF_t, ShapeSetting *sett_t);
    double GetScale();
    double getChi2(double scale);
};
#endif
