#ifndef SHAPEGSF_H
#define SHAPEGSF_H
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"
#include <TGraphErrors.h>

typedef struct gSF_struct{
    double egamma1, egamma2;       //gamma energy at which gSF is calcualted for level1 and level2
    double value1, value2;            //gamma ray strength of level 1 and level2
    double dvalue1, dvalue2;           //error gamma ray strength of level 1 and level 2
}gSF_str;

typedef struct gSF_sorting{
    double egamma;               //gamma energy
    double value;              //gamma ray strength
    double dvalue;             //error gamma ray strength
}gSF_sor;


class ShapeGSF {
    
private:
    ShapeSetting *sett;
    ShapeMatrix *gSF_matrix;
   
public:

    std::vector <gSF_sor> gSF_sort;      //vector of energies and gSF used for sorting and printing results
    std::vector <gSF_str> gSF;          //vector of gamma ray strength functions (egamma and value with error for each level)
    ShapeGSF(ShapeSetting* setting_t, ShapeMatrix* matrix_t);
    void GetgSF() {return gSF;}
    void FillgSF();                         //calculates the gSF
    double getBgRatio(int bin, int level);      //return the ratio of peak to background counts
    TGraphErrors* gSF_Histo();
    void DoInterpol();
    int binRange[2];                        //first and last bin taken into account in analysis of gSF
    std::vector <double> InterpolValue(double ene);
    void Scale(double factor);               //multiplies the gSF values with scale
    void Update();                    //refills the gSF values
    TGraphErrors* plotLit();            //plots literature data of gSF
    void gSF_Print();
    void gSF_Collect();
};
#endif
