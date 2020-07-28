#ifndef SHAPEGSF_H
#define SHAPEGSF_H
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"
#include <TGraphErrors.h>

class ShapeGSF {
    
private:
    ShapeSetting *m_sett;
    ShapeMatrix *m_gsfMatrix;

    int nOfLevel = 1;                           // number of levels used in analysis
    double m_B = 1;                           //transfomration gSF scaling
    double m_alpha = 0;                       //trnasformation gSF exponential slope
   
public:

    
    ShapeGSF(ShapeSetting* t_setting, ShapeMatrix* t_matrix);
    ShapeGSF(ShapeSetting* t_setting);
    
    std::vector<TGraphErrors*> levGraph;            //TGraphs to contain the gSF data for level1 and level2
    
    TGraphErrors* mergeGraph;
    
    void FillgSF();                                 //calculates the gSF values from the m_gSF_matrix
    
    void DoInterpol();
    
    void Update();                                    //refills the gSF values
    
    void gsfPrint();
    
    void Transform(double t_B, double t_alpha);          //change trnasformation parameters and transform gSF via B*exp(alpha E_gamma)

    
};
#endif
