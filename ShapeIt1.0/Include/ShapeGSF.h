#ifndef SHAPEGSF_H
#define SHAPEGSF_H
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>

#include <TObjArray.h>
#include <TH1.h>


//#include <TContainer.h>


class ShapeGSF {
    
private:
    ShapeSetting *m_sett;
    ShapeMatrix *m_matrix;
    double scale_bak;                           //stores the scaling value applied to gSF
    double m_B = 1;                           //transfomration gSF scaling
    double m_alpha = 0;                       //trnasformation gSF exponential slope
   
public:

    TGraphErrors* litGraph;
    ShapeGSF(ShapeSetting* t_setting, ShapeMatrix* t_matrix);
    
    void ReadLit();
    void FillgSF();                                 //calculates the gSF values from the m_gSF_matrix
    TMultiGraph* getMultGraph();                        //returns the multGrap
            
    double Slope(int i);
    void Merge();
    double Norm(TGraphErrors* T1, TGraphErrors* T2 );

    void DoInterpol();
    double getBgRatio(int bin, int level);

    void Print();
    void Transform(double t_B, double t_alpha);          //change trnasformation parameters and transform gSF via B*exp(alpha E_gamma)
    void Scale(int j, double factor);
    void ScaleAll(double factor);
    TGraphAsymmErrors* Smooth(int res);
    void MergeAll();

    std::vector<TGraphErrors*> levGraph_1;            //TGraphs to contain the gSF data for level1
    std::vector<TGraphErrors*> levGraph_2;            //TGraphs to contain the gSF data for level2
    std::vector<TGraphErrors*> levGraph;            //TGraphs to contain the gSF data for both levels
    TGraphErrors *levGraphAll;                      //contains all gSF data points from all iterations
    
};
#endif
