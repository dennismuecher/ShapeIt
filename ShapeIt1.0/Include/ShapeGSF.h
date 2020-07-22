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
	int peakID;				//peakID =1 for peak1, peakID=2 for peak2 	
}gSF_sor;


class ShapeGSF {
    
private:
    ShapeSetting *sett;
    ShapeMatrix *gSF_matrix;
    double B = 1;                           //transfomration gSF scaling
    double alpha = 0;                       //trnasformation gSF exponential slope
   
public:

    std::vector <gSF_sor> gSF_sort;      //vector of energies and gSF used for sorting and printing results
    std::vector <gSF_str> gSF;          //vector of gamma ray strength functions (egamma and value with error for each level)
    std::vector <double> p_ratio;       //vector of ratios of counts in peak1 / peak2
    ShapeGSF(ShapeSetting* setting_t, ShapeMatrix* matrix_t);
    ShapeGSF(ShapeSetting* setting_t);
    void GetgSF() {return gSF;}
	void ScaleAll(double scale);				//scale all results of gSF with factor scale
    void FillgSF();                         //calculates the gSF
    double getBgRatio(int bin, int level);      //return the ratio of peak to background counts
    void DoInterpol();
    
    int binRange[2];                        //first and last bin taken into account in analysis of gSF
    std::vector <double> InterpolValue(double ene);
    std::vector <double> InterpolValueSort(double ene);
    void ScaleSort(double factor);            //multiplies the gSF values in sort vector with scale
    void Scale(double factor);               //multiplies the gSF values with scale
    void Update();                    //refills the gSF values
    TGraphErrors* plotLit();            //plots literature data of gSF
    void gSF_Print();
    void gSF_Collect();
    TMultiGraph* gSF_SortHisto(bool colour);          //plots all gSF results using two colours, if colour = true
    TGraph* getRatioGraph();                //plot of peak area ratios
    void Transform(double B_t, double alpha_t);          //change trnasformation parameters and transform gSF via B*exp(alpha E_gamma)

    
};
#endif
