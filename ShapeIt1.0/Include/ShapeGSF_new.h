#ifndef SHAPEGSF_H
#define SHAPEGSF_H
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TObjArray.h>

//#include <TContainer.h>


class ShapeGSF_new {
    
private:
    ShapeSetting *m_sett;
    ShapeMatrix *m_matrix;

    double m_B = 1;                           //transfomration gSF scaling
    double m_alpha = 0;                       //trnasformation gSF exponential slope
   
public:

    TMultiGraph* multGraph;                         //the graph containing all gSF data of all levGraph
    TGraphErrors* litGraph;
    ShapeGSF_new(ShapeSetting* t_setting, ShapeMatrix* t_matrix);
    
    void ReadLit();
    void FillgSF();                                 //calculates the gSF values from the m_gSF_matrix
    TMultiGraph* getMultGraph();                        //returns the multGrap
            
    double Slope(int i);
    void Merge();

    void DoInterpol(int m_i);
    double getBgRatio(int bin, int level);

    void Print();
    void Transform(double t_B, double t_alpha);          //change trnasformation parameters and transform gSF via B*exp(alpha E_gamma)
    void Scale(double factor);
  
    std::vector<TGraphErrors*> levGraph_1;            //TGraphs to contain the gSF data for level1
    std::vector<TGraphErrors*> levGraph_2;            //TGraphs to contain the gSF data for level2
    std::vector<TGraphErrors*> levGraph;            //TGraphs to contain the gSF data for both levels
};
#endif
