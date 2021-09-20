/* ***********************************************************************
* Copyright (C) 2019-2021, Dennis Muecher.                               *
* All rights reserved.                                                   *
*                                                                        *
* This program is free software: you can redistribute it and/or modify   *
* it under the terms of the GNU General Public License as published by   *
* the Free Software Foundation, either version 3 of the License, or      *
* (at your option) any later version.                                    *
* You should have received a copy of the GNU General Public License      *
* along with this program. If not, see  http://www.gnu.org/licenses/.    *
*************************************************************************/

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
    TGraphAsymmErrors* getSmoothGraph()     {return (levGraphSmooth);}            
    double Slope(int i);
    void Merge();
    double Norm(TGraphErrors* T1, TGraphErrors* T2 );

    void DoInterpol();
    double getBgRatio(int bin, int level);

    void Print();
    void Transform(double t_B, double t_alpha);          //change trnasformation parameters and transform gSF via B*exp(alpha E_gamma)
    void Scale(int j, double factor);
    void ScaleAll(double factor);
    void Smooth(int res);
    void MergeAll();
    double getChi2(int n);

    std::vector<TGraphErrors*> levGraph_1;            //TGraphs to contain the gSF data for level1
    std::vector<TGraphErrors*> levGraph_2;            //TGraphs to contain the gSF data for level2
    std::vector<TGraphErrors*> levGraph;            //TGraphs to contain the gSF data for both levels
    TGraphErrors *levGraphAll;                      //contains all gSF data points from all iterations
    TGraphAsymmErrors* levGraphSmooth;          //contains all smoothed gSF data points from all iterations
    
};
#endif
