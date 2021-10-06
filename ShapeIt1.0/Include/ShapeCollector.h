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

#ifndef SHAPECOLLECTOR_H
#define SHAPECOLLECTOR_H
#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"
#include "../Include/ShapeGSF.h"

#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>

#include <TObjArray.h>
#include <TH1.h>

class ShapeCollector {
    
private:
    ShapeSetting                    *m_sett;
    ShapeMatrix                     *m_matrix;
    std::vector<ShapeGSF*>          gSFCollector;       //gSF data from the experimental matrix
    ShapeGSF                        *litCollector;       //gSF data from literature
    TGraphErrors                    *gSFGraph;            //contains all gSF data
    TGraphAsymmErrors               *gSFGraphSmooth;       //contains all gSF data, smoothed 
    int                             kmax = 10;            //maximum number of steps in sliding window variation
    
    TMultiGraph*                    getMultGraph();
    double                          Norm(ShapeGSF* T1, ShapeGSF* T2);
    void                            Merge();
    void                            Smooth(int res);
    void                            NormCollect();


public:

    ShapeCollector(ShapeSetting* t_setting, ShapeMatrix* t_matrix);
    void                            Draw();
    void                            Print();
    void                            Collect();
};
#endif

