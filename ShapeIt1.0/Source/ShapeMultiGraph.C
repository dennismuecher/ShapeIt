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

#include "../Include/ShapeMultiGraph.h"

//fills the fillGraph
void ShapeMultiGraph::doFill(int min_index, int max_index) {
    std::cout <<"cool so far!" <<std::endl;
    
    fillGraph = new TGraphErrors();
    fillGraph->Set(0);
    TIter next(fGraphs);
    TObject *obj;
    vector <TGraph*> gr;
    gr.clear();

    double x_step = 0.01;
    
    double y_min, y_max;
    
    while ((obj = next())) {
        gr.push_back((TGraph*)obj);
    }
    
    //check if min_index and max_index are within range
    if (min_index <1) {
        std::cout <<"Error: min_index cannot me smaller than one!" <<std::endl;
        return;
    }
    else if (max_index > gr.size()){
        std::cout <<"Error: max_index cannot me larger than number of elements in multiGraph!" <<std::endl;
        return;
    }
        
    //for each x point, determine smallest and largest y value and fill into fillGraph
    for (double x = lower_ene; x < higher_ene; x+=x_step) {
        for (int i = min_index-1; i < max_index ; i++) {
            double y = gr[i]->Eval(x);
            if (i == min_index-1) {
                y_min = y;
                y_max = y;
            }
            else if (y < y_min)
                y_min = y;
            else if (y > y_max)
                y_max = y;
        }
        fillGraph->SetPoint(fillGraph->GetN(), x, (y_min+y_max)/2.);
        fillGraph->SetPointError(fillGraph->GetN()-1, 0, (y_max - y_min)/2.);
        //std::cout <<"y: " <<y_min << " " << y_max << std::endl;
    }
}
