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


/*
 
This class is derived from TMultiGraph. Once TGraphs have been added, the function doFill(...) can be used to calculate the TGraphErrors object "fillGraph", which "fills" the area between any TGraphs stored in this object, offering to display an error band. The TGraphs which are considered can be specified with the min_index and max_index parameters of doFill(...): the TGraphs with index between min_index and max_index are considered. The index starts with 1, i.e. the first TGraph which has been added to this object has index 1, the second TGraph has index 2 etc.
 
 Once fillGraph has been calculated, one can either add it to this multigraph to display the error band (together with all the lines), or just use it as a normal TGraphsError object; In case multiple error bands are required, simply add the fillGraph object via e.g.
 
 Add((TGraphErrors*)fillGraph->Clone(),"3");

 to this (or any other) mulitgraph object.
 
 */

#ifndef SHAPEMULTIGRAPH_H
#define SHAPEMULTIGRAPH_H

#include <iostream>

class ShapeMultiGraph: public TMultiGraph {

private:
    double                  lower_ene = 1.0;
    double                  higher_ene = 10.0;
    bool                    visible = true; //if true, the lines will be drawn

public:
    TGraphErrors*           fillGraph;
    TGraph*                 fillGraphUpper;     //upper bounds of fillGraph; used to plot boundary lines
    TGraph*                 fillGraphLower;     //lower bounds of fillGraph; used to plot boundary lines
    void                    doFill(int min_index, int max_index);
    
    void                    getLowerEne() {return lower_ene;}
    void                    getHigherEne() {return higher_ene;}
    void                    setLowerEne(double p_lower_ene) {lower_ene = p_lower_ene;}
    void                    setHigherEne(double p_higher_ene) {higher_ene = p_higher_ene;}



};

#endif
