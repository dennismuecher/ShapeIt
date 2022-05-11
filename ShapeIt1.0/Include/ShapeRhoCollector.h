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

#ifndef SHAPERHOCOLLECTOR_H
#define SHAPERHOCOLLECTOR_H

#include <iostream>
#include "ShapeSetting.h"
#include "ShapeRho.h"


class ShapeRhoCollector {
    
private:
    ShapeSetting*           sett;                //the settings file
    ShapeRho*               rho;                  //stores the experimental rho values
    TGraphAsymmErrors*      rhoTrafo;             // transformed experimental rho
    std::vector <TGraphErrors*> litGraphs;        //any literature rho values can be added here
    ShapeTalys*             ldmodel[6];           //ldmodel graphs from settings file
    TGraph*                 ldmodelGraph[6];      //contain ldmodel results
    int                     ldmodelCounter = 0;     //number of ldmodels used (1-6)
    bool                    ldmodelDisplay[6];    //display switch for ldmodel graph
    bool                    ldmodelChi2Display[6]; //defines if the chi2 MC should be run
    bool                    rhoOrigDisplay = false;  //if true, display untransformed level density
    bool                    litGraphDisplay = false; //if ture, display literature level denisties
    TH1F*                   discreteLevel;         //graph of discrete levels
    TGraphErrors*           discBand;               //band of discrete level fits
    TGraph*                 discBandUpper;
    TGraph*                 discBandLower;
    TGraph*                 discFitGraph[30];
    double                  norm = 1;               //normalization of absolute level density
    TColor*                 color;
    void                    FillTrafoGraph();   //fills rhoTrafo and adds to collGraph
    void                    FillTalys();          //fills the ldmodel objects
    void                    ReadDiscrete();         //reads in discrete level data
    void                    ExpMC();                //runs MC to find exp error band
    void                    ExpNorm();              //finds normalization for level density
    void                    DiscreteBand();         //provides a fit of discrete levels
    double                  Chi2DiscreteExp(double eMin, double eMax, double norm); //chi2 between discrete levels and data
    string                  talysNames[6]={"ld1: Const. temp. + Fermi","ld2: Back-shift. Fermi","ld3: Gen. superf.","ld4: HFB, Skyrme","ld5: HFB, Skyrme + comb.","ld6: temp. dep. HFB, Gogny"};
    
    
public:
    
    ShapeRhoCollector(ShapeSetting* p_sett);    //constructor
    
    TMultiGraph*            collGraph;          //the multigraph containing all elements
    void                    Draw();             //draws all collector components
    void                    AddLitGraph(string litName);  //adds literature values stored in filename litName
   
};
#endif
