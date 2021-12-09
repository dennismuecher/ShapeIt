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

#ifndef SHAPETALYS_H
#define SHAPETALYS_H

#include <iostream>
#include "ShapeSetting.h"

class ShapeTalys {
    
private:
    static const int            nOfSpins = 9;             //number of spins in talys output file
    TGraphAsymmErrors*          rhoGraph;                        //pointer to the experimental level densities
    ShapeSetting*              sett;
    TGraph*                     denTotGraph;                //graph of total lev density
    TGraph*                     denTotGraphTrans;                //graph of total lev density transformed via ptable and ctable
    
    TGraph*                     denPartialGraph;            //graph of partial level density using spinLow and spinHigh boundaries
    TGraph*                     denPartialGraphTrans;            //graph of partial level density using spinLow and spinHigh boundaries, transformed via ptable and ctable
    
    TGraph*                     denSpinGraph[nOfSpins];     //graphs of level densities for each spin
    TGraph*                     denSpinGraphTrans[nOfSpins];     //graphs of level densities for each spin transformed via ptable and ctable
    vector<double>              p_energy;                         //vector of energies for each parity
    vector<double>              energy;                         //vector of energies
    vector<double>              p_densityTot;                  //total level density for each parity
    vector<double>              densityTot;                  //total level density for each parity
    vector<double>              densityPartial;              //partial level density
    vector<double>              p_densitySpin[nOfSpins];       //level density for each spin for each parity
    vector<double>              densitySpin[nOfSpins];       //level density for each spin for each parity
    bool                        parityFlag;                            //talys output contains separate level densities for parities (parityFlag = 1) or not (parityFlag 0)
    bool                        format;                         //set to zero if rows 2 and 3 are empty in talys output, otherwise set to 1
    int                         spinLow = 0;                            //lower spin limit for partial level density
    int                         spinHigh = 3;                           //upper spin limit for partial level density
    
    std::string                 talysOutFile;                   //talys output file
    std::string                 discreteLevelFile;                  //file with discrete levels
    
    double                      chi2_min = 1E5;                  //minimum chi2 fit result
    int                         nOfDegFreedom=1;                   //noumber of degrees of freedom for chi2 fit (=#data points -1)

    
public:
    
    ShapeTalys(ShapeSetting* p_sett, std::string p_talysOutFile, TGraphAsymmErrors* p_rhoGraph, bool p_parityFlag, bool p_format, double p_ptable, double p_ctable);
    double                      ptable = 0;
    double                      ctable = 0;
    double                      ptablePartial = 0;      //best fit result to partial level density
    double                      ctablePartial = 0;      //best fit result to partial level density
    int                         NewReadTree();
    TGraph*                     getDenTotGraph() {return denTotGraph;}
    TGraph*                     getDenPartialGraph() {return denPartialGraph;}
    TGraph*                     getDenPartialGraphTrans() {return denPartialGraphTrans;}
    TGraph*                     getDenSpinGraph(int spin) {return denSpinGraph[spin];}
    
    TH1F*                       discreteHist;                      //graph of discrete levels

    void                        SetPCTable(double m_ptable, double m_ctable);
    void                        SetPCTable();
    void                        ReadDiscrete();
    double                      GetChi2Partial(double lower_ene, double higher_ene);
    void                        Chi2PartialLoop(double lower_ene, double higher_ene);
    void                        BestFitPartial();



    
};
#endif
