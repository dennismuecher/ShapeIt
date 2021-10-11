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

class ShapeTalys {
    
private:
    static const int            nOfSpins = 9;             //number of spins in talys output file
    TGraph*                     denTotGraph;                //graph of total lev density
    TGraph*                     denPartialGraph;            //graph of partial level density using spinLow and spinHigh boundaries
    TGraph*                     denSpinGraph[nOfSpins];     //graphs of level densities for each spin
    
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

public:
    
    ShapeTalys(std::string p_talysOutFile, bool p_parityFlag, bool p_format, double p_ptable);
    double                      ptable = 0;
    int NewReadTree();
    TGraph*                     getDenTotGraph() {return denTotGraph;}
    TGraph*                     getDenPartialGraph() {return denPartialGraph;}
    TGraph*                     getDenSpinGraph(int spin) {return denSpinGraph[spin];}

    
};
#endif
