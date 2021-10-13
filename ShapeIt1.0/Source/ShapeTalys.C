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

#include "../Include/ShapeTalys.h"

//consstructor
ShapeTalys::ShapeTalys(std::string p_talysOutFile, bool p_parityFlag, bool p_format, double p_ptable) {
    parityFlag = p_parityFlag;
    ptable = p_ptable;
    format = p_format;
    talysOutFile = p_talysOutFile;
    NewReadTree();
}


int ShapeTalys::NewReadTree()
{
    ifstream file(talysOutFile.c_str());
    //ptable=0;
    if (!file)
    {
        cerr << "cannot read the file"
        << strerror(errno) << endl;
        return -1;
    }
    
    float e, denTot;
    float den[nOfSpins];
    float a, sigma;
    if (format == 0) {
        while (file >> e >> denTot >>den[0] >>den[1] >>den[2] >>den[3]>>den[4] >>den[5]>>den[6] >>den[7] >>den[8] ) {
            p_energy.push_back(e);
            p_densityTot.push_back(denTot);
            for (int i =0; i < nOfSpins; i++)
                p_densitySpin[i].push_back(den[i]);
        }
    }
    else {
        while (file >> e >> a >> sigma >> denTot >>den[0] >>den[1] >>den[2] >>den[3]>>den[4] >>den[5]>>den[6] >>den[7] >>den[8] ) {
            p_energy.push_back(e);
            p_densityTot.push_back(denTot);
            for (int i =0; i < nOfSpins; i++)
                p_densitySpin[i].push_back(den[i]);
        }
    }
    
    //std::cout <<"Found " <<p_energy.size() <<" rows!" <<std::endl;
    
    // in case of files with two parities, add level densities, accordingly; otherwise multiply with 2
    if (parityFlag) {
        int n =(int) p_energy.size()/2;
        for (int i =0; i < n; i++) {
            densityTot.push_back(p_densityTot[i]+p_densityTot[i+n]);
            energy.push_back(p_energy[i] + ptable);
            double partial = 0;
            for (int j =0; j < nOfSpins; j++) {
                densitySpin[j].push_back(p_densitySpin[j][i]+p_densitySpin[j][i+n]);
                if (j >=spinLow && j <=spinHigh)
                    partial +=densitySpin[j][i];
            }
            densityPartial.push_back(partial);
        }
    }
    else {
        int n =(int) p_energy.size();
        for (int i =0; i < n; i++) {
            densityTot.push_back(2*p_densityTot[i]);
            energy.push_back(p_energy[i] + ptable);
            double partial = 0;
            for (int j =0; j < nOfSpins; j++) {
                densitySpin[j].push_back(2*p_densitySpin[j][i]);
                if (j >=spinLow && j <=spinHigh)
                    partial +=densitySpin[j][i];
            }
            densityPartial.push_back(partial);
        }
    }
    
    //create TGraphs
    int n =densityTot.size();
    denTotGraph = new TGraph(n, &energy[0],&densityTot[0]);
    denPartialGraph = new TGraph(n, &energy[0],&densityPartial[0]);

    for (int i = 0; i <nOfSpins; i++) {
        denSpinGraph[i] = new TGraph(n, &energy[0],&densitySpin[i][0]);
        denSpinGraph[i]->SetLineColor(i);
    }
    
    return 0;
}



