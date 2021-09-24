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

#include "../Include/ShapeSetting.h"

ShapeSetting::ShapeSetting(void)
{
    ResetWidth();
}

//resets the width calibration to zero
void ShapeSetting::ResetWidth() {
    widthCal[0][0] = 0; widthCal[0][1] = 0;
    widthCal[1][0] = 0; widthCal[1][1] = 0;
}

//sets the background energies for level 1
void ShapeSetting::setBgEne1(double ene[4]) {
    
    for (int i =0; i < 4; i++)
        bgEne[0][i] = ene[i];
    
}

//sets the background energies for level 2
void ShapeSetting::setBgEne2(double ene[4]) {
    
    for (int i =0; i < 4; i++)
        bgEne[1][i] = ene[i];
}

//calculates the number of bins
int ShapeSetting::SizeToBin() {
    int diff = exiEne[1] - exiEne[0];
    int bins = (int) diff / exi_size[0];
    if (exi_size[0] * bins < diff)
        bins++;
    return bins;
}

//calculates the number of bins; the last bin might be smaller than sett->exi_size!
int ShapeSetting::SizeToBin(double size) {
    double diff = exiEne[1] - exiEne[0];
    
    int bins = (int) diff / size;
    if (exi_size[0] * bins < diff)
        bins++;
    return bins;
}

int ShapeSetting::BinToSize() {
    int size = (int) (exiEne[1] - exiEne[0] ) / nOfBins;
    return size;
}

int ShapeSetting::BinToSize(int n) {
    int size = (int) (exiEne[1] - exiEne[0] ) / n;
    return size;
}

void ShapeSetting::readEffi()
{
    
    //read file data into eGraph
    if (effiFileName == "") {
        std::cout << "No Efficiency File loaded!"<<std::endl;
        return NULL;
    }
    
    if (verbose)
        std::cout <<"\nReading Efficiency DATA... " <<endl;
    
    ifstream inp;
    inp.open(effiFileName.c_str());
    
    if (inp.is_open() ) {
        
        double e;
        double eff;
        
        //stores the minimum and maximum energies from the input file;
        double emin = -1;
        double emax = -1;
        
        int i = 0;
        eGraph = new TGraph();
        
        while ( !inp.eof() ) {
            inp >> e >> eff;
            eGraph->SetPoint(i,e,eff);
            i++;
            if (verbose)
                std::cout<< i<< " " << e << " "<< eff <<endl;
        }
        eGraph->SetMarkerStyle(4);
        eGraph->SetMarkerColor(kRed);
        eGraph->SetTitle("Energy dependend scaling factor; E_{#gamma} (keV); scaling factor");
    }
}

//returns an energy-dependend efficiency factor
double ShapeSetting::getEffCor(double ene, int level) {
	
    double c = 1;
    
    if (doEffi) {
        
        //find minimum and maximum x values of efficiency factors; add 50 keV margin 
        double xmin = TMath::MinElement(eGraph->GetN(),eGraph->GetX()) -50;
        double xmax = TMath::MaxElement(eGraph->GetN(),eGraph->GetX()) +50;
        
        if (verbose  > 1)
            std::cout <<"Minimum and maximum values in efficiency data file: " <<xmin <<" "<<xmax <<std::endl;
        
        //if ene not inside graph, return 1, otherwise return interpolated value
        if ( (ene >= xmin) && (ene <=xmax) )
            c = eGraph->Eval(ene);
        }
    //apply efficiency factor stored in eff_cor for level 2 in any case
    if (level == 2)
        c = c * eff_corr;
    return c;
}

void ShapeSetting::SaveSettings() {
    ofstream outfile;
    outfile.open (settFileName.c_str());
    outfile << dataFileName << "\n";
    outfile << osloFileName << "\n";
	outfile << effiFileName << "\n";
    outfile << "matrixName " << matrixName<<"\n";
    outfile << "MeV: " << MeV<<"\n";
    outfile << "mode " << mode << "\n";
    outfile << "doInterpol " << doInterpol<<"\n";
    outfile << "doOslo " << doOslo<<"\n";
    outfile << "doEffi " << doEffi<<"\n";
    outfile << "doAutoScale " << doAutoScale<<"\n";
	outfile << "colour " << colour<<"\n";
	outfile << "doGRF " << doGRF<<"\n";
    outfile << "displayAverage " << displayAvg<<"\n";
    outfile << "displaySingle " << displaySingle<<"\n";
    outfile << "doSlidingWindow " << doSlidingWindow <<"\n";
    outfile << "doBinVariation " << doBinVariation <<"\n";
    outfile << "doBackground " << doBackground <<"\n";
    outfile << "verbose " << verbose<<"\n";
    outfile << "gSF_norm " << gSF_norm<<"\n";
    outfile << "lit_norm " << lit_norm<<"\n";
    outfile << "lit_alpha " << lit_alpha<<"\n";
    outfile << "level1 " << levEne[0] <<" "<<levEne[1] <<"\n";
    outfile << "level1_2 " << levEne_2[0] <<" "<<levEne_2[1] <<"\n";
    outfile << "bg_level1 " << bgEne[0][0] <<" "<< bgEne[0][1] <<" "<< bgEne[0][2] <<" "<< bgEne[0][3] <<"\n";
    outfile << "bg_level2 " << bgEne[1][0] <<" "<< bgEne[1][1] <<" "<< bgEne[1][2] <<" "<< bgEne[1][3] <<"\n";
    outfile << "level2 " << levEne[2] <<" "<<levEne[3] <<"\n";
    outfile << "level2_2 " << levEne_2[2] <<" "<<levEne_2[3] <<"\n";
	outfile << "excitation " << exiEne[0] <<" "<<exiEne[1] <<"\n";
    outfile << "gamma_limits " << gammaEne[0] <<" "<<gammaEne[1] <<"\n";
    outfile << "excitation_bin_1 " << exi_size[0] <<"\n";
    outfile << "excitation_bin_2 " << exi_size[1] <<"\n";
    outfile << "nOfBins " << nOfBins <<"\n";
    outfile << "eff_corr " << eff_corr <<"\n";
    outfile << "minCounts " << minCounts <<"\n";
    outfile << "doWidthCal " << doWidthCal <<"\n";
    outfile << "widthCal " << widthCal[0][0] <<" "<< widthCal[0][1] <<" "<< widthCal[1][0] <<" "<< widthCal[1][1] <<"\n";
    outfile << "rhoFileName " << rhoFileName <<"\n";
	
    outfile.close();
    if (verbose)
        std::cout <<"Successfully saved Settings to file " << settFileName <<std::endl <<std::endl;
}

void ShapeSetting::ReadSettings() {
    if (verbose)
        std::cout <<"\nREADING FROM INPUT FILE: " <<endl;
    
    //clear fileNames; currently for rhoFileName, only, due to backwards compatibility
    rhoFileName.clear();
    
    std::ifstream inp (settFileName.c_str());
    if (inp.is_open()) {
        string word, line;
        //get the root matrix filename; this goes extra because of issues wih absolute paths containing white spaces
        getline(inp,line);
        dataFileName = line;
        //get the literature data filename; this goes extra because of issues wih absolute paths containing white spaces
        getline(inp,line);
        osloFileName = line;
        
		//get the efficiency data filename; due to backwards compatibility, check if this line contains the matrixName information (old file format)
        getline(inp,line);
        word.clear();
        istringstream isstr(line);
        isstr >> word;
		if (word == "matrixName" ) {
			isstr >> matrixName;
		}
		else {
			effiFileName = line;
		}
		
        //now get everything else
        while ( getline(inp,line) ) {
            word.clear();
            istringstream isstr(line);
            isstr >> word;
            if (word == "matrixName" ) isstr >> matrixName;
            if (word == "MeV:" ) isstr >> MeV;
            if (word == "mode" ) isstr >> mode ;
            if (word == "doInterpol" ) isstr >> doInterpol;
            if (word == "doOslo" ) isstr >> doOslo;
            if (word == "doEffi" ) isstr >> doEffi;
            if (word == "doAutoScale" ) isstr >> doAutoScale;
			if (word == "colour" ) isstr >> colour;
            if (word == "doGRF" ) isstr >> doGRF;
            if (word == "displayAverage" ) isstr >> displayAvg;
            if (word == "displaySingle" ) isstr >> displaySingle;
            if (word == "doSlidingWindow" ) isstr >> doSlidingWindow ;
            if (word == "doBinVariation" ) isstr >> doBinVariation ;
            if (word == "doBackground" ) isstr >> doBackground ;
            if (word == "verbose" ) isstr >> verbose;
            if (word == "gSF_norm" ) isstr >> gSF_norm;
            if (word == "lit_norm" ) isstr >> lit_norm;
            if (word == "lit_alpha" ) isstr >> lit_alpha;
            if (word == "level1" ) { isstr >> levEne[0]; isstr >>levEne[1];}
			if (word == "level1_2" ) { isstr >> levEne_2[0]; isstr >>levEne_2[1];}
            if (word == "level2" ) { isstr >> levEne[2]; isstr >>levEne[3];}
			if (word == "level2_2" ) { isstr >> levEne_2[2]; isstr >>levEne_2[3];}
            if (word == "bg_level1" ){ isstr >> bgEne[0][0]; isstr >> bgEne[0][1]; isstr >> bgEne[0][2]; isstr >> bgEne[0][3];}
            if (word == "bg_level2" ){ isstr >> bgEne[1][0]; isstr >> bgEne[1][1]; isstr >> bgEne[1][2]; isstr >> bgEne[1][3];}
            if (word == "excitation" ) { isstr >> exiEne[0]; isstr >>exiEne[1];}
            if (word == "gamma_limits" ) { isstr >> gammaEne[0]; isstr >>gammaEne[1];}
            if (word == "excitation_bin_1" ) isstr >> exi_size[0] ;
            if (word == "excitation_bin_2" ) isstr >> exi_size[1] ;
            if (word == "nOfBins" ) isstr >> nOfBins ;
            if (word == "eff_corr" ) isstr >> eff_corr ;
            if (word == "minCounts" ) isstr >> minCounts ;
            if (word == "doWidthCal" ) isstr >> doWidthCal ;
            if (word == "widthCal" ){ isstr >> widthCal[0][0]; isstr >> widthCal[0][1]; isstr >> widthCal[1][0]; isstr >> widthCal[1][1];}
            if (word == "rhoFileName" ) {
                string pName;
                rhoFileName.clear();
                while (isstr >>pName)
                    rhoFileName+=" "+pName;
            }
        }
        if (verbose)
                PrintSettings();
        
        if (doEffi)
            readEffi();
    }
}


void ShapeSetting::PrintSettings(){
    std::cout  << "root file name " << dataFileName<<"\n";
    std::cout  << "matrixName " << matrixName<<"\n";
    std::cout  << "Literature values gSF " << osloFileName<<"\n";
    std::cout  << "Literature values level density " << rhoFileName<<"\n";
    std::cout  << "Efficiency correction " << effiFileName<<"\n";
    std::cout  << "MeV: " << MeV<<"\n";
    std::cout  << "mode " << mode << "\n";
    std::cout  << "doInterpol " << doInterpol<<"\n";
    std::cout  << "doOslo " << doOslo<<"\n";
    std::cout  << "doEffi " << doEffi<<"\n";
    std::cout  << "doAutoScale " << doAutoScale<<"\n";
	std::cout  << "colour " << colour<<"\n";
    std::cout  << "doGRF " << doGRF<<"\n";
    std::cout  << "displayAverage " << displayAvg<<"\n";
    std::cout  << "displaySingle " << displaySingle<<"\n";
    std::cout  << "doSlidingWindow " << doSlidingWindow <<"\n";
    std::cout  << "doBinVariation " << doBinVariation <<"\n";
    std::cout  << "doBackground " << doBackground <<"\n";
    std::cout  << "verbose " << verbose<<"\n";
    std::cout  << "gSF norm " << gSF_norm<<"\n";
    std::cout  << "lit norm " << lit_norm<<"\n";
    std::cout  << "lit alpha " << lit_alpha<<"\n";
    std::cout  << "level1 " << levEne[0] <<" "<<levEne[1] <<"\n";
	std::cout  <<  "level1_2 " << levEne_2[0] <<" "<<levEne_2[1] <<"\n";
    std::cout  << "left background level1 " << bgEne[0][0] <<"-" << bgEne[0][1] <<"\n";
    std::cout  << "right background level1 " << bgEne[0][2] <<"-" << bgEne[0][3] <<"\n";
    std::cout  << "left background level2 " << bgEne[1][0] <<"-" << bgEne[1][1] <<"\n";
    std::cout  << "right background level2 " << bgEne[1][2] <<"-" << bgEne[1][3] <<"\n";
    std::cout  << "level2 " << levEne[2] <<" "<<levEne[3] <<"\n";
	std::cout  <<  "level2_2 " << levEne_2[2] <<" "<<levEne_2[3] <<"\n";
    std::cout  << "excitation " << exiEne[0] <<" "<<exiEne[1] <<"\n";
    std::cout  << "gamma limits " << gammaEne[0] <<" "<<gammaEne[1] <<"\n";
    std::cout  << "excitation_bin_1 " << exi_size[0] <<"\n";
    std::cout  << "excitation_bin_2 " << exi_size[1] <<"\n";
    std::cout  << "nOfBins " << nOfBins <<"\n";
    std::cout  << "eff_corr " << eff_corr <<"\n";
    std::cout  << "minCounts " << minCounts <<"\n";
    std::cout  << "doWidthCal " << doWidthCal <<"\n";
    std::cout  << "widthCal " << widthCal[0][0] <<" "<< widthCal[0][1] <<" "<< widthCal[1][0] <<" "<< widthCal[1][1] <<"\n";
  
}

