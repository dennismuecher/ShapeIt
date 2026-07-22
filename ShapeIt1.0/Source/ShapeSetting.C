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
    eGraph = nullptr;  // Initialize efficiency graph pointer
    for (int i = 0; i < 6; i++) {
        ldFileName[i] = "";
        pTable[i] = 0;
        cTable[i] = 0;
        parityFlag[i] = false;
        formatFlag[i] = false;
    }
    
    // Initialize peak position settings
    fixPeakPos[0] = false;
    fixPeakPos[1] = false;
    peakPos[0] = 0;
    peakPos[1] = 0;
    
    // Initialize doublet peak position settings
    fixDoubletPeakPos[0] = false;
    fixDoubletPeakPos[1] = false;
    doubletPeakPos[0] = 0;
    doubletPeakPos[1] = 0;
    
    // Initialize triplet peak position settings
    fixTripletPeakPos[0] = false;
    fixTripletPeakPos[1] = false;
    tripletPeakPos[0] = 0;
    tripletPeakPos[1] = 0;
}

ShapeSetting::~ShapeSetting(void)
{
    // Clean up efficiency graph if it exists
    if (eGraph) {
        delete eGraph;
        eGraph = nullptr;
    }
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
        return;
    }
    
    if (verbose)
        std::cout <<"\nReading Efficiency DATA... " <<std::endl;
    
    std::ifstream inp;
    inp.open(effiFileName.c_str());
    
    if (!inp.is_open()) {
        std::cout << "ERROR: Could not open efficiency file: " << effiFileName << std::endl;
        return;
    }
    
    // Delete old graph if it exists (prevent memory leak)
    if (eGraph) {
        delete eGraph;
        eGraph = nullptr;
    }
    
    double e;
    double eff;
    int i = 0;
    eGraph = new TGraph();
    
    std::string line;
    // Read line by line, skip comments
    while (std::getline(inp, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        if (iss >> e >> eff) {
            eGraph->SetPoint(i, e, eff);
            i++;
            if (verbose)
                std::cout << i << " " << e << " " << eff << std::endl;
        }
    }
    
    inp.close();
    
    if (i == 0) {
        std::cout << "WARNING: No data points read from efficiency file!" << std::endl;
        delete eGraph;
        eGraph = nullptr;
        doEffi = false;
        return;
    }
    
    eGraph->SetMarkerStyle(4);
    eGraph->SetMarkerColor(kRed);
    eGraph->SetTitle("Energy dependend scaling factor; E_{#gamma} (keV); scaling factor");
    
    if (verbose)
        std::cout << "Successfully loaded " << i << " efficiency correction points." << std::endl;
}

//returns an energy-dependend efficiency factor
double ShapeSetting::getEffCor(double ene, int level) {
	
    double c = 1;
    
    // Energy-dependent correction only applies to Level 2
    if (level == 2) {
        if (doEffi && eGraph && eGraph->GetN() > 0) {
            // Use energy-dependent correction from file (REPLACES constant eff_corr)
            double xmin = TMath::MinElement(eGraph->GetN(),eGraph->GetX()) -50;
            double xmax = TMath::MaxElement(eGraph->GetN(),eGraph->GetX()) +50;
            
            if (verbose  > 1)
                std::cout <<"Minimum and maximum values in efficiency data file: " <<xmin <<" "<<xmax <<std::endl;
            
            //if ene not inside graph, return 1, otherwise return interpolated value
            if ( (ene >= xmin) && (ene <=xmax) )
                c = eGraph->Eval(ene);
            // else: c = 1.0 (no correction if outside range)
        }
        else {
            // No file loaded or disabled: use constant efficiency correction
            c = eff_corr;
        }
    }
    // Level 1: no efficiency correction applied (always returns 1.0)
    
    return c;
}

void ShapeSetting::SaveSettings() {
    std::ofstream outfile;
    outfile.open (settFileName.c_str());
    if (isVersion2Format) {
        outfile << "# ShapeIt 2.0 settings file\n";
    }
    outfile << dataFileName << "\n";
    outfile << osloFileName << "\n";
	outfile << effiFileName << "\n";
    outfile << "matrixName " << matrixName<<"\n";
    outfile << "MeV: " << MeV<<"\n";
    outfile << "mode " << mode << "\n";
    outfile << "doInterpol " << doInterpol<<"\n";
    outfile << "doOslo " << doOslo<<"\n";
    outfile << "alphaLimit " << alphaLimit[0]<<" "<<alphaLimit[1] << "\n";
    outfile << "alphaIter " << alphaIter<<"\n";
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
    outfile << "lit_alpha_error " << lit_alpha_error[0] <<" "<<lit_alpha_error[1] <<"\n";
    outfile << "level1 " << levEne[0] <<" "<<levEne[1] <<"\n";
    outfile << "level1_2 " << levEne_2[0] <<" "<<levEne_2[1] <<"\n";
    outfile << "level1_3 " << levEne_3[0] <<" "<<levEne_3[1] <<"\n";
    outfile << "doDoublet1 " << doDoublet[0] <<"\n";
    outfile << "doTriplet1 " << doTriplet[0] <<"\n";
    outfile << "fixDoubletWidth1 " << fixDoubletWidth[0] <<"\n";
    outfile << "fixTripletWidth1 " << fixTripletWidth[0] <<"\n";
    outfile << "bg_level1 " << bgEne[0][0] <<" "<< bgEne[0][1] <<" "<< bgEne[0][2] <<" "<< bgEne[0][3] <<"\n";
    outfile << "bg_level2 " << bgEne[1][0] <<" "<< bgEne[1][1] <<" "<< bgEne[1][2] <<" "<< bgEne[1][3] <<"\n";
    outfile << "level2 " << levEne[2] <<" "<<levEne[3] <<"\n";
    outfile << "level2_2 " << levEne_2[2] <<" "<<levEne_2[3] <<"\n";
    outfile << "level2_3 " << levEne_3[2] <<" "<<levEne_3[3] <<"\n";
    outfile << "doDoublet2 " << doDoublet[1] <<"\n";
    outfile << "doTriplet2 " << doTriplet[1] <<"\n";
    outfile << "fixDoubletWidth2 " << fixDoubletWidth[1] <<"\n";
    outfile << "fixTripletWidth2 " << fixTripletWidth[1] <<"\n";
	outfile << "excitation " << exiEne[0] <<" "<<exiEne[1] <<"\n";
    outfile << "excitationMCLimits " << exiEneMC[0] <<" "<<exiEneMC[1] <<"\n";
    outfile << "excitation_bin_1 " << exi_size[0] <<"\n";
    outfile << "excitation_bin_2 " << exi_size[1] <<"\n";
    outfile << "nOfBins " << nOfBins <<"\n";
    outfile << "eff_corr " << eff_corr <<"\n";
    outfile << "minCounts " << minCounts <<"\n";
    outfile << "doWidthCal " << doWidthCal <<"\n";
    outfile << "widthCal " << widthCal[0][0] <<" "<< widthCal[0][1] <<" "<< widthCal[1][0] <<" "<< widthCal[1][1] <<"\n";
    outfile << "fixPeakPos1 " << fixPeakPos[0] <<"\n";
    outfile << "peakPos1 " << peakPos[0] <<"\n";
    outfile << "fixPeakPos2 " << fixPeakPos[1] <<"\n";
    outfile << "peakPos2 " << peakPos[1] <<"\n";
    outfile << "fixDoubletPeakPos1 " << fixDoubletPeakPos[0] <<"\n";
    outfile << "doubletPeakPos1 " << doubletPeakPos[0] <<"\n";
    outfile << "fixDoubletPeakPos2 " << fixDoubletPeakPos[1] <<"\n";
    outfile << "doubletPeakPos2 " << doubletPeakPos[1] <<"\n";
    outfile << "fixTripletPeakPos1 " << fixTripletPeakPos[0] <<"\n";
    outfile << "tripletPeakPos1 " << tripletPeakPos[0] <<"\n";
    outfile << "fixTripletPeakPos2 " << fixTripletPeakPos[1] <<"\n";
    outfile << "tripletPeakPos2 " << tripletPeakPos[1] <<"\n";
    outfile << "rhoFileName " << rhoFileName <<"\n";
    outfile << "rhoScale " << rhoScale <<"\n";
    outfile << "spinLow " << spinLow <<"\n";
    outfile << "spinHigh " << spinHigh <<"\n";

    outfile << "discreteLevelFile " << discreteLevelFile << " " << discreteBins << "\n";
    
    for (int i = 0; i < 6; i++) {
        if (ldFileName[i]!="")
            outfile << "ldmodel" <<i+1<<" " << ldFileName[i] <<" " <<pTable[i] <<" "<< cTable[i]<< " " << parityFlag[i] << " " << formatFlag[i] << "\n";
    }
    outfile.close();
    if (verbose)
        std::cout <<"Successfully saved Settings to file " << settFileName <<std::endl <<std::endl;
}

void ShapeSetting::ReadSettings() {
    if (verbose)
        std::cout <<"\nREADING FROM INPUT FILE: " <<std::endl;
    
    //clear fileNames; currently for rhoFileName, only, due to backwards compatibility
    rhoFileName.clear();
    
    std::ifstream inp (settFileName.c_str());
    if (inp.is_open()) {
        std::string word, line;
        //get the root matrix filename; this goes extra because of issues wih absolute paths containing white spaces
        getline(inp,line);
        
        // Detect file format version
        if (line.length() > 0 && line[0] == '#') {
            // Version 2.0 format detected
            isVersion2Format = true;
            // Skip all comment lines and read the actual data filename
            while (line.length() > 0 && line[0] == '#') {
                getline(inp,line);
            }
            dataFileName = line;
        } else {
            // Old format (ShapeIt 1.0)
            isVersion2Format = false;
            dataFileName = line;
        }
        
        if (verbose) {
            std::cout << "Settings file format: " << (isVersion2Format ? "ShapeIt 2.0" : "ShapeIt 1.0") << std::endl;
        }
        //get the literature data filename; this goes extra because of issues wih absolute paths containing white spaces
        getline(inp,line);
        osloFileName = line;
        
		//get the efficiency data filename; due to backwards compatibility, check if this line contains the matrixName information (old file format)
        getline(inp,line);
        word.clear();
        std::istringstream isstr(line);
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
            std::istringstream isstr(line);
            isstr >> word;
            if (word == "matrixName" ) isstr >> matrixName;
            if (word == "MeV:" ) isstr >> MeV;
            if (word == "mode" ) isstr >> mode ;
            if (word == "doInterpol" ) isstr >> doInterpol;
            if (word == "doOslo" ) isstr >> doOslo;
            if (word == "alphaLimit" ) { isstr >> alphaLimit[0]; isstr >>alphaLimit[1];}
            if (word == "alphaIter" ) isstr >> alphaIter;
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
            if (word == "lit_alpha_error" ) { isstr >> lit_alpha_error[0]; isstr >> lit_alpha_error[1];}
            if (word == "level1" ) { isstr >> levEne[0]; isstr >>levEne[1];}
			if (word == "level1_2" ) { isstr >> levEne_2[0]; isstr >>levEne_2[1];}
			if (word == "level1_3" ) { isstr >> levEne_3[0]; isstr >>levEne_3[1];}
            if (word == "doDoublet1" ) isstr >> doDoublet[0];
            if (word == "doTriplet1" ) isstr >> doTriplet[0];
            if (word == "fixDoubletWidth1" ) isstr >> fixDoubletWidth[0];
            if (word == "fixTripletWidth1" ) isstr >> fixTripletWidth[0];
            if (word == "level2" ) { isstr >> levEne[2]; isstr >>levEne[3];}
			if (word == "level2_2" ) { isstr >> levEne_2[2]; isstr >>levEne_2[3];}
			if (word == "level2_3" ) { isstr >> levEne_3[2]; isstr >>levEne_3[3];}
            if (word == "doDoublet2" ) isstr >> doDoublet[1];
            if (word == "doTriplet2" ) isstr >> doTriplet[1];
            if (word == "fixDoubletWidth2" ) isstr >> fixDoubletWidth[1];
            if (word == "fixTripletWidth2" ) isstr >> fixTripletWidth[1];
            if (word == "bg_level1" ){ isstr >> bgEne[0][0]; isstr >> bgEne[0][1]; isstr >> bgEne[0][2]; isstr >> bgEne[0][3];}
            if (word == "bg_level2" ){ isstr >> bgEne[1][0]; isstr >> bgEne[1][1]; isstr >> bgEne[1][2]; isstr >> bgEne[1][3];}
            if (word == "excitation" ) { isstr >> exiEne[0]; isstr >>exiEne[1];}
            if (word == "excitationMCLimits" ) { isstr >> exiEneMC[0]; isstr >>exiEneMC[1];}
            if (word == "excitation_bin_1" ) isstr >> exi_size[0] ;
            if (word == "excitation_bin_2" ) isstr >> exi_size[1] ;
            if (word == "nOfBins" ) isstr >> nOfBins ;
            if (word == "eff_corr" ) isstr >> eff_corr ;
            if (word == "minCounts" ) isstr >> minCounts ;
            if (word == "doWidthCal" ) isstr >> doWidthCal ;
            if (word == "widthCal" ){ isstr >> widthCal[0][0]; isstr >> widthCal[0][1]; isstr >> widthCal[1][0]; isstr >> widthCal[1][1];}
            if (word == "fixPeakPos1" ) isstr >> fixPeakPos[0];
            if (word == "peakPos1" ) isstr >> peakPos[0];
            if (word == "fixPeakPos2" ) isstr >> fixPeakPos[1];
            if (word == "peakPos2" ) isstr >> peakPos[1];
            if (word == "fixDoubletPeakPos1" ) isstr >> fixDoubletPeakPos[0];
            if (word == "doubletPeakPos1" ) isstr >> doubletPeakPos[0];
            if (word == "fixDoubletPeakPos2" ) isstr >> fixDoubletPeakPos[1];
            if (word == "doubletPeakPos2" ) isstr >> doubletPeakPos[1];
            if (word == "fixTripletPeakPos1" ) isstr >> fixTripletPeakPos[0];
            if (word == "tripletPeakPos1" ) isstr >> tripletPeakPos[0];
            if (word == "fixTripletPeakPos2" ) isstr >> fixTripletPeakPos[1];
            if (word == "tripletPeakPos2" ) isstr >> tripletPeakPos[1];
            if (word == "discreteLevelFile" ) { isstr >>discreteLevelFile; isstr >>discreteBins;}
            if (word == "rhoScale" ) isstr >> rhoScale ;
            if (word == "spinLow" ) isstr >> spinLow ;
            if (word == "spinHigh" ) isstr >> spinHigh ;

            if (word == "rhoFileName" ) {
                std::string pName;
                rhoFileName.clear();
                while (isstr >>pName)
                    //rhoFileName+=" "+pName;
                    rhoFileName+=pName;
            }
            for (int i = 0; i < 6; i++) {
                std::string ldname = "ldmodel"+to_string(i+1);
                if (word == ldname) {
                    isstr >>ldFileName[i];
                    isstr >>pTable[i];
                    isstr >>cTable[i];
                    isstr >>parityFlag[i];
                    isstr >>formatFlag[i];
                }
            }
        }
        
        // Backward compatibility: Calculate peak positions from levEne if not set
        bool needsUpdate = false;
        if (peakPos[0] == 0 && levEne[0] != 0 && levEne[1] != 0) {
            peakPos[0] = (levEne[0] + levEne[1]) / 2.0;
            std::cout << "\n*** ShapeIt 1.0 compatibility mode ***" << std::endl;
            std::cout << "Peak position 1 (Level 1) was calculated from level energy range: " << peakPos[0] << " keV" << std::endl;
            needsUpdate = true;
        }
        
        if (peakPos[1] == 0 && levEne[2] != 0 && levEne[3] != 0) {
            peakPos[1] = (levEne[2] + levEne[3]) / 2.0;
            if (!needsUpdate) {
                std::cout << "\n*** ShapeIt 1.0 compatibility mode ***" << std::endl;
            }
            std::cout << "Peak position 2 (Level 2) was calculated from level energy range: " << peakPos[1] << " keV" << std::endl;
            needsUpdate = true;
        }
        
        if (doubletPeakPos[0] == 0 && levEne_2[0] != 0 && levEne_2[1] != 0 && doDoublet[0]) {
            doubletPeakPos[0] = (levEne_2[0] + levEne_2[1]) / 2.0;
            if (!needsUpdate) {
                std::cout << "\n*** ShapeIt 1.0 compatibility mode ***" << std::endl;
            }
            std::cout << "Doublet peak position 1 (Level 1, Peak 2) was calculated from level energy range: " << doubletPeakPos[0] << " keV" << std::endl;
            needsUpdate = true;
        }
        
        if (doubletPeakPos[1] == 0 && levEne_2[2] != 0 && levEne_2[3] != 0 && doDoublet[1]) {
            doubletPeakPos[1] = (levEne_2[2] + levEne_2[3]) / 2.0;
            if (!needsUpdate) {
                std::cout << "\n*** ShapeIt 1.0 compatibility mode ***" << std::endl;
            }
            std::cout << "Doublet peak position 2 (Level 2, Peak 2) was calculated from level energy range: " << doubletPeakPos[1] << " keV" << std::endl;
            needsUpdate = true;
        }
        
        if (tripletPeakPos[0] == 0 && levEne_3[0] != 0 && levEne_3[1] != 0 && doTriplet[0]) {
            tripletPeakPos[0] = (levEne_3[0] + levEne_3[1]) / 2.0;
            if (!needsUpdate) {
                std::cout << "\n*** ShapeIt 1.0 compatibility mode ***" << std::endl;
            }
            std::cout << "Triplet peak position 1 (Level 1, Peak 3) was calculated from level energy range: " << tripletPeakPos[0] << " keV" << std::endl;
            needsUpdate = true;
        }
        
        if (tripletPeakPos[1] == 0 && levEne_3[2] != 0 && levEne_3[3] != 0 && doTriplet[1]) {
            tripletPeakPos[1] = (levEne_3[2] + levEne_3[3]) / 2.0;
            if (!needsUpdate) {
                std::cout << "\n*** ShapeIt 1.0 compatibility mode ***" << std::endl;
            }
            std::cout << "Triplet peak position 2 (Level 2, Peak 3) was calculated from level energy range: " << tripletPeakPos[1] << " keV" << std::endl;
            needsUpdate = true;
        }
        
        if (needsUpdate) {
            std::cout << "\nThese peak positions will be saved to the settings file when you save your settings." << std::endl;
            std::cout << "Future versions will use the peak position values directly.\n" << std::endl;
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
    std::cout  << "Discrete Level File " << discreteLevelFile << " " << discreteBins <<"\n";

    for (int i =0; i < 6; i++) {
        if (ldFileName[i]!="") {
            std::cout  << "ldmodel" <<i+1 <<" " << ldFileName[i]<<"\n";
            std::cout  << "ptable " << pTable[i]<<"\n";
            std::cout  << "ctable " << cTable[i]<<"\n";
            std::cout  << "parityFlag " << parityFlag[i]<<"\n";
            std::cout  << "formatFlag " << formatFlag[i]<<"\n";


        }

    }
    std::cout  << "Efficiency correction " << effiFileName<<"\n";
    std::cout  << "MeV: " << MeV<<"\n";
    std::cout  << "mode " << mode << "\n";
    std::cout  << "doInterpol " << doInterpol<<"\n";
    std::cout  << "doOslo " << doOslo<<"\n";
    std::cout  << "alphaLimit " << alphaLimit[0]<<" " <<alphaLimit[1]<<"\n";
    std::cout  << "alphaIter " << alphaIter<<"\n";
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
    std::cout  << "lit alpha_error " << lit_alpha_error[0] << " "<<lit_alpha_error[1] <<"\n";
    std::cout  << "level1 " << levEne[0] <<" "<<levEne[1] <<"\n";
	std::cout  <<  "level1_2 " << levEne_2[0] <<" "<<levEne_2[1] <<"\n";
	std::cout  <<  "level1_3 " << levEne_3[0] <<" "<<levEne_3[1] <<"\n";
    std::cout  <<  "doDoublet1 " << doDoublet[0] <<"\n";
    std::cout  <<  "doTriplet1 " << doTriplet[0] <<"\n";
    std::cout  <<  "fixDoubletWidth1 " << fixDoubletWidth[0] <<"\n";
    std::cout  <<  "fixTripletWidth1 " << fixTripletWidth[0] <<"\n";
    std::cout  << "left background level1 " << bgEne[0][0] <<"-" << bgEne[0][1] <<"\n";
    std::cout  << "right background level1 " << bgEne[0][2] <<"-" << bgEne[0][3] <<"\n";
    std::cout  << "left background level2 " << bgEne[1][0] <<"-" << bgEne[1][1] <<"\n";
    std::cout  << "right background level2 " << bgEne[1][2] <<"-" << bgEne[1][3] <<"\n";
    std::cout  << "level2 " << levEne[2] <<" "<<levEne[3] <<"\n";
	std::cout  <<  "level2_2 " << levEne_2[2] <<" "<<levEne_2[3] <<"\n";
	std::cout  <<  "level2_3 " << levEne_3[2] <<" "<<levEne_3[3] <<"\n";
    std::cout  <<  "doDoublet2 " << doDoublet[1] <<"\n";
    std::cout  <<  "doTriplet2 " << doTriplet[1] <<"\n";
    std::cout  <<  "fixDoubletWidth2 " << fixDoubletWidth[1] <<"\n";
    std::cout  <<  "fixTripletWidth2 " << fixTripletWidth[1] <<"\n";
    std::cout  << "excitation " << exiEne[0] <<" "<<exiEne[1] <<"\n";
    std::cout  << "excitation range for lower energies in MC " << exiEneMC[0] <<" "<<exiEneMC[1] <<"\n";
    std::cout  << "excitation_bin_1 " << exi_size[0] <<"\n";
    std::cout  << "excitation_bin_2 " << exi_size[1] <<"\n";
    std::cout  << "nOfBins " << nOfBins <<"\n";
    std::cout  << "eff_corr " << eff_corr <<"\n";
    std::cout  << "minCounts " << minCounts <<"\n";
    std::cout  << "rhoScale " << rhoScale <<"\n";
    std::cout  << "spinLow " << spinLow <<"\n";
    std::cout  << "spinHigh " << spinHigh <<"\n";
    std::cout  << "doWidthCal " << doWidthCal <<"\n";
    std::cout  << "widthCal " << widthCal[0][0] <<" "<< widthCal[0][1] <<" "<< widthCal[1][0] <<" "<< widthCal[1][1] <<"\n";
    std::cout  << "fixPeakPos1 " << fixPeakPos[0] <<"\n";
    std::cout  << "peakPos1 " << peakPos[0] <<"\n";
    std::cout  << "fixPeakPos2 " << fixPeakPos[1] <<"\n";
    std::cout  << "peakPos2 " << peakPos[1] <<"\n";
    std::cout  << "fixDoubletPeakPos1 " << fixDoubletPeakPos[0] <<"\n";
    std::cout  << "doubletPeakPos1 " << doubletPeakPos[0] <<"\n";
    std::cout  << "fixDoubletPeakPos2 " << fixDoubletPeakPos[1] <<"\n";
    std::cout  << "doubletPeakPos2 " << doubletPeakPos[1] <<"\n";
    std::cout  << "fixTripletPeakPos1 " << fixTripletPeakPos[0] <<"\n";
    std::cout  << "tripletPeakPos1 " << tripletPeakPos[0] <<"\n";
    std::cout  << "fixTripletPeakPos2 " << fixTripletPeakPos[1] <<"\n";
    std::cout  << "tripletPeakPos2 " << tripletPeakPos[1] <<"\n";
  
}

