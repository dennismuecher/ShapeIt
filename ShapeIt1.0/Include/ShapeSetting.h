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


#ifndef SHAPESETTING_H
#define SHAPESETTING_H

class ShapeSetting {

private:
	


public:
	ShapeSetting();                         //the constructor
    
	std::string dataFileName; 				//name of file containing matrix
	std::string matrixName;  				//name of matrix in dataFileName
	std::string osloFileName = "";				 //name of file containing reference data for gSF
    std::string rhoFileName = "";                 //name of file containing reference data for level density
	std::string effiFileName = "";				 //name of file containing efficiency correction
    std::string settFileName = "";                  //name of settings save file
    TFile *settFile;                         //the file used to safe the current settings
    int mode = 1;                           //mode: 1 = integration; 2 = autofit
	
    double MeV = 1.;
    bool doOslo = false;					//if true, plot reference data
    bool doEffi = false;					//if true, apply efficiency correction
  
	bool doInterpol = true;                 //if true, do interpolation to get gSF
    bool doBackground = true;               //if true, do background subtraction
    bool doWidthCal = false;                 //if true, use results from width calibrations
    bool colour = true;						//if true, plot gSF points in two different colours for peak1 and peak2
	bool doGRF = false;                      //if true, display fit of gSF via giant resonance formula from RAINIER
    bool displayAvg = false;                    //sets if average gSF is displayed
    bool displaySingle = true;                    //sets if single (both) gSF data is displayed
	int verbose = 0;					//prints extra information if true; verbose =1 is basic information; verbose > 1 gives full detail
	double gSF_norm = 1 ;				//normalization factor for gSF
    double lit_norm = 1 ;                //normalization factor for Oslo literature value
    double lit_alpha = 0 ;                //transformation "slope" parameter for literature values

	double levEne[4] = {0, 0, 0, 0};						//energies for level 1 and level 2 (upper, lower)	
	double levEne_2[4] = {0, 0, 0, 0};						//energies for doublet peaks level 1 and level 2 (upper, lower)	
	double bgEne[2][4];                      //background regions for level 1 and level 2
    double bgWidth = 100;                   //width of each background window
    double exiEne[2] = {0, 0};						//lower and upper excitation energy considered
    double gammaEne[2] = {0, 0};                        //lower and upper gamma energy considered
    double exi_size[2] = {0, 0};					//bin size for shape method
    double widthCal[2][2];                  //parameters for width calibration
    int nOfBins = 1;
    double eff_corr = 1;					//correction factor for level 2;
    int minCounts = 0;                      //minimum number of counts required in a bin to be considered
	
	void SetMeV(bool b) {MeV = (b) ? 1000 : 1;}		//set MeV=1000 in case of b true
    void SetFileName(string str) {dataFileName = str;}
    std::string GetFileName() {return dataFileName;}
    void PrintSettings();                       //prints all settings
    int SizeToBin();                        //calculates the nr of integration bins based on the size of each bin
    int SizeToBin(double size);             //calculates the nr of integration bins based on size "size" of each bin
    int BinToSize();                    //calculates the size of integratio bin based on nr of integration bins
    int BinToSize(int n);                    //calculates the size of integratio bin based on nr of integration bins
    
    bool doSlidingWindow = false;           //says if sliding window variation should be performed
    bool doBinVariation = false;
    bool doAutoScale = false;               //if true, automatically scale gSF to Oslo data
    
    void SaveSettings();
    void ReadSettings();
    void setBgEne1(double ene[4]);
    void setBgEne2(double ene[4]);
    void ResetWidth();                      //resets the width calibration parameters for re-fitting
	double getEffCor(double ene, int level);			//returns an energy-dependend efficiency factor
    void readEffi();                        //reads the energy-depend efficiency factors
    TGraph* eGraph;                         //contains the data of the efficiency correction

};

#endif
