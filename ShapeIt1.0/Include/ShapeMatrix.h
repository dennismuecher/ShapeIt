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

#ifndef SHAPEMATRIX_H
#define SHAPEMATRIX_H

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "ShapeSetting.h"

#include <TMath.h>
#include <TObject.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TKey.h>

class ShapeMatrix : public TObject {
	
private:
    ShapeSetting * sett;
    TFile* dataFile;                           //the data file containing the matrix
    TH2* inputMatrix;                          //the input matrix in dataFile
    TH1F* diag;	                                //the matrix of (excitation energy - gamma-ray energy)
    TH2F* diagEx;                               ///diag matrix vs excitation energy
    
    TH2F* diagExSquare;                           //diag matrix with factor E_g^2 for each entry vs excitation energy; used for calculating the average E_gamma for each bin
    TH2F* diagExCube;                           //diag matrix with factor E_g^3 for each entry vs excitation energy
    void diagMatrix();	                       //creates the diagonalized matrices
    void openMatrix();                         //opens the matrix defined in settings
    void openMatrix(std::string fname);              //open root file with name fname
    std::vector <std::string> matrixName;    //names of all matrices found in dataFile
   
    int xbins, ybins;							//number of x and y bins for diagonalized matrices 
	double ene0, ene1;							//lowest and highest excitation energies to be used in shape method
	double esize;								//size of excitation energy bin used in shape method
	int XNum, YNum;
    std::vector <double> width1;                //peak width for level 1 from all autofits
    std::vector <double> width2;                //peak width for leve2 1 from all autofits
    std::vector <double> allgamma1;                //corresponding excitation energy for width1 vector
    std::vector <double> allgamma2;                //corresponding excitation energy for width2 vector
    
public:
	ShapeMatrix(ShapeSetting* setting);
    
	double GetEne0() {return ene0;}
	double GetEne1() {return ene1;}
	double GetESize() {return esize;}
	
	void SetEne0(double nene0) {ene0 = nene0;}
	void SetEne1(double nene1) {ene1 = nene1;}
	void SetESize(double nesize) {esize = nesize;}
	
	int GetYBins() {return ybins;}
    double eMin_x, eMax_x, eMin_y, eMax_y;
    double eMin_diag;                           //minimum value of Excitation Energy and Gamma ray energy
    
    void Diag();
    void BrowseRootFile();
    std::vector <std::string> GetMatrixName() {return matrixName;}
    void SetMatrix(int mNr);                  //assigns inputMatrix to mNr matrix found in dataFile
    
	TH2* GetInputMatrix (std::string title);
    TH1F* GetDiag(std::string title) {diag->SetTitle(title.c_str());  return diag;}
    TH2F* GetDiagEx(std::string title) {diagEx->SetTitle(title.c_str()); return diagEx;}
    TH2F* GetDiagExCube(std::string title) {diagExCube->SetTitle(title.c_str()); return diagExCube;}
    TH1D* GetDiagEx(int bbin, std::string title);                  //returns projection of DiagEx for bin bbin, starting at bin 1

    int energyToBinX (double e);
    int energyToBinY (double e);
    
	void Reset();
    void Integrate();
    void IntegrateBg();
    void IntegrateSquare();
    void IntegrateCube();
    void FitIntegral();                             //calculates peak areas of level 1 and level 2 for all bins
                           
    std::vector <double> integral1;           //integral values of level1
    std::vector <double> integral2;           //integral values of level2
    std::vector <double> integral1Bg;           //integral values of background level1
    std::vector <double> integral2Bg;           //integral values of background level2
    std::vector <double> integral1Net;           //integral values of level1, background subtracted
    std::vector <double> integral2Net;           //integral values of level2, background subtracted
    
    std::vector <double> fit_integral1;           //integral values of level1
    std::vector <double> fit_integral2;           //integral values of level2
    std::vector <double> fit_integral1Bg;           //integral values of background level1
    std::vector <double> fit_integral2Bg;           //integral values of background level2
    std::vector <double> fit_integral1Net;           //integral values of level1, background subtracted
    std::vector <double> fit_integral2Net;           //integral values of level2, background subtracted
    
    std::vector <double> integral1Square;           //integral values of level1 in diagExSquare
    std::vector <double> integral2Square;           //integral values of level2 in diagExSquare
    
    std::vector <double> integral1Cube;           //integral values of level1 in diagExCube
    std::vector <double> integral2Cube;           //integral values of level2 in diagExCube
    void FitGauss(TH1D *histo, int bin , int level);             //performs a gauss fit for bin and level
    TF1* fit_result[2];                                 //stores the autofit results
    TGraph* getFitWidthGraph(int level);                           //returns a graph containing the width results of the autofit for level 1
    //TGraph* getFitWidthGraph2;                      //returns a graph containing the width results of the autofit for level 2
};
#endif
