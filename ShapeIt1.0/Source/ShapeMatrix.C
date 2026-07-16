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

#include "../Include/ShapeMatrix.h"
#include "../Source/ShapeFitFunction.C"

ShapeMatrix::ShapeMatrix(ShapeSetting* setting) {
    sett = setting;
	ene0 = sett->exiEne[0];
	ene1 = sett->exiEne[1];
	esize = sett->exi_size[0];
    openMatrix();
    BrowseRootFile();
}

ShapeMatrix::~ShapeMatrix() {
    delete fit_result[0];
    delete fit_result[1];
    delete fit_func[0];
    delete fit_func[1];
    delete diag;
    delete diagEx;
    delete diagExSquare;
    delete diagExCube;
}


//cleaning up before new "ShapeIt" run (but not while running different iterations during bin size and sliding window variation)

void ShapeMatrix::Reset() {
    width1.clear();
    width2.clear();
    allgamma1.clear();
    allgamma2.clear();
}

TH2* ShapeMatrix::GetInputMatrix (std::string title) {
    inputMatrix->SetTitle(title.c_str());
    return inputMatrix;
}


void ShapeMatrix::openMatrix(){
    dataFile = new TFile(sett->dataFileName.c_str(),"read");
    if (dataFile->IsZombie()) {
          cout << "Error opening root file " <<sett->dataFileName.c_str() << endl;
          exit(-1);
       }
}

void ShapeMatrix::SetMatrix(int mNr) {

    dataFile->GetObject(matrixName[mNr-1].c_str(), inputMatrix);
    //inputMatrix->SetTitle(sett->dataFileName.c_str());
    sett->matrixName = matrixName[mNr-1];

}

//convert energy to bin for excitation energies
int ShapeMatrix::energyToBinY (double e) {
    TAxis *yaxis = diagEx->GetYaxis();
    return ( yaxis->FindBin(e) );
}


//convert energy to bin in projection histo for X
int ShapeMatrix::energyToBinX (double e) {
    TAxis *xaxis = diag->GetXaxis();
    return ( xaxis->FindBin(e) );
}

//returns the projection of diagonalized matrix for bin bbin
TH1D* ShapeMatrix::GetDiagEx(int bbin, std::string title) {
    if (bbin <= ybins) {
        char name[50];
        char title[50];
        snprintf(name,50,"diag_bin%d",bbin);
        snprintf(title,50,"projection %d - %d ",(int) ( (bbin -1) * esize + ene0 ), (int) ( bbin * esize + ene0 ) );
        
		//delet eprevious object
		delete gROOT->FindObject(name);
        
		//create projection of diagEx
        hist = new TH1D(name,title,xbins,eMin_diag, eMax_y);
        hist = diagEx->ProjectionX(name,bbin,bbin,"o");
        hist->GetXaxis()->SetTitle("E_{x} - E_{#gamma} (keV)");
        hist->GetYaxis()->SetTitle("counts");
        if (sett->mode == 2) {
            FitGauss(hist, bbin, 0);
            FitGauss(hist, bbin, 1);
        }
        return hist;
    }
    else
        return NULL;
}

//background integration and filling of net integral vector
void ShapeMatrix::IntegrateBg() {
    integral1Bg.clear();
    integral2Bg.clear();
    integral1Net.clear();
    integral2Net.clear();
    
    double a1, a2;
    //background integration for level 1
    int integBin[4];
    for (int i =0; i < 4; i++) {
        integBin[i] = energyToBinX(sett->bgEne[0][i]);
    }
    for (int i = 0; i < ybins; i++) {
        //left region
        a1 =diagEx->Integral(integBin[0],integBin[1], i+1, i+1);
        //right region
        a2 =diagEx->Integral(integBin[2],integBin[3], i+1, i+1);
        integral1Bg.push_back((a1+a2)/2);
    }
    
    //background integration for level 2
    for (int i =0; i < 4; i++) {
        integBin[i] = energyToBinX(sett->bgEne[1][i]);
    }
    for (int i = 0; i < ybins; i++) {
        //left region
        a1 =diagEx->Integral(integBin[0],integBin[1], i+1, i+1);
        //right region
        a2 =diagEx->Integral(integBin[2],integBin[3], i+1, i+1);
        integral2Bg.push_back((a1+a2)/2);
        if (sett->verbose)
            std::cout <<"Background Integral of bin " <<i+1 <<": Level 1 = " <<integral1Bg[i] << ", Level 2 = " <<integral2Bg[i] <<std::endl;
        
    }
}
void ShapeMatrix::Integrate() {
    integral1.clear();
    integral2.clear();
    int integBin[4];
    for (int i =0; i < 4; i++) {
        integBin[i] = energyToBinX(sett->levEne[i]);
    }
    for (int i = 0; i < ybins; i++) {
        integral1.push_back(diagEx->Integral(integBin[0],integBin[1], i+1, i+1) );
        integral2.push_back(diagEx->Integral(integBin[2],integBin[3], i+1, i+1) );
        if (sett->verbose)
            std::cout <<"Integral of bin " <<i+1 <<": Level 1 = " <<integral1[i] << ", Level 2 = " <<integral2[i] <<std::endl;
    }
}

void ShapeMatrix::IntegrateSquare() {
    integral1Square.clear();
    integral2Square.clear();
    int integBin[4];
    for (int i =0; i < 4; i++) {
        integBin[i] = energyToBinX(sett->levEne[i]);
    }
    for (int i = 0; i < ybins; i++) {
        integral1Square.push_back(diagExSquare->Integral(integBin[0],integBin[1], i+1, i+1) );
        integral2Square.push_back(diagExSquare->Integral(integBin[2],integBin[3], i+1, i+1) );
    }
}
void ShapeMatrix::IntegrateCube() {
    integral1Cube.clear();
    integral2Cube.clear();
    int integBin[4];
    for (int i =0; i < 4; i++) {
        integBin[i] = energyToBinX(sett->levEne[i]);
    }
    for (int i = 0; i < ybins; i++) {
        integral1Cube.push_back(diagExCube->Integral(integBin[0],integBin[1], i+1, i+1) );
        integral2Cube.push_back(diagExCube->Integral(integBin[2],integBin[3], i+1, i+1) );
    }
}

//determines the peak areas using autofit
void ShapeMatrix::FitIntegral(){
    //clearing all previous results
    fit_integral1.clear();
    fit_integral2.clear();
    fit_integral1Bg.clear();
    fit_integral2Bg.clear();
    fit_integral1Net.clear();
    fit_integral2Net.clear();
   
    //normalization of integral to bin width: the integration of the fit result is in keV whereas the matrix is integrated by bins
    double integral_norm = diagEx->GetXaxis()->GetBinWidth(1);
    
    //level 1
    for (int i = 0; i < ybins; i++) {
        //if (integral1[i] == 0) {
        if (integral1[i] < sett->minCounts) {
            fit_integral1.push_back(0);
            fit_integral1Bg.push_back(0);
            fit_integral1Net.push_back(0);
            width1.push_back(0);
            allgamma1.push_back(0);
        }
        else {
            //Gauss fit
            //NOTE: ProjectionX() allocates a new histogram every call. It used to be passed
            //straight into FitGauss() as a temporary with no variable ever holding the
            //pointer, so it (and everything FitGauss()'s three Fit(...,"+",...) calls attach
            //to it) was never freed -- once per bin, per level, per iteration.
            TH1D *projLevel1 = diagEx->ProjectionX("fit",i+1,i+1,"o");
            FitGauss(projLevel1, i, 0);
            
            //integration
            double integral_tot =fit_result[0]->Integral(sett->levEne[0], sett->levEne[1]) / integral_norm;
            
            fit_integral1.push_back(integral_tot);
            //integral of background, only: setting amplitude of Gauss to zero
            fit_result[0]->FixParameter(3,0);
            double integral_bg =fit_result[0]->Integral(sett->levEne[0], sett->levEne[1]) / integral_norm;
            fit_integral1Bg.push_back(integral_bg);
            fit_integral1Net.push_back(integral_tot - integral_bg);
            
            //store fit result and corresponding gamma ray energy
            if (fit_integral1Net[i] >= sett->minCounts) {
                width1.push_back(TMath::Abs(fit_result[0]->GetParameter(5)));
                allgamma1.push_back(integral1Square[i] /integral1Cube[i]);
            }
            
            if (sett->verbose)
                std::cout <<"Autofit results for level 1 in bin " << i <<"(tot, bg, net): " <<integral_tot <<" " <<integral_bg <<" " <<fit_integral1Net[i] <<std::endl;

            delete projLevel1;
        }
    }
    //level 2
    for (int i = 0; i < ybins; i++) {
        //if (integral2[i] == 0) {
        if (integral2[i] < sett->minCounts) {

            fit_integral2.push_back(0);
            fit_integral2Bg.push_back(0);
            fit_integral2Net.push_back(0);
            width2.push_back(0);
            allgamma2.push_back(0);
        }
        else {
            //Gauss fit
            TH1D *projLevel2 = diagEx->ProjectionX("fit",i+1,i+1,"o");
            FitGauss(projLevel2, i, 1);
            
            //integration
            double integral_tot =fit_result[1]->Integral(sett->levEne[2], sett->levEne[3]) / integral_norm;
            fit_integral2.push_back(integral_tot);
            //integral of background, only: setting amplitude of Gauss to zero
            fit_result[1]->FixParameter(3,0);
            double integral_bg =fit_result[1]->Integral(sett->levEne[2], sett->levEne[3]) / integral_norm;
            fit_integral2Bg.push_back(integral_bg);
            fit_integral2Net.push_back(integral_tot - integral_bg);
            
            if (fit_integral2Net[i] >= sett->minCounts) {
                //store fit result and corresponding gamma ray energy
                width2.push_back(TMath::Abs(fit_result[1]->GetParameter(5)));
                allgamma2.push_back(integral2Square[i] /integral2Cube[i]);
            }
                
            if (sett->verbose)
                std::cout <<"Autofit results for level 2 in bin " << i <<"(tot, bg, net): " <<integral_tot <<" " <<integral_bg <<" " <<fit_integral2Net[i] <<std::endl;

            delete projLevel2;
        }
        
    }
    
}

//performs a gauss fit to histo of bin "bin" for level 1 (level =0) or level 2 (level =1)
void ShapeMatrix::FitGauss(TH1D *histo, int bin, int level) {
    char name[50];
    // Determine multiplet type based on which peaks are enabled
    // multiplet_type: 0=single peak, 1=doublet (peak 2 only), 2=triplet (peak 3 only), 3=both doublet and triplet
    int multiplet_type = 0;
    
    // Check if doublet should be included (checkbox enabled and position is valid/non-zero)
    bool validDoublet = sett->doDoublet[level] && sett->doubletPeakPos[level] > 0;
    
    // Check if triplet should be included (checkbox enabled and position is valid/non-zero)
    bool validTriplet = sett->doTriplet[level] && sett->tripletPeakPos[level] > 0;
    
    // Encode both states: 1=doublet, 2=triplet, 3=both
    if (validDoublet && validTriplet) {
        multiplet_type = 3;  // both peaks
    } else if (validTriplet) {
        multiplet_type = 2;  // triplet only (no doublet)
    } else if (validDoublet) {
        multiplet_type = 1;  // doublet only
    }
    
    if (sett->verbose) {
        std::cout << "FitGauss: level=" << level << " bin=" << bin 
                  << " multiplet_type=" << multiplet_type 
                  << " doDoublet=" << sett->doDoublet[level]
                  << " doTriplet=" << sett->doTriplet[level]
                  << " peakPos=" << sett->peakPos[level]
                  << " doubletPeakPos=" << sett->doubletPeakPos[level]
                  << " tripletPeakPos=" << sett->tripletPeakPos[level] << std::endl;
    }
    
    bool fix_multiplet_width = ((multiplet_type == 1) || (multiplet_type == 3)) ? sett->fixDoubletWidth[level] : 
                               (multiplet_type == 2) ? sett->fixTripletWidth[level] : true;
    
    //background regions
    double bgRange[4];
    for (int i =0; i <4; i++)
        bgRange[i] = sett->bgEne[level][i];
    
    //peak region - SINGLE fit region for ALL peaks for this level
    double peakRange[2];
    peakRange[0] = sett->levEne[2*level];
    peakRange[1] = sett->levEne[2*level+1];
			
    //peak energy - use explicit peak positions from UI as starting values
    double p = sett->peakPos[level];
    double p_2 = sett->doubletPeakPos[level];
    double p_3 = sett->tripletPeakPos[level];
	

    //initial guess for sigma: depends on number of peaks to maintain reasonable width
    // Single peak: 1/6 of range, doublet: 1/12, triplet: 1/18
    double dp;
    if (multiplet_type == 0) {
        dp = (sett->levEne[2*level+1] - sett->levEne[2*level]) / 6;
    } else if (multiplet_type == 1) {
        dp = (sett->levEne[2*level+1] - sett->levEne[2*level]) / 12;
    } else {
        dp = (sett->levEne[2*level+1] - sett->levEne[2*level]) / 18;
    }

    //define fit function and set ranges
    //free the previous fit's functor and TF1 before overwriting fit_result[level]/fit_func[level] --
    //this used to just overwrite the pointers, leaking both every single call (i.e. every bin, every level, every iteration)
    delete fit_func[level];
	ShapeFitFunction *fitfunc = new ShapeFitFunction(multiplet_type, fix_multiplet_width);
    fit_func[level] = fitfunc;
    fitfunc->SetPeakRanges(peakRange);
    fitfunc->SetBgRanges(bgRange);
    
    //TF1 object for the fit
    snprintf(name,50,"fit_level%d_bin%d",level+1, bin);
    delete fit_result[level];
    
    // Calculate number of parameters based on multiplet type
    // Background: 3 params (0-2)
    // Main peak: 3 params (3-5: amplitude, position, width)
    // 
    // multiplet_type: 0=single, 1=doublet only, 2=triplet only, 3=both
    //
    // If fix_multiplet_width=true (widths equal):
    //   Doublet only (1): 2 params (6-7: amp, pos) - uses par[5] for width
    //   Triplet only (2): 2 params (6-7: amp, pos) - uses par[5] for width
    //   Both (3): 4 params (6-9: doublet amp/pos, triplet amp/pos) - uses par[5] for width
    //   Total: 6, 8, or 10 params
    //
    // If fix_multiplet_width=false (widths independent):
    //   Doublet only (1): 3 params (6-8: amp, pos, width)
    //   Triplet only (2): 3 params (6-8: amp, pos, width)
    //   Both (3): 6 params (6-11: doublet amp/pos/width, triplet amp/pos/width)
    //   Total: 6, 9, or 12 params
    int nParams = 6;  // bg (3) + main peak (3)
    if (multiplet_type == 1 || multiplet_type == 2) {
        // Doublet only or triplet only
        nParams += fix_multiplet_width ? 2 : 3;
    } else if (multiplet_type == 3) {
        // Both doublet and triplet
        nParams += fix_multiplet_width ? 4 : 6;
    }
    
    fit_result[level] = new TF1(name,fitfunc ,eMin_y,eMax_y, nParams );
	
    //get bin content at peak and bgRange[1] as starting value for fit amplitude
    double amplitude_init = histo->GetBinContent(energyToBinX(p));
    double amplitude_init_2 = histo->GetBinContent(energyToBinX(p_2));
    double amplitude_init_3 = histo->GetBinContent(energyToBinX(p_3));
    
	double amplitude_init_bg = histo->GetBinContent(energyToBinX(bgRange[1]));
    
    //set initial fit values
    fit_result[level]->SetParameter(1, 0);
    fit_result[level]->SetParameter(2, amplitude_init_bg);
	
    //do linear background
    fit_result[level]->FixParameter(0,0);
    
    //perform fit of background, first
	fit_result[level]->SetLineColor(kCyan-6);
	fit_result[level]->SetLineWidth(3);
	fit_result[level]->FixParameter(3, 0);
	if (multiplet_type == 1 || multiplet_type == 3) 
			fit_result[level]->FixParameter(6, 0);
	if (multiplet_type == 2 || multiplet_type == 3) {
		// Triplet amplitude: par[6] if triplet-only, par[8] if both (widths fixed), par[9] if both (widths free)
		int triplet_amp_idx = (multiplet_type == 2) ? 6 : (fix_multiplet_width ? 8 : 9);
		fit_result[level]->FixParameter(triplet_amp_idx, 0);
	}
	//set peakrange to zero
	double peakRange_temp[2]={0,0};
	fitfunc->SetPeakRanges(peakRange_temp);
	//fit
	histo->Fit(name,"RQ+","",bgRange[0],bgRange[3]);
	
	//fix background fit parameter
	fit_result[level]->FixParameter(1, fit_result[level]->GetParameter(1));
	fit_result[level]->FixParameter(2, fit_result[level]->GetParameter(2));

	//set peak range for peak fit
	fitfunc->SetPeakRanges(peakRange);
	
	fit_result[level]->SetParameter(3, amplitude_init );
    fit_result[level]->SetParameter(4, p);
    fit_result[level]->SetParameter(5, dp);
    
	if (multiplet_type == 1 || multiplet_type == 3) {
		// Doublet peak
		fit_result[level]->SetParameter(6, amplitude_init_2);
		fit_result[level]->SetParameter(7, p_2);
		if (!fix_multiplet_width) {
			fit_result[level]->SetParameter(8, dp);
		}
	}
	
	if (multiplet_type == 2 || multiplet_type == 3) {
		// Triplet peak
		if (multiplet_type == 2) {
			// Triplet only - params 6-7 (fixed) or 6-8 (free)
			fit_result[level]->SetParameter(6, amplitude_init_3);
			fit_result[level]->SetParameter(7, p_3);
			if (!fix_multiplet_width) {
				fit_result[level]->SetParameter(8, dp);
			}
		} else {
			// Both peaks - params 8-9 (fixed) or 9-11 (free)
			if (fix_multiplet_width) {
				fit_result[level]->SetParameter(8, amplitude_init_3);
				fit_result[level]->SetParameter(9, p_3);
			} else {
				fit_result[level]->SetParameter(9, amplitude_init_3);
				fit_result[level]->SetParameter(10, p_3);
				fit_result[level]->SetParameter(11, dp);
			}
		}
	}
	
    //set fit boundaries - all peaks search within the same fit region
    fit_result[level]->SetParLimits(3,0.01*amplitude_init, 100*amplitude_init);
    fit_result[level]->SetParLimits(4, sett->levEne[2*level], sett->levEne[2*level+1]);
    fit_result[level]->SetParLimits(5, 0.5*dp, 2*dp);
	
	if (multiplet_type == 1 || multiplet_type == 3) {
		// Doublet peak
		fit_result[level]->SetParLimits(6,0.01*amplitude_init_2, 100*amplitude_init_2);
		fit_result[level]->SetParLimits(7, sett->levEne[2*level], sett->levEne[2*level+1]);
		if (!fix_multiplet_width) {
			fit_result[level]->SetParLimits(8, 0.5*dp, 2*dp);
		}
	}
	
	if (multiplet_type == 2 || multiplet_type == 3) {
		// Triplet peak
		if (multiplet_type == 2) {
			// Triplet only - params 6-7 (fixed) or 6-8 (free)
			fit_result[level]->SetParLimits(6, 0.01*amplitude_init_3, 100*amplitude_init_3);
			fit_result[level]->SetParLimits(7, sett->levEne[2*level], sett->levEne[2*level+1]);
			if (!fix_multiplet_width) {
				fit_result[level]->SetParLimits(8, 0.5*dp, 2*dp);
			}
		} else {
			// Both peaks - params 8-9 (fixed) or 9-11 (free)
			if (fix_multiplet_width) {
				fit_result[level]->SetParLimits(8, 0.01*amplitude_init_3, 100*amplitude_init_3);
				fit_result[level]->SetParLimits(9, sett->levEne[2*level], sett->levEne[2*level+1]);
			} else {
				fit_result[level]->SetParLimits(9, 0.01*amplitude_init_3, 100*amplitude_init_3);
				fit_result[level]->SetParLimits(10, sett->levEne[2*level], sett->levEne[2*level+1]);
				fit_result[level]->SetParLimits(11, 0.5*dp, 2*dp);
			}
		}
	}

    //if doFitWidth is set, fix width according to calibration
    if (sett->doWidthCal)
    {
        double eGamma = 0;
        if (level == 0 && integral1Cube[bin] > 0)
            eGamma = integral1Square[bin] /integral1Cube[bin];
        if (level == 1 && integral2Cube[bin] > 0)
            eGamma = integral2Square[bin] /integral2Cube[bin];
        double widthFix = sett->widthCal[level][1] * eGamma + sett->widthCal[level][0];
        if (!(eGamma == 0)) {
            fit_result[level]->FixParameter(5, widthFix);
            //fit_result[level]->SetParameter(5, widthFix);
            //fit_result[level]->SetParLimits(5, widthFix-5, widthFix +5);
            
        }
    }
    
    //if fixPeakPos is set, fix peak position according to user setting
    if (sett->fixPeakPos[level])
    {
        fit_result[level]->FixParameter(4, sett->peakPos[level]);
    }
    
    //if fixDoubletPeakPos is set, fix doublet peak position according to user setting
    if ((multiplet_type == 1 || multiplet_type == 3) && sett->fixDoubletPeakPos[level])
    {
        fit_result[level]->FixParameter(7, sett->doubletPeakPos[level]);
    }
    
    //if fixTripletPeakPos is set, fix triplet peak position according to user setting
    if ((multiplet_type == 2 || multiplet_type == 3) && sett->fixTripletPeakPos[level])
    {
        // Triplet position: par[7] if triplet-only, par[9] if both (widths fixed), par[10] if both (widths free)
        int triplet_pos_idx = (multiplet_type == 2) ? 7 : (fix_multiplet_width ? 9 : 10);
        fit_result[level]->FixParameter(triplet_pos_idx, sett->tripletPeakPos[level]);
    }
    
    //perform fit
    snprintf(name,50,"fit_level%d_bin%d",level+1, bin);
    histo->Fit(name,"RQ+","",bgRange[0],bgRange[3]);
	
	//fix all fit parameters and refit function without rejection for a nicer display
	for (int i = 0; i < nParams; i++) 
		fit_result[level]->FixParameter(i, fit_result[level]->GetParameter(i));
	fitfunc->SetReject(false);
    fit_result[level]->SetLineColor(6);
	if (level == 0)
        histo->Fit(name,"RQ+","",bgRange[0],bgRange[3]);
    else
        histo->Fit(name,"RQ+","",bgRange[0],bgRange[3]);
	

	//set amplitude of doublet and triplet peaks to zero as we only want the peak content of the peak of interest
	if (multiplet_type == 1 || multiplet_type == 3) 
		fit_result[level]->FixParameter(6,0);
	if (multiplet_type == 2 || multiplet_type == 3) {
		// Triplet amplitude: par[6] if triplet-only, par[8] if both (widths fixed), par[9] if both (widths free)
		int triplet_amp_idx = (multiplet_type == 2) ? 6 : (fix_multiplet_width ? 8 : 9);
		fit_result[level]->FixParameter(triplet_amp_idx,0);
	}
}


//creates the matrices for shape method for excitation energies between ene0 and ene1, with bin size esize for the excitation energies
void ShapeMatrix::Diag(){
	
    double MeV = sett->MeV;
    XNum  = inputMatrix->GetNbinsX();
    YNum  = inputMatrix->GetNbinsY();
 
    if (sett->verbose) {
        std::cout <<"\nNumber of x bins in matrix: "<< XNum <<std::endl;
        std::cout <<"\nNumber of y bins in matrix: "<< YNum <<std::endl;
    }

	//number of bins used for x-axis (Ex - Eg) in diagonalized matrix; I choose the number of bins equal to the number of bins in the gamma ray axis of the input matrix
	//xbins=XNum;
	
    eMax_x  =  ( inputMatrix->GetXaxis()->GetBinCenter(XNum) * MeV) + (inputMatrix->GetXaxis()->GetBinWidth(XNum) * MeV /2 ) ;
    eMax_y  =  ( inputMatrix->GetYaxis()->GetBinCenter(YNum) * MeV) + (inputMatrix->GetYaxis()->GetBinWidth(YNum) * MeV /2 ) ;
    eMin_x = inputMatrix->GetXaxis()->GetXmin();
    eMin_y = inputMatrix->GetYaxis()->GetXmin();
    
    //minimum value of Excitation energy - gamma ray energy
    eMin_diag = eMin_y - eMin_x;
    //creating some extra space at low energies to allow proper fits at low energies 
	eMin_diag = eMin_diag -100;
    
    if (sett->verbose) {
        std::cout <<"\nMaximum gamma ray energy detected in input matrix: "<< eMax_x <<std::endl;
        std::cout <<"\nMaximum excitation energy detected in input matrix: "<< eMax_y <<std::endl;
        std::cout <<"\nMinimum gamma ray energy detected in input matrix: "<< eMin_x <<std::endl;
        std::cout <<"\nMinimum excitation energy detected in input matrix: "<< eMin_y <<std::endl;
    }
   
	//calculate number of xbins
	//xbins = (int) (eMax_y - eMin_diag) / inputMatrix->GetXaxis()->GetBinWidth(50);
    xbins = (int) ( (eMax_y - eMin_diag) * XNum / (eMax_x - eMin_x));
    //calculate number of ybins
	ybins = ( ene1 - ene0 ) / esize;
	//increase number of ybins by one if we haven't reached the chosen maximum excitation energy, yet
	if (  (int)(ene1 - ene0) % (int)esize != 0)
		ybins++;
 	//update highest energy
    //ene1_diag is the value of the highest energy used to create the actual matricies; this value is equal to or larger than sett->ene1 to have a full-sized last bin; the value in sett->ene1 remains unchanged
	double ene1_diag = esize * ybins + ene0;
	
    if (sett->verbose) {
        std::cout <<"\nx binnings for diag matrix:"<< xbins <<" " <<eMin_diag << " " <<eMax_y <<std::endl;
        std::cout <<"\ny binnings for diag matrix:"<< ybins <<" " <<ene0 << " " <<ene1_diag <<std::endl;
	}
	
    //cleaning up before generating matrices.
    //NOTE: this used to be gROOT->FindObject("diag") etc. followed by delete -- but new
    //histograms register to whatever TDirectory is "current" at construction time, which
    //while a matrix TFile is open is usually that file, not gROOT. So the lookup would
    //silently find nothing and the delete was a no-op, leaking a full set of these
    //histograms on every call (i.e. every Monte Carlo iteration, since this runs once per
    //ShapeGSF construction). Deleting the actual member pointers is robust regardless of
    //the current directory.
    delete diag;
    delete diagEx;
    delete diagExSquare;
    delete diagExCube;
    
	//creating histograms;
    diag  = new TH1F("diag","Diagonal Projection",xbins, eMin_diag, eMax_y);
    diag->GetXaxis()->SetTitle("E_{x} - E_{#gamma} (keV)");
    diag->GetYaxis()->SetTitle("counts");
    diagEx= new TH2F("diag_ex","Diagonal Projection vs excitation energy",xbins,eMin_diag, eMax_y, ybins, ene0, ene1_diag);
    
    
    diagExSquare = new TH2F("diag_ex_square","Diagonal Projection vs excitation energy squared", xbins,eMin_diag, eMax_y, ybins, ene0, ene1_diag);
	diagExCube = new TH2F("diag_ex_cube","Diagonal Projection vs excitation energy cubed", xbins,eMin_diag, eMax_y, ybins, ene0, ene1_diag);
    
    double X; // corresponds to gamma ray energy
    double Y; // corresponds to excitation energy
    
     //Loop over input matrix
    double diag_ene;
    
    TAxis *xaxis = inputMatrix->GetXaxis();
    TAxis *yaxis = inputMatrix->GetYaxis();
    
    //fill the diag matrices and histos by looping over inputMatrix
    for (int i = 1; i < XNum + 1; i++) {
         X = xaxis->GetBinCenter(i);
        for (int j = 1; j < YNum + 1; j++) {
            Y = yaxis->GetBinCenter(j);
            diag_ene = Y-X;
            diag->Fill( diag_ene, inputMatrix->GetBinContent(i,j)*MeV);
            diagEx->Fill( diag_ene,Y, inputMatrix->GetBinContent(i,j)*MeV);
            if (X > 0) {
                diagExSquare->Fill( diag_ene,Y, inputMatrix->GetBinContent(i,j) * MeV / TMath::Power( X, 2 ) );
                diagExCube->Fill( diag_ene,Y, inputMatrix->GetBinContent(i,j) * MeV / TMath::Power( X, 3 ) );
            }
        }
    }
    
    //fill vectors with integration results
    Integrate();
    IntegrateBg();
    IntegrateSquare();
    IntegrateCube();
}



TGraph* ShapeMatrix::getFitWidthGraph(int level) {
    
    const std::vector<double>& width = (level == 0) ? width1 : width2;
    const std::vector<double>& allgamma = (level == 0) ? allgamma1 : allgamma2;
    
    double peakEne = ( sett->levEne[2*level] + sett->levEne[2*level+1] ) /2;

    std::vector<double> x;
    std::vector<double> y;
    
    for (int i = 0; i < width.size(); i++) {
        if (allgamma[i] >= sett->exiEne[0] - peakEne && allgamma[i] <= sett->exiEne[1] - peakEne && width[i] > 0) {
           x.push_back(allgamma[i]);
            y.push_back(width[i]);
        }
    }
    
    TGraph *T = new TGraph(x.size(), x.data(), y.data());
    
    if (level == 0)
        T->SetMarkerColor(kRed);
    else
        T->SetMarkerColor(kBlue);
    
    T->SetMarkerStyle(21);
    
    return T;
}

//loop through all objects in root file and store all TH2D object names in matrixName
void ShapeMatrix::BrowseRootFile()
{
   
    TIter keyList(dataFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)keyList())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH2")) continue;
        TH2 *h = (TH2*)key->ReadObj();
        matrixName.push_back(h->GetName());
    
    }
  
}
