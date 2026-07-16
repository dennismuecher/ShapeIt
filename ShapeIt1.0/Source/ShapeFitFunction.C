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


//the class defining the peak fit function

#include "TF1.h"

class  ShapeFitFunction {
    
private:
    int multip;								//number of multiplet peaks (1=doublet, 2=triplet)
    bool fix_width;							//keep the width of all peaks the same during fit
    double bgRanges[4];						//ranges of the background (left and right region)
    double peakRanges[2];					//ranges of the peak
    double multi_gaus[4];					//gauss of multiplet peaks (supports up to 4 additional peaks)
    double gauss;							//gauss of main peak
    bool do_reject;							// flag to control if events outside the background ranges should be rejected
    
	bool Reject (double xx)	{					//returns true if point is not within bgRange or peakRange and do_reject is set to "true"
        
        if (xx >= bgRanges[0] && xx <= bgRanges[1] )
            return false;
        else if (xx >= bgRanges[2] && xx <= bgRanges[3] )
            return false;
        else if (xx >= peakRanges[0] && xx <= peakRanges[1] )
			  return false;
        else if ( do_reject )
            return true;
        else return false;
    }
    
public:
    
    //Constructor
    
    ShapeFitFunction (int multiplet_type, bool fix_multiplet_width = true ) {
		// multiplet_type: 0 = single peak, 1 = doublet only, 2 = triplet only, 3 = both doublet and triplet
		multip = multiplet_type;
		
        for (int j = 0; j < 4; j++) {
            bgRanges[j] = 0;
            multi_gaus[j] = 0;
        }
        gauss = 0;
        peakRanges[0] = 0; peakRanges[1] = 0;
        do_reject = true;
        fix_width = fix_multiplet_width;
    }
    
    void SetReject (bool ddo_reject) {do_reject = ddo_reject;}
    
    void SetBgRanges ( double bbgRanges[4]) {
        for (int i = 0; i < 4; i++)
            bgRanges[i] = bbgRanges[i];
    }
    
    void SetPeakRanges ( double ppeakRanges[2]) {
        for (int i = 0; i < 2; i++)
            peakRanges[i] = ppeakRanges[i];
    }
    
    
    double fitFunction_bg(double *x, double *par) {
        
        return ( par[0] * x[0] *x[0] ) + ( par[1] *x[0] ) + par[2];
        
    }
    
    double operator() (double *x, double *par) {
        // perform fit only if point is within background or peak region
        if ( this->Reject( x[0] ) ) {
            TF1::RejectPoint();
            return 0;
        }
        
        // Background: par[0-2]
        // Main peak: par[3] (amplitude), par[4] (position), par[5] (width)
        gauss = par[3]*exp(-0.5*TMath::Power(((x[0]-par[4])/par[5]),2));
        
        // Single peak (multip=0): no additional peaks
        if (multip == 0) {
            return this->fitFunction_bg(x, par) + gauss;
        }
        
        // Parameter layout depends on fix_width and which peaks are enabled:
        // multip=1: doublet only
        // multip=2: triplet only (skip doublet params, triplet uses 6-7 or 6-8)
        // multip=3: both doublet and triplet
        //
        // If fix_width=true:
        //   Doublet: par[6] (amp), par[7] (pos), uses par[5] for width
        //   Triplet: par[8] (amp), par[9] (pos), uses par[5] for width
        // If fix_width=false:
        //   Doublet: par[6] (amp), par[7] (pos), par[8] (width)
        //   Triplet: par[9] (amp), par[10] (pos), par[11] (width)
        
        // Doublet peak (if multip==1 or multip==3)
        if (multip == 1 || multip == 3) {
            if (fix_width)
                gauss += par[6]*exp(-0.5*TMath::Power(((x[0]-par[7])/par[5]),2));
            else
                gauss += par[6]*exp(-0.5*TMath::Power(((x[0]-par[7])/par[8]),2));
        }
        
        // Triplet peak (if multip==2 or multip==3)
        if (multip == 2 || multip == 3) {
            if (multip == 2) {
                // Triplet only - use params 6-7 (fixed width) or 6-8 (free width)
                if (fix_width)
                    gauss += par[6]*exp(-0.5*TMath::Power(((x[0]-par[7])/par[5]),2));
                else
                    gauss += par[6]*exp(-0.5*TMath::Power(((x[0]-par[7])/par[8]),2));
            } else {
                // Both peaks - triplet uses params 8-9 (fixed width) or 9-11 (free width)
                if (fix_width)
                    gauss += par[8]*exp(-0.5*TMath::Power(((x[0]-par[9])/par[5]),2));
                else
                    gauss += par[9]*exp(-0.5*TMath::Power(((x[0]-par[10])/par[11]),2));
            }
        }
        
        return this->fitFunction_bg(x, par) + gauss;
    }
};
