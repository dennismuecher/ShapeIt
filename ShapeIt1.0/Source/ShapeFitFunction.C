//the class defining the fit function

#include "TF1.h"

class  ShapeFitFunction {
    
private:
    int multip;								//number of multiplet peaks
    bool fix_width;							//keep the width of all peaks the same during fit
    double bgRanges[4];						//ranges of the background (left and right region)
    double peakRanges[2];					//ranges of the peak
    double multi_gaus[4];					//guass of multiplet peaks
    double gauss;							//gauss of main peak
    bool do_reject;							// flag to control if events outside the background ranges should be rejected
    bool Reject (double xx)	{					//returns true if point is not within bgRange or peakRange
        
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
    
    ShapeFitFunction ( ) {
        multip = 0;
        for (int j = 0; j < 4; j++) {
            bgRanges[j] = 0;
            multi_gaus[j] = 0;
        }
        gauss = 0;
        peakRanges[0] = 0; peakRanges[1] = 0;
        do_reject = false;
        fix_width = false;
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
        //define one gaussian for each multiplet peak
        //the first gaussian has parameters par[3], par[4], par[5], the second peak has par[6], par[7]
        
        gauss = par[3]*exp(-0.5*TMath::Power(((x[0]-par[4])/par[5]),2));
        for (int j = 0; j < multip; j++) {
            if (fix_width)
                multi_gaus[j] = par[3*j + 6]*exp(-0.5*TMath::Power((2.35*(x[0]-par[3*j + 7])/par[5]),2));
            else
                multi_gaus[j] = par[3*j + 6]*exp(-0.5*TMath::Power((2.35*(x[0]-par[3*j + 7])/par[3*j + 8]),2));
            gauss = gauss + multi_gaus[j];
            
        }
        return  this->fitFunction_bg(x, par) + gauss;
        
    }
};
