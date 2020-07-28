#include "../Include/ShapeGSF.h"
#include <algorithm>

//constructor using setting file and a matrix
ShapeGSF::ShapeGSF(ShapeSetting* t_sett, ShapeMatrix* t_matrix):m_sett(t_sett), m_matrix(t_matrix)
{
    nOfLevels = 2;
    levGraph[0] = new TGraph();
    levGraph[1] = new TGraph();
    mergeGraph = new TGraph();
}

//constructor using setting file; in this case use literature values to fill levGraph[0]
ShapeGSF::ShapeGSF(ShapeSetting* t_sett, ShapeMatrix* t_matrix):m_sett(t_sett), m_matrix(t_matrix)
{
    nOfLevels = 1;
    levGraph[0] = new TGraph();
    mergeGraph = new TGraph();
}




//scales the gSF pair at position i to the one at position i-1
void ShapeGSF::SewNext(Int_t i) {
    if (i < levGraph[0]->GetN()-1) {
        double scale = levGraph[1]->GetPointY(i) / levGraph[0]
    }
        
}


void ShapeGSF::DoInterpol() {
    if (sett->verbose)
        std::cout <<"\nINTERPOLATION OF gSF DATA... " <<endl;
    
    //slope of the last two points
    double slope = 0;
    
    //interpolate from interEne to lower energies
      for (int i = gSF_matrix->energyToBinY(sett->sewingEne); i > (-1) ; i--) {
        //do nothing for the first loop
        if ( slope != 0 ) {
            
            //this is the interpolated value for gSF
            double y_norm = ( gSF[i].egamma1 -gSF[i+1].egamma1 ) * slope  + gSF[i+1].value1;
            if (gSF[i].value1 <= 0) {
                //in case the count rates are zero (or neagtive) the gSF value could still be defined through interpolation, but its error bar would not be defined anymore. One also cannot interpolate gSF[i].value2 anymore, so the entire bin, and all following bins, are ignored
                //erase all elements in gSF from current bin onwards
                binRange[0] = i+2; //will not include the current bin
                if (sett->verbose > 1)
                    std::cout <<"Setting lowest active bin to "<< i+2 <<std::endl;
                break;
            }
            else if (gSF[i].value2 <= 0) {
                //in this case, gSF[i].value1 could still be useful, but interpolation beyond this bin will not be useful, so all earlier bins will be ignored
                gSF[i].value2 = 0; gSF[i].egamma2 = 0;
                gSF[i].dvalue1 = gSF[i].dvalue1 * y_norm / gSF[i].value1;
                gSF[i].value1 = y_norm;
                //binRange[0] = i+1;//will include the current bin and terminate the interpolation
                binRange[0] = i+2;//will NOT include the current bin and terminate the interpolation
                if (sett->verbose > 1)
                    std::cout <<"Setting lowest active bin to "<< i+1 <<std::endl;
                break;
            }
            
            gSF[i].value2 = gSF[i].value2 * y_norm / gSF[i].value1;
            gSF[i].dvalue2 = gSF[i].dvalue2 * y_norm / gSF[i].value1;
            
            gSF[i].dvalue1 = gSF[i].dvalue1 * y_norm / gSF[i].value1;
            gSF[i].value1 = y_norm;
            
            if (sett->verbose > 1)
                std::cout <<"\nValues of gSF after interpolation: "<< gSF[i].value1 <<" "<< gSF[i].value2 << endl;
        }
        
        //calculate slope for next iteration
      
        slope = ( gSF[i].value1-gSF[i].value2 ) / ( gSF[i].egamma1 - gSF[i].egamma2 );
        
        if (sett->verbose > 1)
            std::cout << "\nSlope for interpolation:" <<slope<<endl;
    }
    
    //interpolate from interEne to higher energies
    slope = 0;
      for (int i = gSF_matrix->energyToBinY(sett->sewingEne); i < gSF_matrix->GetYBins(); i++) {
        //do nothing for the first loop
        if ( slope != 0 ) {
            
            //this is the interpolated value for gSF
            double y_norm = ( gSF[i].egamma2 -gSF[i-1].egamma2 ) * slope  + gSF[i-1].value2;
            
            if (gSF[i].value2 <= 0) {
                //in case the count rates are zero (or neagtive) the gSF value could still be defined through interpolation, but its error bar would not be defined anymore. One also cannot interpolate gSF[i].value1 anymore, so the entire bin, and all following bins, are ignored
                binRange[1] = i; //will not include the current bin
                if (sett->verbose > 1)
                    std::cout <<"Setting number of active bins to "<< i <<std::endl;
                break;
            }
            else if (gSF[i].value1 <= 0) {
                //in this case, gSF[i].value2 could still be useful, but interpolation beyond this bin will not be useful, so all following bins will be ignored
                gSF[i].value1 = 0; gSF[i].egamma1 = 0;
                gSF[i].dvalue2 = gSF[i].dvalue2 * y_norm / gSF[i].value2;
                gSF[i].value2 = y_norm;
                //binRange[1] = i+1; //will include the current bin
                binRange[1] = i; //will NOT include the current bin
                if (sett->verbose > 1)
                    std::cout <<"Setting number of active bins to "<< i+1 <<std::endl;
                break;
            }
            else {
                //both gSF values are > 0
                gSF[i].value1 = gSF[i].value1 * y_norm / gSF[i].value2;
                gSF[i].dvalue1 = gSF[i].dvalue1 * y_norm / gSF[i].value2;
                gSF[i].dvalue2 = gSF[i].dvalue2 * y_norm / gSF[i].value2;
                gSF[i].value2 = y_norm;
                
            }
            
            if (sett->verbose > 1)
                std::cout << "\nSlope for interpolation:" <<slope<<endl;
            if (sett->verbose > 1)
                std::cout <<"\nValues of gSF after interpolation: "<< gSF[i].value1 <<" "<< gSF[i].value2 << endl;
        }
        //calculate slope for next iteration
        slope = ( gSF[i].value1-gSF[i].value2 ) / ( gSF[i].egamma1 - gSF[i].egamma2 );

    }
}

//calculates the ratio of Peak to background for bin i

double ShapeGSF::getBgRatio(int bin, int level) {
    
    if (!sett->doBackground)
        return 1;
    
    if (level == 1 ) {
        
        double peak = gSF_matrix->integral1[bin];
        double backgr = gSF_matrix->fit_integral1Bg[bin];
        if (peak > 0)
            return (peak - backgr) / peak;
        else
            return 0;
    }
    
    if (level == 2 ) {
        
        double peak = gSF_matrix->integral2[bin];
        double backgr =gSF_matrix->fit_integral2Bg[bin];
        if (peak > 0)
            return (peak - backgr) / peak;
        else
            return 0;
    }
    return 0;
}


void ShapeGSF::FillgSF() {
    
    gSF_matrix->Integrate();
    gSF_matrix->IntegrateBg();
    gSF_matrix->IntegrateSquare();
    gSF_matrix->IntegrateCube();
    
    //if autofit or background subtraction is set, perform autofit using Gauss
    if (sett->mode == 2 || sett->doBackground)
        gSF_matrix->FitIntegral();
    
    double elevel1 = ( sett->levEne[0] + sett->levEne[1] ) / 2;
    double elevel2 = ( sett->levEne[2] + sett->levEne[3] ) / 2;
    double egamma1, egamma2;
    
    for (int i = 0; i < gSF_matrix->integral1Cube.size(); i++ ) {
        
        double egamma = 0;
        double value = 0;
        double dvalue = 0;
        if (getBgRatio(i, 1) * gSF_matrix->integral1[i] < sett->minCounts || gSF_matrix->integral1[i] == 0 ) {
                
            egamma = gSF_matrix->GetEne0() + ( i + 0.5) * gSF_matrix->GetESize()  - elevel1;
            value = 0;
            dvalue = 0;
        }
        else {
			
            //calculate E_gamma
            //this is the better way of calculating E_gamma by calculating the average E_g, weighted by gSF, i.e. egamma_avg = SUM (E_g*gSF) / SUM(gSF). E_g*gSF is contained in the gSFSquare matrix
            egamma = gSF_matrix->integral1Square[i] / gSF_matrix->integral1Cube[i];
            //this calculates the average E_gamma just using the middle of each excitation bin; not the preferred way as it doesn't take into account of feeding as a function iof excitation energy
            //gSF_t.egamma1 = ( gSF_matrix->ybins[i] + gSF_matrix->ybins[i+1] ) / 2 - elevel1;
            
            //calculate gamma ray strength
            value =  getBgRatio(i, 1) * gSF_matrix->integral1Cube[i];
            dvalue = TMath::Power(1./gSF_matrix->integral1[i], 0.5) * value;
            
            //recalculate gSF in case autofit is activated
            if (sett->mode == 2 && sett->doBackground) {
                value = gSF_matrix->fit_integral1Net[i] * gSF_matrix->integral1Cube[i]  / gSF_matrix->integral1[i];
                
                if (gSF_matrix->fit_integral1Net[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    value = 0;
                }
            }
            else if( sett->mode == 2 && !sett->doBackground) {
                value = gSF_matrix->fit_integral1[i] * gSF_matrix->integral1Cube[i]  / gSF_matrix->integral1[i];
                if (gSF_matrix->fit_integral1[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    value = 0;
                }
            }
        }
        
        levGraph[0]->SetPoint(i, egamma, value);
        levGraph[0]->SetPointErrors(i, 0, dvalue);
        
        if (getBgRatio(i, 2) * gSF_matrix->integral2[i] < sett->minCounts || gSF_matrix->integral2[i] == 0 ) {

            egamma = gSF_matrix->GetEne0() + ( i + 0.5) * gSF_matrix->GetESize() - elevel2;
            value = 0;
            dvalue = 0;
        }
        else {
            
            //calculate E_gamma
            //this is the better way of calculating E_gamma by calculating the average E_g, weighted by gSF, i.e. egamma_avg = SUM (E_g*gSF) / SUM(gSF). E_g*gSF is contained in the gSFSquare matrix
            egamma = gSF_matrix->integral2Square[i] / gSF_matrix->integral2Cube[i];

            //this calculates the average E_gamma just using the middle of each excitation bin; not the preferred way as it doesn't take into account of feeding as a function iof excitation energy
            //gSF_t.egamma1 = ( gSF_matrix->ybins[i] + gSF_matrix->ybins[i+1] ) / 2 - elevel1;
            
            //calculate gamma ray strength
            
            value = getBgRatio(i, 2) * gSF_matrix->integral2Cube[i] * sett->getEffCor(egamma);
            dvalue = TMath::Power(1./gSF_matrix->integral2[i], 0.5) * value;
            
            //recalculate gSF in case autofit is activated
            
            if (sett->mode == 2 && sett->doBackground) {
                
                //gSF_t.value2 = sett->eff_corr * gSF_matrix->fit_integral2Net[i] * gSF_matrix->integral2Cube[i]  / gSF_matrix->integral2[i];
                value = sett->getEffCor(egamma) * gSF_matrix->fit_integral2Net[i] * gSF_matrix->integral2Cube[i] / gSF_matrix->integral2[i];
                if (gSF_matrix->fit_integral2Net[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    value = 0;
                }
            }
            else if (sett->mode == 2 && !sett->doBackground) {
                value = sett->getEffCor(egamma) * gSF_matrix->fit_integral2[i] * gSF_matrix->integral2Cube[i]  / gSF_matrix->integral2[i];
                if (gSF_matrix->fit_integral2[i] < 0) {
                    std::cout <<"Fit result smaller than zero, skipping this event" <<std::endl;
                    value= 0;
                }
            }
        }
        
        levGraph[1]->SetPoint(i, egamma, value, dvalue);
        levGraph[1]->SetPointErrors(i, 0, dvalue);

        if (sett->verbose) {
            std::cout <<"Bin: " <<i+1<<std::endl;
            std::cout <<"energies: " <<levGraph[0]->GetPointX(i) << " " << levGraph[1]->GetPointX(i) <<std::endl;
            std::cout <<"gSF1: " <<levGraph[0]->GetPointY(i) << "+- " << levGraph[0]->GetPointEY(i) <<std::endl;
            std::cout <<"gSF2: " <<levGraph[1]->GetPointY(i) << "+- " << levGraph[1]->GetPointEY(i) <<std::endl;

        }
    }
}

void ShapeGSF::GetMergeGraph() {
    for (int i = 0; i < levGraph[0]->GetN(); i++) {
        for (int l = 0 ; l < nOfLevels; l++) {
            mergeGraph->SetPoint(i, levGraph[l]->GetPointX(i), levGraph[l]->GetPointY(i) );
            mergeGraph->SetPointErrors(i, levGraph[l]->GetPointEX(i), levGraph[l]->GetPointEY(i) );
        }
    }
}

//prints the gSF values, sorted for egamma
void ShapeGSF::gSF_Print() {
    
    std::cout << "\n\nResults for gamma ray strength function: " <<std::endl;
    mergeGraph()->Sort();
    for (int i = 0; i < mergeGraph->GetN(); i++ )
        std::cout << mergeGraph->GetPointX(i) <<"     "<< sett->gSF_norm * ( mergeGraph->GetPointY(i) + mergeGraph->GetPointEY(i) ) <<"     " << sett->gSF_norm *  (mergeGraph->GetPointY(i) - mergeGraph->GetPointEY(i) ) << std::endl;
}

//transforms all gSF values via B*exp(alpha E_gamma)
void ShapeGSF::Transform(double B_t, double alpha_t) {
    
    //gSF was previously transformed via B and Alpha, so only transform according to the change in B_t and Alpha_t
    for (int i = 0; i < gSF_sort.size() ; i++ ) {
        gSF_sort[i].value =  B_t / B * TMath::Exp(( alpha_t - alpha) * gSF_sort[i].egamma / 1000.) *gSF_sort[i].value;
        gSF_sort[i].dvalue = B_t / B * TMath::Exp(( alpha_t - alpha) * gSF_sort[i].egamma / 1000.) *gSF_sort[i].dvalue;
    }
    B = B_t;
    alpha = alpha_t;
}

//provides a plot of gSF read from a file
TGraphErrors* ShapeGSF::plotLit() {
    
    int nOfPoints = gSF_sort.size();
    
    double x[nOfPoints], y[nOfPoints];
    double dx[nOfPoints], dy[nOfPoints];

    for (int i = 0; i < nOfPoints ; i++ ) {
        x[i] = gSF_sort[i].egamma;
        y[i] = gSF_sort[i].value;
        dx[i] = 1;
        dy[i] = gSF_sort[i].dvalue;
    }
    TGraphErrors *lit_data = new TGraphErrors(nOfPoints,x,y,dx,dy);
    lit_data->SetFillColor(4);
    lit_data->SetFillStyle(3010);
    return lit_data;
}


//creats TGraphError using the sorted data and drawing option for plotitng a band
TMultiGraph* ShapeGSF::gSF_SortHisto(bool colour) {
    std::sort(gSF_sort.begin(), gSF_sort.end(), compare);
    int nOfPoints = gSF_sort.size();
    
	//nOfPoints must always be an even number, or something is terribly wrong
	
	if (nOfPoints % 2 != 0) {
		std::cout <<"This is a bug! The number of gSF data points is not an even number! " <<std::endl;
		exit(0);
	}
	
	TMultiGraph *mg = new TMultiGraph();
	
	nOfPoints = nOfPoints / 2;
		
	double x1[nOfPoints], y1[nOfPoints];
	double dx1[nOfPoints], dy1[nOfPoints];
    	
	double x2[nOfPoints], y2[nOfPoints];
	double dx2[nOfPoints], dy2[nOfPoints];
    int id_1 = 0;
	int id_2 = 0;	
	for (int i = 0; i < 2*nOfPoints ; i++ ) {
        	
		if ( gSF_sort[i].peakID == 1 ) {
			x1[id_1] = gSF_sort[i].egamma;
			y1[id_1] = gSF_sort[i].value * sett->gSF_norm;
			dx1[id_1] = 1;
			dy1[id_1] = gSF_sort[i].dvalue * sett->gSF_norm;
			id_1++;
		}
		else {
			x2[id_2] = gSF_sort[i].egamma;
			y2[id_2] = gSF_sort[i].value * sett->gSF_norm;
			dx2[id_2] = 1;
			dy2[id_2] = gSF_sort[i].dvalue * sett->gSF_norm;
			id_2++;
		}
	}	
	TGraphErrors *gSFPlot_1 = new TGraphErrors(nOfPoints,x1,y1,dx1,dy1);
    gSFPlot_1->SetMarkerStyle(22);
    gSFPlot_1->SetMarkerSize(2);
    gSFPlot_1->SetMarkerColor(6);
	
	mg->Add(gSFPlot_1,"P");
	 
  	TGraphErrors *gSFPlot_2 = new TGraphErrors(nOfPoints,x2,y2,dx2,dy2);
    gSFPlot_2->SetMarkerStyle(22);
    gSFPlot_2->SetMarkerSize(2);
    if (colour) 
		gSFPlot_2->SetMarkerColor(7);
	else
		gSFPlot_2->SetMarkerColor(6);
    mg->Add(gSFPlot_2,"P");
	
	return mg;
}

//scales the sorted vector of gSF
void ShapeGSF::ScaleSort(double factor){
    
    for (int i = 0; i < gSF_sort.size(); i++) {
        gSF_sort[i].value = gSF_sort[i].value * factor;
        gSF_sort[i].dvalue = gSF_sort[i].dvalue * factor;
    }
}

void ShapeGSF::Scale(double factor){
    
    for (int i = 0; i < gSF_matrix->GetYBins(); i++) {
        gSF[i].value1 = gSF[i].value1 * factor;
        gSF[i].value2 = gSF[i].value2 * factor;
        gSF[i].dvalue1 = gSF[i].dvalue1 * factor;
        gSF[i].dvalue2 = gSF[i].dvalue2 * factor;
    }
}

