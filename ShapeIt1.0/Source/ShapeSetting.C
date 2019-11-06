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

//calculates the number of bins; the last bin might be smaller than sett->exi_size!
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


void ShapeSetting::SaveSettings() {
    ofstream outfile;
    outfile.open (settFileName.c_str());
    outfile << dataFileName << "\n";
    outfile << osloFileName << "\n";
    outfile << "matrixName " << matrixName<<"\n";
    outfile << "MeV: " << MeV<<"\n";
    outfile << "mode " << mode << "\n";
    outfile << "doInterpol " << doInterpol<<"\n";
    outfile << "doOslo " << doOslo<<"\n";
    outfile << "doSlidingWindow " << doSlidingWindow <<"\n";
    outfile << "doBinVariation " << doBinVariation <<"\n";
    outfile << "doBackground " << doBackground <<"\n";
    outfile << "verbose " << verbose<<"\n";
    outfile << "gSF_norm " << gSF_norm<<"\n";
    outfile << "level1 " << levEne[0] <<" "<<levEne[1] <<"\n";
    outfile << "bg_level1 " << bgEne[0][0] <<" "<< bgEne[0][1] <<" "<< bgEne[0][2] <<" "<< bgEne[0][3] <<"\n";
    outfile << "bg_level2 " << bgEne[1][0] <<" "<< bgEne[1][1] <<" "<< bgEne[1][2] <<" "<< bgEne[1][3] <<"\n";
    outfile << "level2 " << levEne[2] <<" "<<levEne[3] <<"\n";
    outfile << "excitation " << exiEne[0] <<" "<<exiEne[1] <<"\n";
    outfile << "excitation_bin_1 " << exi_size[0] <<"\n";
    outfile << "excitation_bin_2 " << exi_size[1] <<"\n";
    outfile << "nOfBins " << nOfBins <<"\n";
    outfile << "eff_corr " << eff_corr <<"\n";
    outfile << "interPoint " << interPoint <<"\n";
    outfile << "minCounts " << minCounts <<"\n";
    outfile << "doWidthCal " << doWidthCal <<"\n";
    outfile << "widthCal " << widthCal[0][0] <<" "<< widthCal[0][1] <<" "<< widthCal[1][0] <<" "<< widthCal[1][1] <<"\n";
   
    outfile.close();
    if (verbose)
        std::cout <<"Successfully saved Settings to file " << settFileName <<std::endl <<std::endl;
}

void ShapeSetting::ReadSettings() {
    if (verbose)
        std::cout <<"\nREADING FROM INPUT FILE: " <<endl;
    std::ifstream inp (settFileName.c_str());
    if (inp.is_open()) {
        string word, line;
        //get the root matrix filename; this goes extra because of issues wih absolute paths containing white spaces
        getline(inp,line);
        dataFileName = line;
        //get the literature data filename; this goes extra because of issues wih absolute paths containing white spaces
        getline(inp,line);
        osloFileName = line;
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
            if (word == "doSlidingWindow" ) isstr >> doSlidingWindow ;
            if (word == "doBinVariation" ) isstr >> doBinVariation ;
            if (word == "doBackground" ) isstr >> doBackground ;
            if (word == "verbose" ) isstr >> verbose;
            if (word == "gSF_norm" ) isstr >> gSF_norm;
            if (word == "level1" ) { isstr >> levEne[0]; isstr >>levEne[1];}
            if (word == "level2" ) { isstr >> levEne[2]; isstr >>levEne[3];}
            if (word == "bg_level1" ){ isstr >> bgEne[0][0]; isstr >> bgEne[0][1]; isstr >> bgEne[0][2]; isstr >> bgEne[0][3];}
            if (word == "bg_level2" ){ isstr >> bgEne[1][0]; isstr >> bgEne[1][1]; isstr >> bgEne[1][2]; isstr >> bgEne[1][3];}
            if (word == "excitation" ) { isstr >> exiEne[0]; isstr >>exiEne[1];}
            if (word == "excitation_bin_1" ) isstr >> exi_size[0] ;
            if (word == "excitation_bin_2" ) isstr >> exi_size[1] ;
            if (word == "nOfBins" ) isstr >> nOfBins ;
            if (word == "eff_corr" ) isstr >> eff_corr ;
            if (word == "interPoint" ) isstr >> interPoint ;
            if (word == "minCounts" ) isstr >> minCounts ;
            if (word == "doWidthCal" ) isstr >> doWidthCal ;
            if (word == "widthCal" ){ isstr >> widthCal[0][0]; isstr >> widthCal[0][1]; isstr >> widthCal[1][0]; isstr >> widthCal[1][1];}
        
        }
        if (verbose)
                PrintSettings();
                
    }
}


int ShapeSetting::BinToSize() {
    int size = (int) (exiEne[1] - exiEne[0] ) / nOfBins;
    return size;
}

int ShapeSetting::BinToSize(int n) {
    int size = (int) (exiEne[1] - exiEne[0] ) / n;
    return size;
}


void ShapeSetting::PrintSettings(){
    std::cout  << "root file name " << dataFileName<<"\n";
    std::cout  << "matrixName " << matrixName<<"\n";
    std::cout  << "Literature values " << osloFileName<<"\n";
    std::cout  << "MeV: " << MeV<<"\n";
    std::cout  << "mode " << mode << "\n";
    std::cout  << "doInterpol " << doInterpol<<"\n";
    std::cout  << "doOslo " << doOslo<<"\n";
    std::cout  << "doSlidingWindow " << doSlidingWindow <<"\n";
    std::cout  << "doBinVariation " << doBinVariation <<"\n";
    std::cout  << "doBackground " << doBackground <<"\n";
    std::cout  << "verbose " << verbose<<"\n";
    std::cout  << "gSF norm " << gSF_norm<<"\n";
    std::cout  << "level1 " << levEne[0] <<" "<<levEne[1] <<"\n";
    std::cout  << "left background level1 " << bgEne[0][0] <<"-" << bgEne[0][1] <<"\n";
    std::cout  << "right background level1 " << bgEne[0][2] <<"-" << bgEne[0][3] <<"\n";
    std::cout  << "left background level2 " << bgEne[1][0] <<"-" << bgEne[1][1] <<"\n";
    std::cout  << "right background level2 " << bgEne[1][2] <<"-" << bgEne[1][3] <<"\n";
    std::cout  << "level2 " << levEne[2] <<" "<<levEne[3] <<"\n";
    std::cout  << "excitation " << exiEne[0] <<" "<<exiEne[1] <<"\n";
    std::cout  << "excitation_bin_1 " << exi_size[0] <<"\n";
    std::cout  << "excitation_bin_2 " << exi_size[1] <<"\n";
    std::cout  << "nOfBins " << nOfBins <<"\n";
    std::cout  << "eff_corr " << eff_corr <<"\n";
    std::cout  << "interPoint " << interPoint <<"\n";
    std::cout  << "minCounts " << minCounts <<"\n";
    std::cout  << "doWidthCal " << doWidthCal <<"\n";
    std::cout  << "widthCal " << widthCal[0][0] <<" "<< widthCal[0][1] <<" "<< widthCal[1][0] <<" "<< widthCal[1][1] <<"\n";
  
}

