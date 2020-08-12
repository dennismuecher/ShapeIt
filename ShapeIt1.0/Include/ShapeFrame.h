#ifndef SHAPEFRAME_H
#define SHAPEFRAME_H

#include "../Include/ShapeSetting.h"
#include "../Include/ShapeMatrix.h"
#include "../Include/ShapeGSF.h"
#include "../Include/ShapeChi2.h"
#include "../Include/ShapeRho.h"

#include <TQObject.h>
#include <RQ_OBJECT.h>
#include <TGraphErrors.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TGComboBox.h>

#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <TGNumberEntry.h>
#include <TString.h>
#include <TApplication.h>
#include <TGLayout.h>
#include <TGWindow.h>
#include <TGLabel.h>
#include <TPave.h>
#include <TH2.h>
#include <TGFileDialog.h>
#include <TGDockableFrame.h>
#include <TGMenu.h>
#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <algorithm>


class TGWindow;
class TGMainFrame;
class TRootEmbeddedCanvas;

enum ETestCommandIdentifiers {
    M_FILE_ABOUT,
    M_FILE_OPEN,
    M_FILE_SAVE,
    M_FILE_SAVEAS,
    M_FILE_PRINT,
    M_FILE_EXIT,
    M_SETTING_OPEN,
    M_SETTING_SAVE,
    M_SETTING_SAVEAS,
    M_SETTING_PRINT,
    M_SETTING_OSLO,
    M_SETTING_RHO,
	M_SETTING_EFFI,
    M_SETTING_WIDTHRESET,
    M_DISPLAY_MAT,
    M_DISPLAY_DIAG,
    M_DISPLAY_DIAGCUBE,
    M_DISPLAY_PROJTOT,
    M_DISPLAY_PROJBIN,
    M_DISPLAY_GSF,
    M_DISPLAY_TRAFO,
    M_DISPLAY_FITWIDTH,
    M_DISPLAY_PRINT,
    M_DISPLAY_BAND,
    M_DISPLAY_RATIO,
    M_DISPLAY_GRF,
	M_DISPLAY_COLOUR,
	M_DISPLAY_VERBOSE0,
	M_DISPLAY_VERBOSE1,
	M_DISPLAY_VERBOSE2,
    M_DISPLAY_RHO,
    M_DISPLAY_EFFI

};

const char *filetypes[2] = {"ROOT files", "*.root"};
const char *filetypes_s[2] = {"dat files", "*.dat"};
const char *filetypes_t[2] = {"txt files", "*.txt"};

class ShapeFrame {
    RQ_OBJECT("ShapeFrame")
private:
    TGMainFrame         *fMain;
    const TGPicture *infoPic;               //the info picture ("File->About")
    string absPath;
    TRootEmbeddedCanvas *fEcanvas;
    TGNumberEntry       *energy[4]; //energies of the two discrete states; lower and upper limit
    TGNumberEntry       *exi[3];//excitation energies; lower and upper limit and interpolation energy
    TGNumberEntry       *bin[2];//size of integration bin; lower and upper limit
    TGNumberEntry       *nOfBins[2];//number of integration bins
    TGNumberEntry       *minContent; //mininmum number of counts an integration bin must have to be considered in the analysis
    TGNumberEntry       *scaling;           //scaling factor for gSF
    TGNumberEntry       *effCorr;              //efficiency factor for level 2
    TGNumberEntry       *transB;              //transformation B parameter field
    TGNumberEntry       *transAlpha;              //transformation alpha parameter field
    
    TGCompositeFrame *fBin;
    TGCompositeFrame *fInter;
    TGCompositeFrame *f1;
    TGCompositeFrame *f2;
    TGLayoutHints *fL1;
    TGLayoutHints *fL2;
    TGRadioButton* fR[6];
    TGCheckButton* OB[7];
    TGCheckButton* autoScale;
    TGPopupMenu* fMenuFile;
    TGPopupMenu* fSettingsFile;
    TGPopupMenu* fDisplayFile;
	TGPopupMenu* fVerboseMenu;
    
    TGComboBox *fMatrix;
    TGComboBox *fBinCombo;
	int binSelect = 1;						//the bin selected in fbin combo (projection spectrum)
    TGComboBox *fInterCombo;
    
    TGGroupFrame* fG[7];
    TGCompositeFrame* fEnergy[8];
    TGCompositeFrame* fTransform[2];
    TH1F *displayHisto;                      //points to the 1d histogram on display in the main canvas
    TH2 *inMatrix;							//the input matrix
    TFile* dataFile;							//the root file containing the input matrix
    ShapeSetting *sett;                  //the settings file
    ShapeMatrix *matrix;                    //the matrix object
    void SetupMenu();
    TLine *l[4];
    TBox *bgBox[4];                         //boxes indicating background regions
    TMultiGraph *mg;                        //the graph showing the gSF plots
    int displayMode = 1;
    TRootEmbeddedCanvas* fEMessage;
    string mname="";                           //filename of the root file without the pathname
    int status = 0;                             //tracks the status of ShapeIt: 0: no matrix loaded; 1: matrix loaded; 2: gSF data available
    void MessageBox(string title, string message);              //shows a dialog box with title and message and ok button
    TH1D* diagHisto;
    double histX1, histX2;                      //current selected x1 and x2 coordinates in 1d histograms
    void InfoWindow();                          //displays the Info Window from the file menu
    ShapeGSF *fitGSF;							//contains the gSF results from the data
    ShapeGSF *litGSF;                           //contains the gSF results from literature values
    bool gSF_band = false;                    //if true, display gSF as error band
    
public:
    ShapeFrame(const TGWindow *p,UInt_t w,UInt_t h, const string path);

    virtual ~ShapeFrame();
    int GetBinSelect() {return binSelect;}						//returns binSelect
	void SetBinSelect(int tBinSelect) {binSelect = tBinSelect;}		//set binSelect
	void DoDraw();
    void DoNumberEntry();
    void DoRadio();
    void HandleMenu(Int_t id);
    void CloseWindow();
    int MatrixSelector();										//updates the matrix selector and returns the index of the matrix saved in the current settings file; returns zero if no such matrix exists
    void UpdateSetting(ShapeSetting *sett_t);                   //updates a settings file
    void UpdateGuiSetting(ShapeSetting *sett_t);
    void MatrixSelect(Int_t mnr);
    void BinSelect(Int_t sbin);
	void HandleVerboseMenu(int vLevel);							//takes care about the verbose level menu

    void fBinComboDraw(TGComboBox *combo);
    void UpdateDisplay();
    void UpdateDisplay(int display);
    void DrawMarker();
    void HandleMyCanvas(Int_t a,Int_t b,Int_t c,TObject* obj);
    void ShapeItBaby();
    void PrintMessage();
	void ShowGraph();							//displays gSF results with literature values, resonance fit etc
	void Scale(Double_t scale);					//scale results of gSF and refresh display
    double AutoScale(int mode);                //auto-scales either gSF of data to literature (mode = 0) or literature to data (mode = 1)
    TGFileInfo fi;                              //file containing matrix
    TMultiGraph *wgraph ;
    double lit_chi2 = 0;                    //value of chi2 fit of lit gSF to fit gSF
};
#endif
