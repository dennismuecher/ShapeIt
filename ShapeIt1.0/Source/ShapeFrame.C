#include "../Include/ShapeFrame.h"
#include "../Source/ShapeInfo.C"
#include "ShapeDialogAlpha.C"

double glo(double *x, double *par){
  //par[0]: sigma
  //par[1]: width
  //par[2]: energy
  double Tf = 0.86; //CT model, T at final states (can be adjusted)
  //double Tf = par[3];
  double a = 8.674E-8; //1/3pi^2hbar^2c^2 constant, in MeV^2 mb^-1
  double b = 39.4784; //4pi^2 constant
  double gamma = (par[1]*(pow(x[0] /1000.,2)+b*pow(Tf,2)))/pow(par[2],2);
  double gamma0 = b*pow(Tf,2)*par[1]/pow(par[2],2);
  double term1 = ((x[0]*gamma / 1000.)/(pow(pow(x[0] / 1000. ,2)-pow(par[2],2),2)+pow(x[0] / 1000.,2)*pow(gamma,2)));
  double term2 = (0.7*gamma0)/(pow(par[2],3));
  double f = a*par[0]*par[1]*(term1 + term2);
  return f;
}

ShapeFrame::ShapeFrame(const TGWindow *p,UInt_t w,UInt_t h, const string path) {
    // Create a main frame
    fMain = new TGMainFrame(p,w,h);
  
    //take care about closing the window
    
    fMain->Connect("CloseWindow()", "ShapeFrame", this, "CloseWindow()");
    
    //set status to zero
    status = 0;
    // New setting object
    sett = new ShapeSetting();
    
    //store path name
    absPath = path;
    
    SetupMenu();
    
    //get resolution of display to sset window size
    Int_t qq;
    UInt_t ww;
    UInt_t hh;
    gVirtualX->GetWindowSize(gVirtualX->GetDefaultRootWindow(), qq, qq, ww, hh);
    
    //check if display width resolution makes sense; if so, set ww to 80% of display size
    if ( ww < 800 || ww > 6000)
        ww = 800;
    else
        ww = 0.8*ww;
    
    if ( hh >= 600 || hh < 400 )
        hh = 600;
    
    TGCompositeFrame* fSuper = new TGCompositeFrame(fMain, ww, hh, kHorizontalFrame | kFixedWidth);
    
    f1 = new TGCompositeFrame(fSuper, 300, 700, kVerticalFrame | kFixedWidth);
    f2 = new TGCompositeFrame(fSuper, 300, 700, kVerticalFrame);
    
    fL1 = new TGLayoutHints( kLHintsLeft | kLHintsExpandY,
                            2, 2, 3, 0);
    fL2 = new TGLayoutHints( kLHintsLeft | kLHintsExpandY | kLHintsExpandX,
                            2, 5, 5, 2);
    TGLayoutHints *fR2 = new TGLayoutHints( kLHintsLeft | kLHintsExpandY,
                            2, 5, 5, 2);
    TGLayoutHints* fL3 = new TGLayoutHints( kLHintsLeft, 2, 5, 0, 2);
    
    //left half of main window
    
    TGCompositeFrame* f3 = new TGCompositeFrame(f1, 300, 300, kHorizontalFrame);
    
    fG[0] = new TGGroupFrame(f3, new TGString("Mode"),kVerticalFrame|kRaisedFrame);
    fG[1] = new TGGroupFrame(f1, new TGString("Energies (all in keV)"),kVerticalFrame|kRaisedFrame);
    fG[2] = new TGGroupFrame(f1, new TGString("Options"),kVerticalFrame|kRaisedFrame);
    fG[3] = new TGGroupFrame(f1, new TGString("Integration bin"),kVerticalFrame|kRaisedFrame);
    fG[4] = new TGGroupFrame(f3, new TGString("Input Matrix"),kVerticalFrame|kRaisedFrame);
    fG[5] = new TGGroupFrame(f1, new TGString("Display"),kVerticalFrame|kRaisedFrame);
    
    for (int i = 0; i < 3; i++)
        fEnergy[i] = new TGCompositeFrame(fG[1], 1, 1, kHorizontalFrame);
    for (int i = 3; i < 8; i++)
        fEnergy[i] = new TGCompositeFrame(fG[3], 1, 1, kHorizontalFrame);
    
    //Mode
    
    fR[0] = new TGRadioButton(fG[0], new TGHotString("Integration"), 11);
    fR[0]->SetState(kButtonDown);
    fR[1] = new TGRadioButton(fG[0], new TGHotString("Autofit"), 12);
    
    //Matrix selector
    fMatrix = new TGComboBox(fG[4], 120);
    fMatrix->AddEntry("no file open",1);
    fMatrix->Resize(120,20);
    fMatrix->Select(1);
    fMatrix->SetEnabled(false);
    fG[4]->AddFrame(fMatrix, fL2);
    //Energy Settings
    energy[0] = new TGNumberEntry(fEnergy[0], 410, 9,1, TGNumberFormat::kNESInteger,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,-10000, 99999);
    fEnergy[0]->AddFrame(energy[0], fL2);
    energy[0]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    TGLabel *l1 = new TGLabel(fEnergy[0], "< Level 1 <");
    fEnergy[0]->AddFrame(l1, fL2);
    energy[1] = new TGNumberEntry(fEnergy[0], 660, 9,2, TGNumberFormat::kNESInteger,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,-10000, 99999);
    energy[1]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[0]->AddFrame(energy[1], fL2);
    fG[1]->AddFrame(fEnergy[0], fL2);

    energy[2] = new TGNumberEntry(fEnergy[1], 1000, 9,3, TGNumberFormat::kNESInteger,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,-10000, 99999);
    energy[2]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[1]->AddFrame(energy[2], fL2);
    TGLabel *l2 = new TGLabel(fEnergy[1], "< Level 2 <");
    fEnergy[1]->AddFrame(l2, fL2);
    energy[3] = new TGNumberEntry(fEnergy[1], 1300, 9,4, TGNumberFormat::kNESInteger,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,-10000, 99999);
    energy[3]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[1]->AddFrame(energy[3], fL2);
    fG[1]->AddFrame(fEnergy[1], fL2);
   
    //initialize background regions
    double bg1[4]={260,360,700,800};
    double bg2[4]={850,950,1350,1450};

    sett->setBgEne1(bg1);
    sett->setBgEne2(bg2);
    
    exi[0] = new TGNumberEntry(fEnergy[2], 2500, 9,5, TGNumberFormat::kNESInteger,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,0, 99999);
    exi[0]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[2]->AddFrame(exi[0], fL2);
    TGLabel *l3 = new TGLabel(fEnergy[2], "< Excitation <");
    fEnergy[2]->AddFrame(l3, fL2);
    exi[1] = new TGNumberEntry(fEnergy[2], 7000, 9,6, TGNumberFormat::kNESInteger,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,0, 99999);
    exi[1]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[2]->AddFrame(exi[1], fL2);
    fG[1]->AddFrame(fEnergy[2], fL2);

    TGLabel *l4 = new TGLabel(fEnergy[3], "Bin size [keV]");
    fEnergy[3]->AddFrame(l4, fR2);
    bin[0] = new TGNumberEntry(fEnergy[3], 400, 9,7, TGNumberFormat::kNESInteger,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,1, 99999);
    bin[0]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[3]->AddFrame(bin[0], fR2);
    
    TGLabel *l9 = new TGLabel(fEnergy[3], "-");
    fEnergy[3]->AddFrame(l9, fR2);
    bin[1] = new TGNumberEntry(fEnergy[3], 800, 9,11, TGNumberFormat::kNESInteger,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,1, 99999);
    bin[1]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[3]->AddFrame(bin[1], fR2);
    
    fG[3]->AddFrame(fEnergy[3], fR2);
    
    TGLabel *l5 = new TGLabel(fEnergy[4], "    Nr. of Bins");
    fEnergy[4]->AddFrame(l5, fR2);
    nOfBins[0] = new TGNumberEntry(fEnergy[4], 12, 9,8, TGNumberFormat::kNESInteger,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,1, 99999);
    nOfBins[0]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[4]->AddFrame(nOfBins[0], fR2);
    
    TGLabel *l10 = new TGLabel(fEnergy[4], "-");
    fEnergy[4]->AddFrame(l10, fR2);
    
    nOfBins[1] = new TGNumberEntry(fEnergy[4], 6, 9,12, TGNumberFormat::kNESInteger,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,1, 99999);
    nOfBins[1]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[4]->AddFrame(nOfBins[1], fR2);
    
    fG[3]->AddFrame(fEnergy[4], fR2);

    bin[1]->SetState(false);
    nOfBins[1]->SetState(false);
    
    TGLabel *l6 = new TGLabel(fEnergy[5], "  Min. Counts");
    fEnergy[5]->AddFrame(l6, fR2);
    minContent = new TGNumberEntry(fEnergy[5], 50, 9,9, TGNumberFormat::kNESInteger,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,0, 99999);
    minContent->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[5]->AddFrame(minContent, fR2);
    fG[3]->AddFrame(fEnergy[5], fL2);
    
    TGLabel *l7 = new TGLabel(fEnergy[6], "   gSF scaling");
    fEnergy[6]->AddFrame(l7, fR2);
    scaling = new TGNumberEntry(fEnergy[6], 0.012, 9,10, TGNumberFormat::kNESReal,TGNumberFormat::kNEAPositive,TGNumberFormat::kNELLimitMinMax,0, 9999999);
    scaling->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[6]->AddFrame(scaling, fR2);
    fG[3]->AddFrame(fEnergy[6], fL2);
    autoScale = new TGCheckButton(fEnergy[6], new TGHotString("Auto"), 40);
    autoScale->Connect("Clicked()", "ShapeFrame", this, "DoRadio()");
    fEnergy[6]->AddFrame(autoScale, fR2);
    
    TGLabel *l8 = new TGLabel(fEnergy[7], "      Eff. Corr.");
    fEnergy[7]->AddFrame(l8, fR2);
    effCorr = new TGNumberEntry(fEnergy[7], 1, 9,13, TGNumberFormat::kNESReal,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,0, 9999999);
    effCorr->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[7]->AddFrame(effCorr, fR2);
    fG[3]->AddFrame(fEnergy[7], fL2);
    
    //options
    OB[0] = new TGCheckButton(fG[2], new TGHotString("Sewing Interpolation"), 14);
    OB[1] = new TGCheckButton(fG[2], new TGHotString("Display expectation"), 15);
    OB[2] = new TGCheckButton(fG[2], new TGHotString("Sliding window variation"), 21);
    OB[3] = new TGCheckButton(fG[2], new TGHotString("Bin size variation"), 22);
    OB[4] = new TGCheckButton(fG[2], new TGHotString("Background subtraction"), 23);
    OB[5] = new TGCheckButton(fG[2], new TGHotString("Use Width Calibration"), 24);
    OB[5]->SetEnabled(false);
    OB[6] = new TGCheckButton(fG[2], new TGHotString("Use Efficiency Calibration from File"), 25);
    OB[6]->SetEnabled(false);
    for (int i = 0 ; i < 7; i++) {
        OB[i]->Connect("Clicked()", "ShapeFrame", this, "DoRadio()");
         fG[2]->AddFrame(OB[i], fL3);
    }
    fG[0]->AddFrame(fR[0], fL3);
    fG[0]->AddFrame(fR[1], fL3);
    f3->AddFrame(fG[0], fL3 );
    f3->AddFrame(fG[4], fL3 );
    f1->AddFrame(f3, fL3 );
    f1->AddFrame(fG[2], fL3 );
    f1->AddFrame(fG[1], fL3 );
    f1->AddFrame(fG[3], fL3);
    
    //display settings
    fBin = new TGCompositeFrame(fG[5], 1, 1, kHorizontalFrame);
    fR[2] = new TGRadioButton(fBin, new TGHotString("Diagonal Projection"), 18);
    fBin->AddFrame(fR[2], new TGLayoutHints( kLHintsTop, 2, 2 , 3, 2));
    
    //Combo box for displaying bin spectra
    fBinCombo = new TGComboBox(fBin, 80);
    fBinComboDraw(fBinCombo);
    fBinCombo->SetEnabled(false);
    fBin->AddFrame(fBinCombo, new TGLayoutHints( kLHintsTop, 0, 0, 1, 0));
   
    fG[5]->AddFrame(fBin);
    
    f1->AddFrame(fG[5], fL3);
    for (int i = 0 ; i < 3; i++)
        fR[i]->Connect("Clicked()", "ShapeFrame", this, "DoRadio()");

    //show ShapeIt Button
    TGPictureButton *fPicture = new TGPictureButton(f1,
                                                    gClient->GetPicture("logo2.jpg"), 20);
    fPicture->Connect("Clicked()", "ShapeFrame", this, "ShapeItBaby()");
    f1->AddFrame(fPicture,new TGLayoutHints(kLHintsLeft,
                                            1, 1, 1, 1));
    
    //right half of main window
    
    //the drawing window
    fEcanvas = new TRootEmbeddedCanvas("Ecanvas",f2,100,100);
    fEcanvas->GetCanvas()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",
                              "ShapeFrame", this, "HandleMyCanvas(Int_t,Int_t,Int_t,TObject*)");
    f2->AddFrame(fEcanvas, fL2);
    
    fSuper->AddFrame(f1, fL1);
    fSuper->AddFrame(f2, fL2);
  
    fMain->AddFrame(fSuper, fL2);
    // Set a name to the main frame
    fMain->SetWindowName("unsaved");
    
    // Map all subwindows of main frame
    fMain->MapSubwindows();
    
    // Initialize the layout algorithm
    fMain->Resize(fMain->GetDefaultSize());
    // Map main frame
    fMain->MapWindow();
    // store default settings in settings object
    UpdateSetting(sett);
}


//updates the GUI settings according to the settings file
void ShapeFrame::UpdateGuiSetting(ShapeSetting *sett_t)
{
    for (int i = 0; i < 4; i++)
        energy[i]->SetNumber(sett_t->levEne[i]);
    for (int i = 0; i < 2; i++)
         exi[i]->SetNumber(sett_t->exiEne[i]);
    bin[0]->SetNumber(sett_t->exi_size[0]);
    //acivate bin size variation TGNumberEntries
    if (sett->doBinVariation) {
        bin[1]->SetState(true);
        nOfBins[1]->SetState(true);
    }
    else {
        bin[1]->SetState(false);
        nOfBins[1]->SetState(false);
    }
    bin[1]->SetNumber(sett_t->exi_size[1]);
    nOfBins[0]->SetNumber(sett_t->nOfBins);

    minContent->SetNumber(sett_t->minCounts);
    scaling->SetNumber(sett_t->gSF_norm);
    effCorr->SetNumber(sett_t->eff_corr);
        
    for (int i = 0; i < 2; i++)
        fR[i]->SetState(kButtonUp);
    if (sett_t->mode == 1)
        fR[0]->SetState(kButtonDown);
    if (sett_t->mode == 2)
        fR[1]->SetState(kButtonDown);
    
    for (int i = 0; i < 5; i++)
        OB[i]->SetState(kButtonUp);
    if (status > 2) {
        OB[5]->SetEnabled(true);
         OB[5]->SetState(kButtonUp);
    }
    
    if (sett_t->doInterpol)
        OB[0]->SetState(kButtonDown);
    if (sett_t->doOslo)
        OB[1]->SetState(kButtonDown);
    if (sett_t->doSlidingWindow)
        OB[2]->SetState(kButtonDown);
    if (sett_t->doBinVariation)
        OB[3]->SetState(kButtonDown);
    if (sett_t->doBackground)
        OB[4]->SetState(kButtonDown);
    if (sett_t->doWidthCal)
        OB[5]->SetState(kButtonDown);
    
    if (sett->doAutoScale) {
        autoScale->SetState(kButtonDown);
        scaling->SetState(false);
    }
    else {
        autoScale->SetState(kButtonUp);
        scaling->SetState(true);
    }
    
    if ( sett_t->effiFileName.empty()  )
        OB[6]->SetEnabled(false);
    else
        OB[6]->SetEnabled(true);
    if (sett_t->doEffi)
        OB[6]->SetState(kButtonDown);
    else
        OB[6]->SetState(kButtonUp);
    
    
    
}


//updates the settings object using the current GUI settings
void ShapeFrame::UpdateSetting(ShapeSetting *sett_t)
{

    for (int i = 0; i < 4; i++)
        sett_t->levEne[i] = energy[i]->GetNumber();
    for (int i = 0; i < 2; i++)
        sett_t->exiEne[i] = exi[i]->GetNumber();
    
    sett_t->exi_size[0] = bin[0]->GetNumber();
    sett_t->exi_size[1] = bin[1]->GetNumber();
    sett_t->nOfBins = nOfBins[0]->GetNumber();
    sett_t->minCounts = minContent->GetNumber();
    sett_t->gSF_norm = scaling->GetNumber();
    sett_t->eff_corr = effCorr->GetNumber();
    
    //options
    if (OB[0]->GetState() == kButtonDown)
        sett_t->doInterpol = true;
    else
        sett_t->doInterpol = false;
    
    if (OB[1]->GetState() == kButtonDown)
        sett_t->doOslo = true;
    else
        sett_t->doOslo = false;
   
    if (OB[2]->GetState() == kButtonDown)
        sett_t->doSlidingWindow = true;
    else
        sett_t->doSlidingWindow = false;
    
    if (OB[3]->GetState() == kButtonDown)
        sett_t->doBinVariation = true;
    else
        sett_t->doBinVariation = false;
    
    if (OB[4]->GetState() == kButtonDown)
        sett_t->doBackground = true;
    else
        sett_t->doBackground = false;
    
    if (OB[5]->GetState() == kButtonDown)
        sett_t->doWidthCal = true;
    else
        sett_t->doWidthCal = false;
    if (status < 3)
        OB[5]->SetEnabled(false);
}


void ShapeFrame::SetupMenu() {
    // Create menubar and popup menus.
    TGDockableFrame* fMenuDock = new TGDockableFrame(fMain);
    fMain->AddFrame(fMenuDock, new TGLayoutHints(kLHintsExpandX, 0, 0, 1, 0));
    // fMain->ChangeOptions((fMain->GetOptions() & ~kVerticalFrame) | kHorizontalFrame);
    //File Menu
    fMenuFile = new TGPopupMenu(gClient->GetRoot());
    fMenuFile->AddEntry("&About ShapeIt", M_FILE_ABOUT);
    fMenuFile->AddEntry("&Load Matrix Root File...", M_FILE_OPEN);
    fMenuFile->AddSeparator();
    fMenuFile->AddSeparator();
    fMenuFile->AddEntry("E&xit", M_FILE_EXIT);
    
    TGLayoutHints* fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX);
    TGLayoutHints* fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
    TGLayoutHints* fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);
    
    TGMenuBar* fMenuBar = new TGMenuBar(fMenuDock, 1, 1, kHorizontalFrame);
    fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
    
    //Settings Menu
    fSettingsFile = new TGPopupMenu(gClient->GetRoot());
    fSettingsFile->AddEntry("&Load Settings ...", M_SETTING_OPEN);
    fSettingsFile->AddEntry("&Save Settings", M_SETTING_SAVE);
    fSettingsFile->AddEntry("S&ave Settings as ...", M_SETTING_SAVEAS);
    fSettingsFile->AddEntry("&Print Settings to terminal", M_SETTING_PRINT);
    fSettingsFile->AddSeparator();
    fSettingsFile->AddEntry("L&oad Literature Values gSF", M_SETTING_OSLO);
    fSettingsFile->AddEntry("Load Literature Values &rho", M_SETTING_RHO);
    fSettingsFile->AddEntry("L&oad Efficiency Calibration", M_SETTING_EFFI);
    fSettingsFile->AddSeparator();
    fSettingsFile->AddEntry("&Reset peak width calibration", M_SETTING_WIDTHRESET);
    fSettingsFile->AddSeparator();
    fSettingsFile->AddEntry("&Transformation ...", M_SETTING_TRAFO);

    fMenuBar->AddPopup("&Settings", fSettingsFile, fMenuBarItemLayout);

    //cascade Menu for verbose level
    fVerboseMenu = new TGPopupMenu(gClient->GetRoot());
    fVerboseMenu->AddEntry("&Silent", M_DISPLAY_VERBOSE0);
    fVerboseMenu->AddEntry("Some &Info", M_DISPLAY_VERBOSE1);
    fVerboseMenu->AddEntry("&All the Details", M_DISPLAY_VERBOSE2);
    fVerboseMenu->CheckEntry(M_DISPLAY_VERBOSE0);
    
    //Display Menu
    fDisplayFile = new TGPopupMenu(gClient->GetRoot());
    fDisplayFile->AddEntry("&Input Matrix", M_DISPLAY_MAT);
    fDisplayFile->AddEntry("&Diag vs Excitation energy", M_DISPLAY_DIAG);
    fDisplayFile->AddEntry("Diag &Cube vs Excitation energy", M_DISPLAY_DIAGCUBE);
    fDisplayFile->AddEntry("Diag &Projection (total)", M_DISPLAY_PROJTOT);
    fDisplayFile->AddEntry("Dia&g Projection (bin)", M_DISPLAY_PROJBIN);
    fDisplayFile->AddSeparator();
    fDisplayFile->AddEntry("&Plot Efficiency Data from File", M_DISPLAY_EFFI);
    fDisplayFile->AddEntry("&Show Peak Width", M_DISPLAY_FITWIDTH);
    fDisplayFile->AddEntry("Show Peak &Ratios", M_DISPLAY_RATIO);
    fDisplayFile->AddSeparator();
    fDisplayFile->AddEntry("Show gSF &average", M_DISPLAY_AVG);
    fDisplayFile->AddEntry("Show gSF single data", M_DISPLAY_SINGLE);
    fDisplayFile->CheckEntry(M_DISPLAY_SINGLE);
    fDisplayFile->AddEntry("Use &Two-Colour Display", M_DISPLAY_COLOUR);
    fDisplayFile->CheckEntry(M_DISPLAY_COLOUR);
    
    fDisplayFile->AddSeparator();
    fDisplayFile->AddEntry("&Fit Giant Resoance Formula", M_DISPLAY_GRF);
    fDisplayFile->AddEntry("Plo&t chi2 graph for transformation", M_DISPLAY_TRAFO);
    
    fDisplayFile->AddSeparator();
    fDisplayFile->AddEntry("&Print Table of gSF results", M_DISPLAY_PRINT);
    
    fDisplayFile->AddSeparator();
    fDisplayFile->AddEntry("&Draw level density", M_DISPLAY_RHO);
    
    fDisplayFile->DisableEntry(M_DISPLAY_FITWIDTH);
    fDisplayFile->DisableEntry(M_DISPLAY_RATIO);
    
    //cascade Menu for verbose level
    fDisplayFile->AddPopup("&Verbose Level", fVerboseMenu);
    
    fMenuBar->AddPopup("&Display", fDisplayFile, fMenuBarItemLayout);

    
    fMenuDock->AddFrame(fMenuBar, fMenuBarLayout);
    
    fMenuFile->Connect("Activated(Int_t)", "ShapeFrame", this,"HandleMenu(Int_t)");
    fSettingsFile->Connect("Activated(Int_t)", "ShapeFrame", this,"HandleMenu(Int_t)");
    fDisplayFile->Connect("Activated(Int_t)", "ShapeFrame", this,"HandleMenu(Int_t)");

}

//fills a combo box with the number of bins
void ShapeFrame::fBinComboDraw(TGComboBox *combo)
{
    combo->SetEnabled(true);
    char name[50];
    combo->RemoveAll();
    for (int i = 1; i < sett->nOfBins+1 ; i++) {
        sprintf(name,"bin %d",i);
        combo->AddEntry(name,i);
    }
    combo->Resize(80,20);
    combo->Select(1);
}


void ShapeFrame::BinSelect(Int_t sbin)
{
    binSelect = sbin;
    diagHisto = matrix->GetDiagEx(sbin, mname);
    diagHisto->GetXaxis()->SetRangeUser(histX1, histX2);
    //display options for paper
    gStyle->SetOptStat(0);
    diagHisto->GetXaxis()->SetTitle("");
    diagHisto->GetYaxis()->SetTitle("");
    //diagHisto->GetXaxis()->SetRangeUser(200,1490);
    //diagHisto->GetXaxis()->SetLabelSize(2);
    diagHisto->GetXaxis()->SetNdivisions(5,5,0);
    diagHisto->GetYaxis()->SetNdivisions(5,5,0);
    diagHisto->Draw("Y+");
    
    //standard display
    //diagHisto->Draw();

    DrawMarker();
    TCanvas *fCanvas = fEcanvas->GetCanvas();
    fCanvas->cd();
    fCanvas->Modified();
    fCanvas->Update();
}

void ShapeFrame::DoRadio()
{
    // Handle radio buttons.
    
    TGButton *btn = (TGButton *) gTQSender;
    Int_t id = btn->WidgetId();
    if (id >= 11 && id <= 12) {
        for (int i = 0; i < 2; i++)
            if (fR[i]->WidgetId() != id)
                fR[i]->SetState(kButtonUp);
    }
    switch (id) {
        case 11: {
            sett->mode = 1;
            UpdateDisplay(displayMode);
            break;
        }
        case 12:
            sett->mode = 2;
            UpdateDisplay(displayMode);
            break;
        case 14:
            sett->doInterpol = !sett->doInterpol;
            break;
        case 15:
            sett->doOslo = !sett->doOslo;
            break;
        case 18:
            displayMode = 5;
            UpdateDisplay(displayMode);
            break;
        case 21:
            sett->doSlidingWindow = !sett->doSlidingWindow;
            break;

        case 22:
            sett->doBinVariation = !sett->doBinVariation;
            UpdateGuiSetting(sett);
            break;
        case 23:
            sett->doBackground = !sett->doBackground;
            DrawMarker();
            break;
        case 24:
            sett->doWidthCal = !sett->doWidthCal;
            break;
        case 25:
            sett->doEffi = !sett->doEffi;
            break;
        case 40:
            sett->doAutoScale = !sett->doAutoScale;
            ShowGraph();
            UpdateGuiSetting(sett);
            break;

    }
}

void ShapeFrame::MessageBox(string title, string message)
{
    int retval = 0;
    EMsgBoxIcon icontype = kMBIconAsterisk;
    new TGMsgBox(gClient->GetRoot(), fMain, title.c_str(), message.c_str(), icontype, kMBOk, &retval);
    
}

void ShapeFrame::DrawMarker() {
    
        TCanvas *fCanvas = fEcanvas->GetCanvas();
        for (int i = 0; i < 4; i++) {
            fCanvas->GetListOfPrimitives()->Remove(l[i]);
            fCanvas->GetListOfPrimitives()->Remove(bgBox[i]);
        }
        fCanvas->Modified();
        fCanvas->Update();
        if (displayMode == 4 || displayMode == 5) {
        for (int i = 0; i < 4; i++) {
            
            //draw verticl lines; there is a bug when using log-y scale, discussed here:
            // https://root-forum.cern.ch/t/getuymax-has-problems-with-log-scale/10511/4
            // but solution not implemented, yet
        
            l[i] = new TLine(sett->levEne[i],fCanvas->GetUymin(),sett->levEne[i],fCanvas->GetUymax());
            
            l[i]->SetLineColor(kRed);
            l[i]->SetLineWidth(2);

            l[i]->SetLineColorAlpha(kRed, 0.45);
            
            l[i]->Draw();
        }
        bgBox[0] = new TBox(sett->bgEne[0][0],fCanvas->GetUymin(),sett->bgEne[0][1],fCanvas->GetUymax());
        bgBox[1] = new TBox(sett->bgEne[0][2],fCanvas->GetUymin(),sett->bgEne[0][3],fCanvas->GetUymax());
        bgBox[2] = new TBox(sett->bgEne[1][0],fCanvas->GetUymin(),sett->bgEne[1][1],fCanvas->GetUymax());
        bgBox[3] = new TBox(sett->bgEne[1][2],fCanvas->GetUymin(),sett->bgEne[1][3],fCanvas->GetUymax());
        
        for (int i = 0; i < 4; i++) {
            if (i ==0 || i==1)
                bgBox[i]->SetFillColorAlpha(kBlue-9, 0.45);
            else if (i ==2 ||i ==3)
                bgBox[i]->SetFillColorAlpha(kBlue-10, 0.45);
            if (sett->doBackground)
                bgBox[i]->Draw();
        }
        fCanvas->Modified();
        fCanvas->Update();
    }
}

void ShapeFrame::HandleMyCanvas(Int_t a,Int_t b,Int_t c,TObject* obj) {
    //catch left or right mouse click...
    if (status == 0 )
        return;
    if (!(displayMode == 4 || displayMode == 5))
        return;
    if ( (a == kButton1Up || a==kButton2Up || a == kButton3Up )) {
        
        for (int i = 0; i < 4; i++)
            sett->levEne[i] = l[i]->GetX1();
        
        sett->bgEne[0][0] = bgBox[0]->GetX1(); sett->bgEne[0][1] = bgBox[0]->GetX2();
        sett->bgEne[0][2] = bgBox[1]->GetX1(); sett->bgEne[0][3] = bgBox[1]->GetX2();
        sett->bgEne[1][0] = bgBox[2]->GetX1(); sett->bgEne[1][1] = bgBox[2]->GetX2();
        sett->bgEne[1][2] = bgBox[3]->GetX1(); sett->bgEne[1][3] = bgBox[3]->GetX2();
        UpdateGuiSetting(sett);

        DrawMarker();

        //if autofit active and displayed, update fit
        if (sett->mode == 2 && displayMode == 5) {
            TCanvas *fCanvas = fEcanvas->GetCanvas();
            fCanvas->cd();
            fCanvas->Modified();
            fCanvas->Update();
            histX1 = gPad->GetUxmin();
            histX2 = gPad->GetUxmax();
            BinSelect(fBinCombo->GetSelected());

        }
    }
}

void ShapeFrame::DoDraw() {
    if (fMatrix->IsEnabled()) {
        matrix->GetInputMatrix(mname)->Draw("colz");
        DrawMarker();
    }
    TCanvas *fCanvas = fEcanvas->GetCanvas();
    fCanvas->cd();
    fCanvas->Update();
}

void ShapeFrame::DoNumberEntry() {
    TGNumberEntry *nbr = (TGNumberEntry *) gTQSender;
    Int_t id = nbr->WidgetId();
    if (fMatrix->IsEnabled()) {
        if (id > 0 || id < 5 ){
            UpdateSetting(sett);
            DrawMarker();
        }
        if (id == 5 || id == 6){
            UpdateSetting(sett);
            nOfBins[0]->SetNumber(sett->SizeToBin());
            nOfBins[1]->SetNumber(sett->SizeToBin(bin[1]->GetNumber()));
            UpdateSetting(sett);
        }
        if (id == 7){
            UpdateSetting(sett);
            nOfBins[0]->SetNumber(sett->SizeToBin());
            UpdateSetting(sett);
        }
        if (id == 8){
            UpdateSetting(sett);
            bin[0]->SetNumber(sett->BinToSize());
            UpdateSetting(sett);
        }
        
        //scaling factor has changed
        if (id == 10) {
            if (scaling->GetNumber() == 0 ) {
                scaling->SetNumber(0.9*scale_bak);
                std::cout <<"scaling zero event detected!" <<std::endl;
            }
            scale_bak = scaling->GetNumber();
             UpdateSetting(sett);
             ShowGraph();
         }
         
        if (id == 11){
            UpdateSetting(sett);
            if (bin[1]->GetNumber() <= bin[0]->GetNumber() + 50  )
                bin[1]->SetNumber(bin[0]->GetNumber() + 50);
            nOfBins[1]->SetNumber(sett->SizeToBin(bin[1]->GetNumber()));
            UpdateSetting(sett);
        }
        
        if (id == 12){
            UpdateSetting(sett);
             bin[1]->SetNumber(sett->BinToSize(nOfBins[1]->GetNumber()));
            if (bin[1]->GetNumber() <= bin[0]->GetNumber() + 50  ) {
                bin[1]->SetNumber(bin[0]->GetNumber() + 50);
                nOfBins[1]->SetNumber(sett->SizeToBin(bin[1]->GetNumber()));
            }
            if (nOfBins[1]->GetNumber() < 3) {
                nOfBins[1]->SetNumber(3);
                bin[1]->SetNumber(sett->BinToSize(nOfBins[1]->GetNumber()));
            }
            UpdateSetting(sett);
        }
        
        //Trnasformation has changed
        if ( (id == 31) || (id == 32) ) {
            ShowGraph();
        }
        
        if (id == 9 || id == 13 ){
            UpdateSetting(sett);
        }
        if (id >4 && id < 9) {
            UpdateSetting(sett);
            matrix->SetEne0(sett->exiEne[0]);
            matrix->SetEne1(sett->exiEne[1]);
            matrix->SetESize(sett->exi_size[0]);
            matrix->Diag();
            fBinComboDraw(fBinCombo);
            UpdateDisplay(displayMode);
        }
    }
}

void ShapeFrame::ShapeItBaby() {
    
    //check if matrix is loaded
    if (status == 0)
        return;
   
    //clean up matrix
    matrix->Reset();
    
    //set current values of excitation energies for matrix
    matrix->SetEne0( sett->exiEne[0] );
    matrix->SetEne1( sett->exiEne[1] );
    matrix->SetESize( sett->exi_size[0] );
    
    //update matrix
    matrix->Diag();
    
    int kmax = 6;    //number of steps in sliding window; should not be hard-coded....needs fixing
    
    //setting display mode
    UpdateDisplay(6);   //is this good for anything?
  
    //update settings accordings to the GUI settings
    UpdateSetting(sett);
   
    sett->nOfBins = sett->SizeToBin();

    //calculate gamma ray strength function
    gSF = new ShapeGSF(sett, matrix);
    gSF->FillgSF();

    //status update: will have values for gSF
    status = 2;
    
    //enable Menu entry showing fit results for peak width and peak ratios
    if (sett->mode == 2)
        fDisplayFile->EnableEntry(M_DISPLAY_FITWIDTH);
    
    fDisplayFile->EnableEntry(M_DISPLAY_RATIO);
    
    //loop over exi_size
    do {
        //sliding window
        for (int i = 2; i < kmax; i++) {
            if (!sett->doSlidingWindow) break;
            // the sliding window moves from the inital position in kmax steps to the end position, which is one bin to the "left"
            // the high energy is kept at the initital value, at all times
            matrix->SetEne0( sett->exiEne[0] - (double) (i-1) * matrix->GetESize() / (kmax-1));
              matrix->SetEne1( sett->exiEne[1] );
            
            //update matrix
            matrix->Diag();
            
            //recalculate gSF values
            gSF->FillgSF();
        }
        //check if bin size should be varried
        if (!sett->doBinVariation) break;
        
        //change bin size
        matrix->SetESize( matrix->GetESize() + 50 );
        
    } while (matrix->GetESize() <= sett->exi_size[1]);
    
    //display results
    ShowGraph();
}

void ShapeFrame::TransGraph()
{
    //update Literature value transformation settings
   
    sett->lit_alpha = AlphaDialog->GetAlphaTransform();
    sett->lit_norm =  AlphaDialog->GetBTransform();
    
    ShowGraph();
}

//displays the results for gSF, literature values, resonance fit etc.
void ShapeFrame::ShowGraph()
{
    //check status
    if (status < 2)
        return;
        
    TCanvas *fCanvas = fEcanvas->GetCanvas();
    fCanvas->Clear();
    
    //fit giant resonance formula if requested
    
    /*if (sett->doGRF) {
        TF1 *fitGDR = new TF1("fitGDR",glo,3000,6000,3);
        fitGDR->SetParameter(0,300);
        fitGDR->SetParameter(1,5);
        fitGDR->SetParameter(2,15);
        //fitGDR->SetParameter(3,0.5);
        fitGDR->SetLineColor(1);
        g->Fit(fitGDR," "," ",3000,6000);
    }*/
    
    //plot the multigraph
    gSF->getMultGraph()->Draw("A*");
    //gSF->getMultGraph();
    //gSF->Smooth(50)->Draw("P,same");
    fCanvas->Modified();
    fCanvas->Update();
}


void ShapeFrame::MatrixSelect(Int_t mnr)
{
    matrix->SetMatrix(mnr);
    matrix->GetInputMatrix(mname)->Draw("colz");
    matrix->Diag();
    TCanvas *fCanvas = fEcanvas->GetCanvas();
    fCanvas->cd();
    fCanvas->Modified();
    fCanvas->Update();
    std::cout <<"selected(int)!"<<std::endl;
}

int ShapeFrame::MatrixSelector()
{

    int length_max = 15;             //maximum length of string; used to resize selector
    std::vector <std::string> s;
     s.clear();
    s = matrix->GetMatrixName();
    int index = 0;

    if (s.size() > 0) {
        if (!fMatrix->IsEnabled())
            fMatrix->SetEnabled(true);
        fMatrix->RemoveAll();
        for (int i = 0; i < s.size(); i++ ) {
           fMatrix->AddEntry(s[i].c_str(), i+1);
            if ( s[i].length() > length_max ) {
                length_max = s[i].length();
            }
            
           if ( s[i] == sett->matrixName )
               index = i+1;
        }
        if (length_max > 15)
          fMatrix->Resize(length_max*8,20);
        
        fMatrix->Connect("Selected(Int_t)", "ShapeFrame", this, "MatrixSelect(Int_t)");
    }
    
    fMain->MapSubwindows();
    fMain->Resize(fMain->GetDefaultSize());
    fMain->MapWindow();
    return index;
    
}

//takes care about setting the verbose level and displaying the verbose cascade menu correctly
void ShapeFrame::HandleVerboseMenu(int vLevel) {
    
    switch (vLevel) {
        case 0: {
            fVerboseMenu->CheckEntry(M_DISPLAY_VERBOSE0);
            fVerboseMenu->UnCheckEntry(M_DISPLAY_VERBOSE1);
            fVerboseMenu->UnCheckEntry(M_DISPLAY_VERBOSE2);
            sett->verbose = 0;
            break;
        }
        case 1: {
            fVerboseMenu->CheckEntry(M_DISPLAY_VERBOSE1);
            fVerboseMenu->UnCheckEntry(M_DISPLAY_VERBOSE0);
            fVerboseMenu->UnCheckEntry(M_DISPLAY_VERBOSE2);
            sett->verbose = 1;
            break;
        }
        case 2: {
            fVerboseMenu->CheckEntry(M_DISPLAY_VERBOSE2);
            fVerboseMenu->UnCheckEntry(M_DISPLAY_VERBOSE0);
            fVerboseMenu->UnCheckEntry(M_DISPLAY_VERBOSE1);
            sett->verbose = 2;
            break;
        }
    }
}

void ShapeFrame::HandleMenu(Int_t id)
{
    // Handle menu items.
    
    switch (id) {
            
        case M_FILE_ABOUT:
        {
            new ShapeInfo(gClient->GetRoot(), fMain, 600, 300, absPath);
            break;
        }
        case M_FILE_OPEN:
        {
            static TString dir(".");
            
            fi.fFileTypes = filetypes;
            fi.fIniDir    = StrDup(dir);
            
            new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
            if (fi.fFilename) {
                //make sure the previous settings file will not be overwritten
                sett->settFileName="";
                
                mname = fi.fFilename;
                //store absolute pathname
                sett->SetFileName(mname);
                //convert to filename, only (for the display)
                mname = mname.substr(mname.find_last_of("\\/") + 1, mname.length());
                OB[5]->SetEnabled(false);
                OB[6]->SetEnabled(false);
                UpdateSetting(sett);

                //create Matrix object
                matrix = new ShapeMatrix(sett);
                
                //status update
                status = 1;
          
                //update Combo Menu showing the matrices
                MatrixSelector();
                
                //select first matrix in root file
                fMatrix->Select(1);
                matrix->SetMatrix(1);
                 DoDraw();
                
                //initialize histogram zoom
                histX1 = 0;
                histX2 = matrix->GetEne0() + matrix->GetESize();
                
                //create projections
                UpdateSetting(sett);
                matrix->Diag();

            }
        }
        break;
            
        case M_FILE_EXIT:
            CloseWindow();   // terminate theApp no need to use SendCloseMessage()
            break;
        
        case M_SETTING_OPEN:
        {
            static TString dir(".");
            TGFileInfo fi_sett;
            fi_sett.fFileTypes = filetypes_s;
            fi_sett.fIniDir    = StrDup(dir);
            
            new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi_sett);
            if (fi_sett.fFilename) {
                string sname = fi_sett.fFilename;
                //store absolute pathname
                sett->settFileName = sname;
                
                //convert to filename, only (for the display)
                sname = sname.substr(sname.find_last_of("\\/") + 1, sname.length());
                OB[5]->SetEnabled(false);
                OB[6]->SetEnabled(false);
                sett->ReadSettings();
                UpdateGuiSetting(sett);
                HandleVerboseMenu(sett->verbose);
                fMain->SetWindowName(sname.c_str());

                //get nme of root file containing matrix
                mname = sett->dataFileName;
                mname = mname.substr(mname.find_last_of("\\/") + 1, mname.length());
                
                //create Matrix object
                matrix = new ShapeMatrix(sett);
                
                //status update
                status = 1;
           
                //update Combo Menu showing the matrices
                int mIndex = MatrixSelector();
                if (mIndex > 0) {
                
                    fMatrix->Select(mIndex);
                    matrix->SetMatrix(mIndex);
                    DoDraw();
                }
                else {
                    MessageBox("Error", "Cannot find matrix name stored in settings file in current root file!");
                    CloseWindow();
                }
                
                //initialize histogram zoom
                histX1 = 0;
                histX2 = histX2 = matrix->GetEne0() + matrix->GetESize();
                
                //create projections
                matrix->Diag();

            }
            
        }
            break;
        case M_SETTING_SAVEAS:
        {
            UpdateSetting(sett);
            static TString dir(".");
            TGFileInfo fi_sett;
            fi_sett.fFileTypes = filetypes_s;
            fi_sett.fIniDir    = StrDup(dir);
            
            new TGFileDialog(gClient->GetRoot(), fMain, kFDSave, &fi_sett);
            if (fi_sett.fFilename) {
                //store absolute pathname
                string sname = fi_sett.fFilename;
                sett->settFileName = sname;
        
                //convert to filename, only (for the display)
                sname = sname.substr(sname.find_last_of("\\/") + 1, sname.length());
                fMain->SetWindowName(sname.c_str());
                //call save settings
                sett->SaveSettings();
            }
        }
            break;
        case M_SETTING_SAVE:
        {
            UpdateSetting(sett);
            if (sett->settFileName == "")
                HandleMenu(M_SETTING_SAVEAS);
            else
                sett->SaveSettings();
            break;
        }
        case M_SETTING_PRINT:
            sett->PrintSettings();
            break;
        
        case M_SETTING_OSLO: {
            static TString dir(".");
            TGFileInfo fi_sett;
            fi_sett.fFileTypes = filetypes_t;
            fi_sett.fIniDir    = StrDup(dir);
            
            new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi_sett);
            if (fi_sett.fFilename) {
                sett->osloFileName = fi_sett.fFilename;
                sett->doOslo = true;
                UpdateGuiSetting(sett);
            }
            break;
        }
            
        //load energy-dependent effi parameters
        case M_SETTING_EFFI: {
            static TString dir(".");
            TGFileInfo fi_sett;
            fi_sett.fFileTypes = filetypes_t;
            fi_sett.fIniDir    = StrDup(dir);
            
            new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi_sett);
            if (fi_sett.fFilename) {
                sett->effiFileName = fi_sett.fFilename;
                sett->doEffi = true;
                sett->readEffi();
                UpdateGuiSetting(sett);
            }
            break;
        }
            //load level density file
        case M_SETTING_RHO: {
            static TString dir(".");
            TGFileInfo fi_sett;
            fi_sett.fFileTypes = filetypes_t;
            fi_sett.fIniDir    = StrDup(dir);
            
            new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi_sett);
            if (fi_sett.fFilename) {
                sett->rhoFileName = fi_sett.fFilename;
                UpdateGuiSetting(sett);
            }
            break;
        }
        
        
            
        case M_SETTING_WIDTHRESET: {
            sett->doWidthCal = 0;
            sett->ResetWidth();
            OB[5]->SetEnabled(false);
            break;
        }
        
        case M_SETTING_TRAFO: {
            AlphaDialog = new ShapeDialogAlpha(gClient->GetRoot(), fMain, this, 400, 200, sett->lit_norm, sett->lit_alpha);
            break;
        }
            
            
        
        case M_DISPLAY_MAT:
            UpdateDisplay(1);
            break;
        case M_DISPLAY_DIAG:
            UpdateDisplay(2);
            break;
        case M_DISPLAY_DIAGCUBE:
            UpdateDisplay(3);
            break;
        case M_DISPLAY_PROJTOT:
             UpdateDisplay(4);
            break;
        case M_DISPLAY_PROJBIN:
             UpdateDisplay(5);
            break;
        case M_DISPLAY_GSF:
            UpdateDisplay(6);
            break;
        case M_DISPLAY_FITWIDTH:
            UpdateDisplay(7);
            break;
        case M_DISPLAY_RATIO:
            UpdateDisplay(8);
            break;
        case M_DISPLAY_VERBOSE0:
            HandleVerboseMenu(0);
            break;
        case M_DISPLAY_VERBOSE1:
            HandleVerboseMenu(1);
            break;
        case M_DISPLAY_VERBOSE2:
            HandleVerboseMenu(2);
            break;
        case M_DISPLAY_AVG: {
            sett->displayAvg = !sett->displayAvg;
            if (sett->displayAvg)
                fDisplayFile->CheckEntry(M_DISPLAY_AVG);
            else
                fDisplayFile->UnCheckEntry(M_DISPLAY_AVG);
            break;
        }
        case M_DISPLAY_SINGLE: {
            sett->displaySingle = !sett->displaySingle;
            if (sett->displaySingle)
                fDisplayFile->CheckEntry(M_DISPLAY_SINGLE);
            else
                fDisplayFile->UnCheckEntry(M_DISPLAY_SINGLE);
            break;
        }
        
        case M_DISPLAY_COLOUR: {
            sett->colour = !sett->colour;
            if (sett->colour)
                fDisplayFile->CheckEntry(M_DISPLAY_COLOUR);
            else
                fDisplayFile->UnCheckEntry(M_DISPLAY_COLOUR);
            break;
        }
        
        case M_DISPLAY_GRF: {
            sett->doGRF = !sett->doGRF;
            if (sett->doGRF)
                fDisplayFile->CheckEntry(M_DISPLAY_GRF);
            else
                fDisplayFile->UnCheckEntry(M_DISPLAY_GRF);
            ShowGraph();
            break;
        }
        case M_DISPLAY_TRAFO: {
            
            UpdateDisplay(9);
            break;
        }
        case M_DISPLAY_EFFI: {
            UpdateDisplay(10);
            break;
        }
        case M_DISPLAY_RHO: {
            UpdateDisplay(11);
            break;
        }
        case M_DISPLAY_PRINT:
            if (status > 1)
                gSF->Print();
            else
                std::cout <<"ShapeIt, first! No values to show" <<std::endl;
            break;
            
    default:
            break;
    }
}


void ShapeFrame::UpdateDisplay(int display) {
    if (status == 0)
        return;
    TCanvas *fCanvas = fEcanvas->GetCanvas();
    fCanvas->cd();
    displayMode = display;
    if (display !=3) {
        fR[2]->SetState(kButtonUp);
        fBinCombo->SetEnabled(false);
    }
    switch (display) {
        case 1:
            matrix->GetInputMatrix(mname)->Draw("colz");
            break;
        case 2:
            matrix->GetDiagEx(mname)->Draw("hist colz");
            break;
        case 3:
            matrix->GetDiagExCube(mname)->Draw("hist colz");
            break;
        case 4:
            displayHisto = matrix->GetDiag(mname);
            displayHisto->Draw("hist");
            break;
        case 5: {
            fR[2]->SetState(kButtonDown);
            fBinCombo->SetEnabled(true);
            fBinComboDraw(fBinCombo);
            BinSelect(binSelect);
            matrix->GetDiagEx(binSelect, mname)->Draw("hist");
            fBinCombo->Connect("Selected(Int_t)", "ShapeFrame", this, "BinSelect(Int_t)");
            break;
        }
        case 7: {
            
            //get histogram of peak width for level 1 and apply linear fit
            TGraph *T1 = matrix->getFitWidthGraph(0);
            T1->Fit("pol1");
            TF1 *fit1 = T1->GetFunction("pol1");
            fit1->SetLineColor(kRed);
            sett->widthCal[0][0] = fit1->GetParameter(0);
            sett->widthCal[0][1] = fit1->GetParameter(1);
            
            //same for level 2
            TGraph *T2 = matrix->getFitWidthGraph(1);
            T2->Fit("pol1");
            TF1 *fit2 = T2->GetFunction("pol1");
            fit2->SetLineColor(kBlue);
            sett->widthCal[1][0] = fit2->GetParameter(0);
            sett->widthCal[1][1] = fit2->GetParameter(1);
            
            //update status: having autofit results and width calibration
            status = 3;
            UpdateGuiSetting(sett);
            wgraph = new TMultiGraph("wplot", "Widths from Autofit");
            
            wgraph->Add(T1);
            wgraph->Add(T2);
                        
            wgraph->Draw("ap");
            TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
            leg->SetFillColor(0);
            leg->AddEntry(T1, "level 1", "lp");
            leg->AddEntry(T2, "level 2", "lp");
            leg->Draw();
            break;
        }
        case 8: {
            if (status > 1) {
                /*TGraph* T = fitGSF->getRatioGraph();
                T->SetMarkerStyle(23);
                T->SetMarkerSize(2);
                T->SetMarkerColor(2);
                T->SetTitle("Peak area 1 / Peak area 2" );
                T->Draw("AP");*/
            }
            break;
        }
        case 9: {
            AlphaChi2();
            break;
        }
        case 10: {
            sett->eGraph->SetMarkerStyle(4);
            sett->eGraph->SetMarkerColor(kRed);
            sett->eGraph->SetTitle("Energy dependend scaling factor; gamma energy [keV]; scaling factor; non");
            sett->eGraph->Draw("AC*");
            break;
        }
        case 11: {
            ShapeRho *rho = new ShapeRho(sett, matrix);
            //rho->Draw();
            //litGSF->Rho();
            //fitGSF->Rho();
            //double rho_scale = rho->Eval(2) / litGSF->GetRhoDiagram()->Eval(2000);
            //std::cout <<"Scale: " <<rho_scale <<std::endl;
            
            //TGraphErrors* scaleGraph = new TGraphErrors;
            //for (int i=0;i<fitGSF->GetRhoDiagram()->GetN();i++) {
              //  double X =fitGSF->GetRhoDiagram()->GetX()[i]/1000;
                //double Y =fitGSF->GetRhoDiagram()->GetY()[i]*rho_scale;
                //if (X < 3.5) {
                  //  X =litGSF->GetRhoDiagram()->GetX()[i]/1000;
                    //Y =litGSF->GetRhoDiagram()->GetY()[i]*rho_scale;
                //}
                //scaleGraph->SetPoint(i,X,Y);
                //litGSF->GetRhoDiagram()->GetY()[i] *= rho_scale;
            //}
            //scaleGraph->SetMarkerStyle(4);
            //scaleGraph->SetMarkerColor(kBlue);
            
            //88Kr
            //rho->Draw(0.03,-0.045,0.10);
            
            //76Ge
            //rho->Draw(0.045,-0.01,0.10);
            
            //88Kr unfolding in y
            rho->Draw(0.4,0.3,0.5);
            
            rho->Draw();
            //scaleGraph->Draw("C*");
            
            break;
        }
    }
    DrawMarker();
    fCanvas->Modified();
    fCanvas->Update();
}

void ShapeFrame::AlphaChi2()
{
    //some hard-coded stuff here....needs to be fixed
    
    double alpha_bak = sett->lit_alpha;
    double exi_bak = sett->exiEne[0];
    double const alphaRange = 1; //will calcualte chi2 values of lit values and ShapeIt results for alpha values between - allphaRange and +alphaRange
    int const nOfPoints = 40; //number of steps in the chi2 evaluation
    int const nOfExi = 15; //number of excitation energies
    double const exiStep = 100; //step size (in keV) by which the excitation energy is lowered in each iteration
    double alphaX[nOfPoints*nOfExi];
    double chi2Y[nOfPoints*nOfExi];
    double enes[nOfPoints*nOfExi];
    int pointC = 0;
    
    for (int j = 0; j < nOfExi; j++) {
        
    for (int i  = 0; i < nOfPoints; i++) {
        alphaX[pointC] = 2*i*alphaRange / nOfPoints - alphaRange;
        //update alpha in settings file
        sett->lit_alpha = alphaX[pointC];

        //apply transformation
        ShowGraph();
        //chi2Y[pointC] = lit_chi2; NEEDS TO BE CHANGED
        enes[pointC] = sett->exiEne[0];
        pointC++;
        //std::cout <<" i = " << i << " j = " << j << "pointC = " <<pointC <<"alphaX[i]" << alphaX[pointC] << "exi[0] " << sett->exiEne[0] << "lit_chi2 " <<lit_chi2 <<  std::endl;
    }
   
    sett->exiEne[0] = sett->exiEne[0] - exiStep;
    UpdateGuiSetting(sett);
    ShapeItBaby();
    }
        
    TGraph* alphaPlot = new TGraph(nOfPoints, alphaX, chi2Y);
    TGraph2D *dt = new TGraph2D(nOfPoints*nOfExi, alphaX, enes, chi2Y);

    alphaPlot->SetMarkerStyle(4);
    alphaPlot->SetMarkerColor(kRed);
    dt->Draw("colz");
    //alphaPlot->Draw("APC*");
    
    //reset values to before
    sett->exiEne[0] = exi_bak;
    sett->lit_alpha = alpha_bak;
    UpdateGuiSetting(sett);
}


void ShapeFrame::CloseWindow()
{
    // Got close message for this MainFrame. Terminates the application.
    std::cout <<"Goodbye.Thanks for using ShapeIt!" <<std::endl;
    gApplication->Terminate();
}

ShapeFrame::~ShapeFrame() {
    // Clean up used widgets: frames, buttons, layout hints
    delete fMenuFile;
    
    fMain->Cleanup();
    delete fMain;
}


