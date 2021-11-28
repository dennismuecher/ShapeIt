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


#include "../Include/ShapeFrame.h"
#include "ShapeDialogAlpha.C"
#include "ShapeInfo.C"
#include "ShapeMultiGraph.C"

ShapeFrame::ShapeFrame(const TGWindow *p,UInt_t w,UInt_t h, const std::string path) {


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
    
    //Level Energy Settings
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
    
    //excitation energy selection
    exi[0] = new TGNumberEntry(fEnergy[2], 2500, 9,5, TGNumberFormat::kNESInteger,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,0, 99999);
    exi[0]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[2]->AddFrame(exi[0], fL2);
    TGLabel *l3 = new TGLabel(fEnergy[2], "< Excitation <");
    fEnergy[2]->AddFrame(l3, fL2);
    exi[1] = new TGNumberEntry(fEnergy[2], 7000, 9,6, TGNumberFormat::kNESInteger,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,0, 99999);
    exi[1]->Connect("ValueSet(Long_t)", "ShapeFrame", this, "DoNumberEntry()");
    fEnergy[2]->AddFrame(exi[1], fL2);
    fG[1]->AddFrame(fEnergy[2], fL2);
    
    //intgeration bin settings
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
    //OB[6] = new TGCheckButton(fG[2], new TGHotString("Use Efficiency Calibration from File"), 25);
    //OB[6]->SetEnabled(false);
    for (int i = 0 ; i < 6; i++) {
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
    if (sett->doMC)
        fPicture->Connect("Clicked()", "ShapeFrame", this, "MonteCarlo()");
    else
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
    
    if (sett->displayAvg)
        fDisplayFile->CheckEntry(M_DISPLAY_AVG);
    else
        fDisplayFile->UnCheckEntry(M_DISPLAY_AVG);
    
    if (sett->displaySingle)
        fDisplayFile->CheckEntry(M_DISPLAY_SINGLE);
    else
        fDisplayFile->UnCheckEntry(M_DISPLAY_SINGLE);
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
    fSettingsFile->AddEntry("L&oad Literature Values gSF ...", M_SETTING_OSLO);
    fSettingsFile->AddEntry("Load Literature Values &rho ...", M_SETTING_RHO);
    fSettingsFile->AddEntry("L&oad Efficiency Calibration ...", M_SETTING_EFFI);
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
    fDisplayFile->AddSeparator();
    fDisplayFile->AddEntry("Show gSF &average", M_DISPLAY_AVG);
    fDisplayFile->AddEntry("Show gSF single data", M_DISPLAY_SINGLE);
    fDisplayFile->CheckEntry(M_DISPLAY_SINGLE);
    fDisplayFile->AddEntry("Use &Two-Colour Display", M_DISPLAY_COLOUR);
    fDisplayFile->CheckEntry(M_DISPLAY_COLOUR);
    
    fDisplayFile->AddSeparator();
    fDisplayFile->AddEntry("&Print Table of gSF results", M_DISPLAY_PRINT);
    
    fDisplayFile->AddSeparator();
    fDisplayFile->AddEntry("&Draw level density", M_DISPLAY_RHO);
    
    fDisplayFile->DisableEntry(M_DISPLAY_FITWIDTH);
    
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
    //new projection selected? reset y range
    if (binSelect !=sbin) {
        binSelect = sbin;
        histY1 = 0;
        histY2 = 0;
    }
    
    diagHisto = matrix->GetDiagEx(sbin, mname);
    
    //when re-drawing, use the same X and Y ranges as before
    diagHisto->GetXaxis()->SetRangeUser(histX1, histX2);
    if (histY1 !=0 && histY2 !=0) {
        if (gPad->GetLogy())
            diagHisto->GetYaxis()->SetRangeUser(TMath::Power(10,histY1), TMath::Power(10,histY2));
        else
            diagHisto->GetYaxis()->SetRangeUser(histY1, histY2);

    }
    diagHisto->GetXaxis()->SetNdivisions(5,5,0);
    diagHisto->GetYaxis()->SetNdivisions(5,5,0);
    diagHisto->Draw();
    displayMode = 5;
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

void ShapeFrame::MessageBox(std::string title, std::string message)
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
        //draw verticl lines; there is a bug (feature?) when using log-y scale, discussed here:
        // https://root-forum.cern.ch/t/getuymax-has-problems-with-log-scale/10511/4
        //hence the lower and upper y coordinates are calculated differently for lin and log scale
        double y1 = gPad->GetUymin();
        double y2 = gPad->GetUymax();
        
        if (gPad->GetLogy() ) {
            y1 = TMath::Power(10,y1);
            y2 = TMath::Power(10,y2);
        }
        
        for (int i = 0; i < 4; i++) {
            
            l[i] = new TLine(sett->levEne[i], y1,sett->levEne[i],y2);
            l[i]->SetLineColor(kRed);
            l[i]->SetLineWidth(2);
            l[i]->SetLineColorAlpha(kRed, 0.45);
            l[i]->Draw();
        }
        bgBox[0] = new TBox(sett->bgEne[0][0],y1,sett->bgEne[0][1],y2);
        bgBox[1] = new TBox(sett->bgEne[0][2],y1,sett->bgEne[0][3],y2);
        bgBox[2] = new TBox(sett->bgEne[1][0],y1,sett->bgEne[1][1],y2);
        bgBox[3] = new TBox(sett->bgEne[1][2],y1,sett->bgEne[1][3],y2);
        
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

//called whenever the mouse is moved or clicked in the main drawing canvas; used to draw the markers and receive their values
void ShapeFrame::HandleMyCanvas(Int_t a,Int_t b,Int_t c,TObject* obj) {
    
    if (status == 0 )
           return;
       if (!(displayMode == 4 || displayMode == 5))
           return;
    
     TCanvas *fCanvas = fEcanvas->GetCanvas();
    
    //if y-axis has changed (unzoom from context menue, for example) re-draw markers
    if (histY1  != gPad->GetUymin() || histY2  != gPad->GetUymax() ) {
        histY1 = gPad->GetUymin();
        histY2 = gPad->GetUymax();
        DrawMarker();
        return;
    }

    if ( (a == kButton1Up || a==kButton2Up || a == kButton3Up )) {
        //check if any marker was moved
        for (int i = 0; i < 4; i++) {
            sett->levEne[i] = l[i]->GetX1();
        }
        
        sett->bgEne[0][0] = bgBox[0]->GetX1(); sett->bgEne[0][1] = bgBox[0]->GetX2();
        sett->bgEne[0][2] = bgBox[1]->GetX1(); sett->bgEne[0][3] = bgBox[1]->GetX2();
        sett->bgEne[1][0] = bgBox[2]->GetX1(); sett->bgEne[1][1] = bgBox[2]->GetX2();
        sett->bgEne[1][2] = bgBox[3]->GetX1(); sett->bgEne[1][3] = bgBox[3]->GetX2();
        UpdateGuiSetting(sett);

        DrawMarker();

        //if autofit active and displayed, update fit
        if (sett->mode == 2) {
            histX1 = gPad->GetUxmin();
            histX2 = gPad->GetUxmax();
            BinSelect(fBinCombo->GetSelected());

        }
    }
}

void ShapeFrame::DoDraw() {
    if (fMatrix->IsEnabled()) {
        matrix->GetInputMatrix(mname)->Draw("colz");
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

//runs a Monte Carlo simulation of the data analysis
void ShapeFrame::MonteCarlo() {
    
    //check if matrix is loaded
    if (status == 0) {
            MessageBox("Error", "No Matrix loaded!");
        return;
    }
    
    //set status of button in Transformation Dialog window
    AlphaDialog->ChangeStartLabel();
    
    //if the start button was pressed, run MC simulation
    if (sett->doMC) {
        //delete previous graph of alpha values and recreate with up-to-date bounds
        if (mcSlopeGraph != NULL)
            mcSlopeGraph->Delete();
        
        mcSlopeGraph = new TH1F("MC slope results", "MC slope results", 400, sett->alphaLimit[0], sett->alphaLimit[1] );
        
        //clean up matrix
        matrix->Reset();
        
        //set current values of excitation energies for matrix
        matrix->SetEne0( sett->exiEne[0] );
        matrix->SetEne1( sett->exiEne[1] );
        
        UpdateDisplay(6);
        //status update: will have values for gSF
        status = 2;
        
        for (int i = 0; i < AlphaDialog->GetNrOfIter(); i++) {
            
            //start-stop button pressed in the Transformation Dialog window?
            if (!AlphaDialog->GetStartStatus())
               return;
            gSFCollMC = new ShapeCollector(sett, matrix);
            gSFCollMC->Collect();
            
            //display results
            if (i%10 ==0)
                ShowGraph();
            
            //calcualte chi2 value
            mcSlopeGraph->Fill(AlphaChi2());
            
            // make sure to update display
            gSystem->ProcessEvents(); gSystem->ProcessEvents();
            gSystem->Sleep(10);
            
            if (i%10 ==0)
                std::cout <<"MC simulation: " <<i <<std::endl;
        }
    }
    
    // at this point the MC simulation is done, so set the MC button, accordingly
    sett->doMC = false;
    
    if (AlphaDialog->GetStartStatus() )
        AlphaDialog->ChangeStartLabel();

    //draw mcSlopeGraph
    TCanvas *fCanvas = fEcanvas->GetCanvas();
    fCanvas->Clear();
    if (mcSlopeGraph)
        mcSlopeGraph->Draw();
    fCanvas->Modified();
    fCanvas->Update();
}

void ShapeFrame::ShapeItBaby() {
    
    //check if matrix is loaded
    if (status == 0) {
            MessageBox("Error", "No Matrix loaded!");
        return;
    }
    
    //setting display mode
    UpdateDisplay(6);
    
    //update settings accordings to the GUI settings
    UpdateSetting(sett);
    
    sett->nOfBins = sett->SizeToBin();
    
    //this call triggers the calculation of gSF values according to the actual settings stored int he settings file
    gSFColl = new ShapeCollector(sett, matrix);
    gSFColl->Collect();
    
    //status update: will have values for gSF
    status = 2;
    
    //enable Menu entry showing fit results for peak width and peak ratios
    if (sett->mode == 2)
        fDisplayFile->EnableEntry(M_DISPLAY_FITWIDTH);
        
    //display results
    ShowGraph();
}

//provides a TPaveText window with infos on the chi2 minimum for the ShapeIt Display
TPaveText* ShapeFrame::getPaveTextShape() {
    
    //can't calculate chi2 if no literature data aare loaded
    if (!sett->doOslo)
        return NULL;
    
    TPaveText *t=new TPaveText(0.8,0.85,0.95,0.95,"brNDC");
    t->SetTextSize(0.025);
    t->SetTextAlign(13);
    t->SetFillColor(10);
    t->SetTextColor(61);
    t->AddText(Form("slope #alpha: %4.2f ",sett->lit_alpha));

    t->AddText(Form("#chi^{2} value: %4.2f ",gSFColl->getChi2()));
    return t;
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
    if (status == 0) {
            MessageBox("Error", "No Matrix loaded!");
        return;
    }
    
    if (status == 1) {
        MessageBox("Error", "Press ShapeIt button, first! No gamma ray strength values to show");
        return;
    }

    TCanvas *fCanvas = fEcanvas->GetCanvas();
    fCanvas->Clear();
    //draw the gSF collector object
    if (sett->doMC)
        gSFCollMC->Draw();
    else {
        gSFColl->Draw();
        //display Chi2 minimum, in case literature values are shown
        if (getPaveTextShape() != NULL)
            getPaveTextShape()->Draw();
    }
    
    //update canvas
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

//takes care of setting the verbose level and displaying the verbose cascade menu correctly
void ShapeFrame::HandleVerboseMenu(int vLevel) {
    
    sett->verbose = vLevel;
    fVerboseMenu->UnCheckEntry(M_DISPLAY_VERBOSE0);
    fVerboseMenu->UnCheckEntry(M_DISPLAY_VERBOSE1);
    fVerboseMenu->UnCheckEntry(M_DISPLAY_VERBOSE2);
    
    switch (vLevel) {
        case 0: {
            fVerboseMenu->CheckEntry(M_DISPLAY_VERBOSE0);
            break;
        }
        case 1: {
            fVerboseMenu->CheckEntry(M_DISPLAY_VERBOSE1);
            sett->verbose = 1;
            break;
        }
        case 2: {
            fVerboseMenu->CheckEntry(M_DISPLAY_VERBOSE2);
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
            //new ShapeInfo(gClient->GetRoot(), fMain, 600, 300, absPath);
            MessageBox("Welcome to ShapeIt!","ShapeIt Version 1.1 \n © 2021 Dennis Muecher \n Questions? Comments? dmuecher@uoguelph.ca \n This program is free software: you can redistribute it and/or modify it under the \n terms of the GNU General Public License as published by the Free Software Foundation, \n either version 3 of the License, or (at your option) any later version. \n You should have received a copy of the GNU General Public License \n along with this program. If not, see  http://www.gnu.org/licenses/.");
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
                //OB[6]->SetEnabled(false);
                UpdateSetting(sett);

                //create Matrix object
                matrix = new ShapeMatrix(sett);
                
                //status update
                status = 1;
                displayMode = 1;
                //update Combo Menu showing the matrices
                MatrixSelector();
                
                //select first matrix in root file
                fMatrix->Select(1);
                matrix->SetMatrix(1);
                 DoDraw();
                
                //initialize histogram zoom
                histX1 = 0;
                histX2 = matrix->GetEne0() + matrix->GetESize();
                histY1 = 0; histY2 = 0;
                
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
                std::string sname = fi_sett.fFilename;
                //store absolute pathname
                sett->settFileName = sname;
                
                //convert to filename, only (for the display)
                sname = sname.substr(sname.find_last_of("\\/") + 1, sname.length());
                OB[5]->SetEnabled(false);
                //OB[6]->SetEnabled(false);
                displayMode = 1;
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
                histY1 = 0; histY2 = 0;
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
                std::string sname = fi_sett.fFilename;

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
            AlphaDialog = new ShapeDialogAlpha(sett, gClient->GetRoot(), fMain, this, 400, 200, sett->lit_norm, sett->lit_alpha);
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
                gSFColl->Print();
            else
                MessageBox("Error", "Press ShapeIt button, first! No gamma ray strength values to show");
            break;
            
    default:
            break;
    }
}


void ShapeFrame::UpdateDisplay(int display) {
    if (status == 0)
        return;
   
    displayMode = display;
    if (display ==6)
        return;
    TCanvas *fCanvas = fEcanvas->GetCanvas();
       fCanvas->cd();
       fCanvas->Clear();
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
        
        case 10: {
            sett->eGraph->SetMarkerStyle(4);
            sett->eGraph->SetMarkerColor(kRed);
            sett->eGraph->SetTitle("Energy dependend scaling factor; gamma energy [keV]; scaling factor; non");
            sett->eGraph->Draw("AC*");
            break;
        }
        case 11: {
            
            TMultiGraph* m_graph = new TMultiGraph();
            ShapeRho *rho = new ShapeRho(sett);
            double alpha_error = 0.05; //THIS SHOULD NOT BE HARD CODED
            TGraphAsymmErrors* rhoTrafo = rho->rhoTrafoGraph(sett->lit_alpha,sett->lit_alpha - alpha_error, sett->lit_alpha + alpha_error );
            
            ShapeMultiGraph *sMultiGraph = new ShapeMultiGraph();
            
            rhoTrafo->SetMarkerColor(kBlack);
            rhoTrafo->SetLineColor(kBlack);
            
            m_graph->Add(rhoTrafo,"APE");
            
            //create ShapeTalys object using ld5 model for Kr88 data; this uses the recommended values for ptable and ctable
/*
            ShapeTalys* ld4 = new ShapeTalys("../Talys/140Ba/Ba140_ld4.out",rhoTrafo,false,0,0.139,0.0);
            ShapeTalys* ld5 = new ShapeTalys("../Talys/140Ba/Ba140_ld5.out",rhoTrafo,true,0,0.712,0.0);
            ShapeTalys* ld6 = new ShapeTalys("../Talys/140Ba/Ba140_ld6.out",rhoTrafo,true,0,0.4,0.0);
            
            
            ShapeTalys* ld4 = new ShapeTalys("../Talys/88Kr/Kr88_ld4.out",rhoTrafo, false,0,0.06, 0.0);
            ShapeTalys* ld5 = new ShapeTalys("../Talys/88Kr/Kr88_ld5.out",rhoTrafo, true,0,0.78, 0.0);
            ShapeTalys* ld6 = new ShapeTalys("../Talys/88Kr/Kr88_ld6.out",rhoTrafo, true,0,0.68, 0.0);*/
            
            ShapeTalys* ld4 = new ShapeTalys("../../Talys/143Ba/Ba143_ld4.out",rhoTrafo,false,0,-0.06, 0.0);
            ShapeTalys* ld5 = new ShapeTalys("../../Talys/143Ba/Ba143_ld5.out",rhoTrafo,true,0,-0.489, 0.0);
            ShapeTalys* ld6 = new ShapeTalys("../../Talys/140Ba/Ba140_ld6.out",rhoTrafo,true,0,-0.694, 0.0);
            
            
            TGraph* ld4Graph = ld4->getDenPartialGraphTrans();
            ld4Graph->SetTitle("ld4");
            ld4Graph->SetLineColor(kBlack);
            ld4Graph->SetLineWidth(3);
            sMultiGraph->Add((TGraph*)ld4Graph->Clone(),"L");
            
            TGraph* ld5Graph = ld5->getDenPartialGraphTrans();
            ld5Graph->SetTitle("ld5");
            ld5Graph->SetLineColor(61);
            ld5Graph->SetLineWidth(3);
            sMultiGraph->Add((TGraph*)ld5Graph->Clone(),"L");
            
            TGraph* ld6Graph = ld6->getDenPartialGraphTrans();
            ld6Graph->SetTitle("ld6");
            ld6Graph->SetLineColor(41);
            ld6Graph->SetFillColor(kRed);
            ld6Graph->SetLineWidth(3);
            sMultiGraph->Add((TGraph*)ld6Graph->Clone(),"L");
            
            //create the band of theoretical values
            sMultiGraph->doFill(1,3);
            TGraphErrors* theoBand = (TGraphErrors*)sMultiGraph->fillGraph->Clone();
            theoBand->SetFillColor(4);
            theoBand->SetFillStyle(3010);
            //add band to the multigraph
            m_graph->Add(theoBand,"3");
            
            
            //run chi2 minimization to fit experimental partial level density with theoretical model
            ld4->Chi2PartialLoop(2.0, 10.0);
            //set ptable and ctable to this minimum
            ld4->BestFitPartial();
            //add the best fit to the multigraph
            TGraph* bestPartialGraph4 = ld4->getDenPartialGraphTrans();
            bestPartialGraph4->SetTitle("best fit data");
            bestPartialGraph4->SetLineColor(kRed);
            bestPartialGraph4->SetLineWidth(3);
            sMultiGraph->Add((TGraph*)bestPartialGraph4->Clone());
            
            //run chi2 minimization to fit experimental partial level density with theoretical model
            ld5->Chi2PartialLoop(2.0, 10.0);
            //set ptable and ctable to this minimum
            ld5->BestFitPartial();
            //add the best fit to the multigraph
            TGraph* bestPartialGraph = ld5->getDenPartialGraphTrans();
            bestPartialGraph->SetTitle("best fit data");
            bestPartialGraph->SetLineColor(51);
            bestPartialGraph->SetLineWidth(3);
            sMultiGraph->Add((TGraph*)bestPartialGraph->Clone());
            
            
            //run chi2 minimization to fit experimental partial level density with theoretical model
            ld6->Chi2PartialLoop(2.0, 10.0);
            //set ptable and ctable to this minimum
            ld6->BestFitPartial();
            //add the best fit to the multigraph
            TGraph* bestPartialGraph6 = ld6->getDenPartialGraphTrans();
            bestPartialGraph6->SetTitle("best fit data ld6");
            bestPartialGraph6->SetLineColor(51);
            bestPartialGraph6->SetLineWidth(3);
            sMultiGraph->Add((TGraph*)bestPartialGraph6->Clone());
            
            //create the band of experimental values
            sMultiGraph->doFill(4,6);
            TGraphErrors* expBand = (TGraphErrors*)sMultiGraph->fillGraph->Clone();
            expBand->SetFillColorAlpha(kRed, 0.8);
            expBand->SetFillStyle(3010);
            
            //add band to the multigraph
            m_graph->Add(expBand,"3");
            
            m_graph->GetYaxis()->SetRangeUser(0.01,1000);
            m_graph->Draw("apl");
            //sMultiGraph->Draw("apl");
            m_graph->GetYaxis()->SetTitleOffset(1.4);
            m_graph->GetYaxis()->SetTitle("Level density #rho (E) (MeV^{-1})");
            m_graph->GetXaxis()->SetTitle("Energy (MeV)");
            m_graph->SetTitle("Level Density Kr^{88}");
            //m_graph->Draw("apl");
            
            //draw discete levels into same canvas
            ld5->discreteHist->Draw("same hist");
            
            
            /* ShapeTalys* ld1 = new ShapeTalys("../Talys/140Ba/Ba140_ld1.out", rhoTrafo, false, 1, 0.8, 0.4);
             TGraph* ld1Graph = ld1->getDenPartialGraph();
             ld1Graph->SetTitle("ld1");
             ld1Graph->SetLineColor(51);
             ld1Graph->SetLineWidth(2);
             ld1->Chi2PartialLoop();
             m_graph->Add(ld1Graph,"L");
             
             ShapeTalys* ld2 = new ShapeTalys("../Talys/140Ba/Ba140_ld2.out",false,1,0.0);
             TGraph* ld2Graph = ld2->getDenPartialGraph();
             ld2Graph->SetTitle("ld2");
             ld2Graph->SetLineColor(61);
             ld2Graph->SetLineWidth(2);
             m_graph->Add(ld2Graph,"L");
             
             ShapeTalys* ld3 = new ShapeTalys("../Talys/140Ba/Ba140_ld3.out",false,1,0.0);
             TGraph* ld3Graph = ld3->getDenPartialGraph();
             ld3Graph->SetTitle("ld3");
             ld3Graph->SetLineColor(71);
             ld3Graph->SetLineWidth(2);
             m_graph->Add(ld3Graph,"L");
             
             ShapeTalys* ld4 = new ShapeTalys("../Talys/140Ba/Ba140_ld4.out",false,0,0.139);
             TGraph* ld4Graph = ld4->getDenPartialGraph();
             ld4Graph->SetTitle("ld4");
             ld4Graph->SetLineColor(81);
             ld4Graph->SetLineWidth(2);
             m_graph->Add(ld4Graph,"L");
             
             ShapeTalys* ld5 = new ShapeTalys("../Talys/140Ba/Ba140_ld5.out",true,0,0.712);
             TGraph* ld5Graph = ld5->getDenPartialGraph();
             ld5Graph->SetTitle("ld5");
             ld5Graph->SetLineColor(91);
             ld5Graph->SetLineWidth(2);
             m_graph->Add(ld5Graph,"L");
             
             ShapeTalys* ld6 = new ShapeTalys("../Talys/140Ba/Ba140_ld6.out",true,0,0.4);
             TGraph* ld6Graph = ld6->getDenPartialGraph();
             ld6Graph->SetTitle("ld6");
             ld6Graph->SetLineColor(kRed);
             ld6Graph->SetLineWidth(2);
             m_graph->Add(ld6Graph,"L");
             
             ShapeTalys* ld1 = new ShapeTalys("../../Talys/143Ba/Ba143_ld1.out",false,1,0.0);
             TGraph* ld1Graph = ld1->getDenPartialGraph();
             ld1Graph->SetTitle("ld1");
             ld1Graph->SetLineColor(51);
             ld1Graph->SetLineWidth(2);
             m_graph->Add(ld1Graph,"L");
             
             ShapeTalys* ld2 = new ShapeTalys("../../Talys/143Ba/Ba143_ld2.out",false,1,0.0);
             TGraph* ld2Graph = ld2->getDenPartialGraph();
             ld2Graph->SetTitle("ld2");
             ld2Graph->SetLineColor(61);
             ld2Graph->SetLineWidth(2);
             m_graph->Add(ld2Graph,"L");
             
             ShapeTalys* ld3 = new ShapeTalys("../../Talys/143Ba/Ba143_ld3.out",false,1,0.0);
             TGraph* ld3Graph = ld3->getDenPartialGraph();
             ld3Graph->SetTitle("ld3");
             ld3Graph->SetLineColor(71);
             ld3Graph->SetLineWidth(2);
             m_graph->Add(ld3Graph,"L");
             
             ShapeTalys* ld4 = new ShapeTalys("../../Talys/143Ba/Ba143_ld4.out",false,0,-0.06);
             TGraph* ld4Graph = ld4->getDenPartialGraph();
             ld4Graph->SetTitle("ld4");
             ld4Graph->SetLineColor(81);
             ld4Graph->SetLineWidth(2);
             m_graph->Add(ld4Graph,"L");
             
             ShapeTalys* ld5 = new ShapeTalys("../../Talys/143Ba/Ba143_ld5.out",true,0,-0.489);
             TGraph* ld5Graph = ld5->getDenPartialGraph();
             ld5Graph->SetTitle("ld5");
             ld5Graph->SetLineColor(91);
             ld5Graph->SetLineWidth(2);
             m_graph->Add(ld5Graph,"L");
             
             ShapeTalys* ld6 = new ShapeTalys("../../Talys/140Ba/Ba140_ld6.out",true,0,-0.694);
             TGraph* ld6Graph = ld6->getDenPartialGraph();
             ld6Graph->SetTitle("ld6");
             ld6Graph->SetLineColor(kRed);
             ld6Graph->SetLineWidth(2);
             m_graph->Add(ld6Graph,"L");*/
            
            /*ShapeTalys* ld1 = new ShapeTalys("../Talys/76Ge/Ge76_ld1.out",false,1,0.0);
             TGraph* ld1Graph = ld1->getDenPartialGraph();
             ld1Graph->SetTitle("ld1");
             ld1Graph->SetLineColor(51);
             ld1Graph->SetLineWidth(2);
             m_graph->Add(ld1Graph,"L");
             
             ShapeTalys* ld2 = new ShapeTalys("../Talys/76Ge/Ge76_ld2.out",false,1,0.0);
             TGraph* ld2Graph = ld2->getDenPartialGraph();
             ld2Graph->SetTitle("ld2");
             ld2Graph->SetLineColor(61);
             ld2Graph->SetLineWidth(2);
             m_graph->Add(ld2Graph,"L");
             
             ShapeTalys* ld3 = new ShapeTalys("../Talys/76Ge/Ge76_ld3.out",false,1,0.0);
             TGraph* ld3Graph = ld3->getDenPartialGraph();
             ld3Graph->SetTitle("ld3");
             ld3Graph->SetLineColor(71);
             ld3Graph->SetLineWidth(2);
             m_graph->Add(ld3Graph,"L");
             
             ShapeTalys* ld4 = new ShapeTalys("../Talys/76Ge/Ge76_ld4.out",false,0,0.505);
             TGraph* ld4Graph = ld4->getDenPartialGraph();
             ld4Graph->SetTitle("ld4");
             ld4Graph->SetLineColor(81);
             ld4Graph->SetLineWidth(2);
             m_graph->Add(ld4Graph,"L");
             
             ShapeTalys* ld5 = new ShapeTalys("../Talys/76Ge/Ge76_ld5.out",true,0,0.889);
             TGraph* ld5Graph = ld5->getDenPartialGraph();
             ld5Graph->SetTitle("ld5");
             ld5Graph->SetLineColor(91);
             ld5Graph->SetLineWidth(2);
             m_graph->Add(ld5Graph,"L");
             
             ShapeTalys* ld6 = new ShapeTalys("../Talys/76Ge/Ge76_ld6.out",true,0,0.327);
             TGraph* ld6Graph = ld6->getDenPartialGraph();
             ld6Graph->SetTitle("ld6");
             ld6Graph->SetLineColor(kRed);
             ld6Graph->SetLineWidth(2);
             m_graph->Add(ld6Graph,"L"); */
            
            /*ShapeTalys* ld1 = new ShapeTalys("../Talys/88Kr/Kr88_ld1.out",false,1,0.0);
             TGraph* ld1Graph = ld1->getDenPartialGraph();
             ld1Graph->SetTitle("ld1");
             ld1Graph->SetLineColor(51);
             ld1Graph->SetLineWidth(2);
             m_graph->Add(ld1Graph,"L");
             
             ShapeTalys* ld2 = new ShapeTalys("../Talys/88Kr/Kr88_ld2.out",false,1,0.0);
             TGraph* ld2Graph = ld2->getDenPartialGraph();
             ld2Graph->SetTitle("ld2");
             ld2Graph->SetLineColor(61);
             ld2Graph->SetLineWidth(2);
             m_graph->Add(ld2Graph,"L");
             
             ShapeTalys* ld3 = new ShapeTalys("../Talys/88Kr/Kr88_ld3.out",false,1,0.0);
             TGraph* ld3Graph = ld3->getDenPartialGraph();
             ld3Graph->SetTitle("ld3");
             ld3Graph->SetLineColor(71);
             ld3Graph->SetLineWidth(2);
             m_graph->Add(ld3Graph,"L");
             
             ShapeTalys* ld4 = new ShapeTalys("../Talys/88Kr/Kr88_ld4.out",false,0,0.06);
             TGraph* ld4Graph = ld4->getDenPartialGraph();
             ld4Graph->SetTitle("ld4");
             ld4Graph->SetLineColor(81);
             ld4Graph->SetLineWidth(2);
             m_graph->Add(ld4Graph,"L");
             
             ShapeTalys* ld5 = new ShapeTalys("../Talys/88Kr/Kr88_ld5.out",true,0,0.78);
             TGraph* ld5Graph = ld5->getDenPartialGraph();
             ld5Graph->SetTitle("ld5");
             ld5Graph->SetLineColor(91);
             ld5Graph->SetLineWidth(2);
             m_graph->Add(ld5Graph,"L");
             
             ShapeTalys* ld6 = new ShapeTalys("../Talys/88Kr/Kr88_ld6.out",true,0,0.68);
             TGraph* ld6Graph = ld6->getDenPartialGraph();
             ld6Graph->SetTitle("ld6");
             ld6Graph->SetLineColor(kRed);
             ld6Graph->SetLineWidth(2);
             m_graph->Add(ld6Graph,"L");*/
            
            fCanvas->SetLogy();
            fCanvas->BuildLegend();
            
            
        }
    }
    fCanvas->Modified();
    fCanvas->Update();
}

double ShapeFrame::AlphaChi2() {
    
    //check status
    if (status == 0) {
            MessageBox("Error", "No Matrix loaded!");
        return 0;
    }
    if (status < 2) {
        MessageBox("Error", "Press ShapeIt button, first! No gamma ray strength values to show");
        return 0;
    }
        
    ShapeAlpha* frameAlpha;

    if (sett->doMC)
        frameAlpha = new ShapeAlpha(sett,gSFCollMC);
    else
        frameAlpha = new ShapeAlpha(sett,gSFColl);

    //run search for chi2 minimum
    frameAlpha->Chi2Loop();
    
    TGraph *test = frameAlpha->getChi2Graph();
    
    TCanvas *fCanvas = fEcanvas->GetCanvas();
    TPaveText* t = frameAlpha->getPaveTextChi2();
    
    if (!sett->doMC) {
        test->Draw("APC*");
        t->Draw();
        fCanvas->Modified();
        fCanvas->Update();
    }
    
    if (sett->verbose)
        std::cout <<"Minimum chi2 value of "<< frameAlpha->getMinChi2() << " found for alpha = " << frameAlpha->getMinAlpha() <<std::endl;
    
    
    return (frameAlpha->getMinAlpha()) ;
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


