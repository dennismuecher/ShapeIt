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

#include "../Include/ShapeDialogAlpha.h"
#include "../Include/ShapeFrame.h"


ShapeDialogAlpha::ShapeDialogAlpha(ShapeSetting* t_sett, const TGWindow *p, const TGWindow *main, ShapeFrame *caller_obj, UInt_t w, UInt_t h, double b_init, double alpha_init, UInt_t options)
{
    //the settings file
    m_sett = t_sett;
    m_sett->doMC = false;
    //frame of this dialog class
    fMainDialog = new TGTransientFrame(p, main, w, h, options);
    fMainDialog->Connect("CloseWindow()", "ShapeDialogAlpha", this, "CloseWindow()");
    fMainDialog->DontCallClose(); // to avoid double deletions.
    
    // use hierarchical cleaning
    fMainDialog->SetCleanup(kDeepCleanup);
        
    //some layouts
    TGLayoutHints *fL = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);
    TGLayoutHints *fR = new TGLayoutHints( kLHintsLeft | kLHintsExpandY,
                            2, 5, 5, 2);
    TGLayoutHints *fL2 = new TGLayoutHints( kLHintsLeft | kLHintsExpandY | kLHintsExpandX,
                            2, 5, 5, 2);
    
    fGDialog = new TGGroupFrame(fMainDialog, new TGString("Transformation of Oslo Literature values"),kVerticalFrame|kRaisedFrame);
    
    //the composite frame for setting slope and scaling of the literature data
    TGCompositeFrame *fTransform[2];
    for (int i = 0; i < 2; i++)
       fTransform[i] = new TGCompositeFrame(fGDialog, 1, 1, kHorizontalFrame);
    
    //this is the scaling factor for the literature data
    TGLabel *t1 = new TGLabel(fTransform[0], "gSF scale B");
    fTransform[0]->AddFrame(t1, fR);
    transBDialog = new TGNumberEntry(fTransform[0], b_init, 9,31, TGNumberFormat::kNESReal,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,0, 99999);
    
    fTransform[0]->AddFrame(transBDialog, fR);

    TGLabel *t2 = new TGLabel(fTransform[1], " slope Alpha");
    fTransform[1]->AddFrame(t2, fR);
    transAlphaDialog = new TGNumberEntry(fTransform[1], alpha_init, 9,32, TGNumberFormat::kNESReal,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,-100, 100);
    fTransform[1]->AddFrame(transAlphaDialog, fR);
    TGLabel *t3 = new TGLabel(fTransform[1], " [1/MeV]");
    fTransform[1]->AddFrame(t3, fR);
    
    /*
    
   ****** Chi2 search Settings *****
    
    */
    
    //the group frame
    TGGroupFrame *fGChi2 = new TGGroupFrame(fMainDialog, new TGString("Chi2 Search Settings (also applies to MC mode below)"), kVerticalFrame|kRaisedFrame);
    
    TGCompositeFrame *fChi2 = new TGCompositeFrame(fGChi2, 1, 1, kHorizontalFrame);
    
    //alpha range for chi2 search (also used for Monte Carlo search)
    alpha[0] = new TGNumberEntry(fChi2, m_sett->alphaLimit[0], 9,1, TGNumberFormat::kNESReal,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,-999, 999);
    fChi2->AddFrame(alpha[0], fL);
    tMC[1] = new TGLabel(fChi2, "< Alpha <");
    fChi2->AddFrame(tMC[1], fL2);
    alpha[1] = new TGNumberEntry(fChi2, m_sett->alphaLimit[1], 9,1, TGNumberFormat::kNESReal,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,-999, 999);
    fChi2->AddFrame(alpha[1], fL);
    tMC[2] = new TGLabel(fChi2, "in");
    fChi2->AddFrame(tMC[2], fL2);
    alpha[2] = new TGNumberEntry(fChi2, m_sett->alphaIter, 9,1, TGNumberFormat::kNESInteger,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,5, 999);
    fChi2->AddFrame(alpha[2], fL);
    tMC[3] = new TGLabel(fChi2, "steps");
    fChi2->AddFrame(tMC[3], fL2);
    
    //the chi2-search button
    Chi2Button = new TGTextButton(fGChi2, "&Search Chi2 Minimum!");
    Chi2Button->SetToolTipText("During the search, the literare values are transformed using slopes within the set limits, and the chi2 in comparison to the ShapeIt gSF values is computed for each iteration. Then the slope alpha with a minimum chi2 value is computed.");
    fGChi2->AddFrame(fChi2, fL);
    fGChi2->AddFrame(Chi2Button, new TGLayoutHints(kLHintsTop,
                                                4, 3, 3, 4));
    
    /*
    
   ****** Monte Carlo Settings *****
    
    */
    
    //the group frame
    TGGroupFrame *fGMC = new TGGroupFrame(fMainDialog, new TGString("Monte Carlo Simulation of best Slopes"), kVerticalFrame|kRaisedFrame);

    //the horizontal composite frames
    TGCompositeFrame *fMC = new TGCompositeFrame(fGMC, 1, 1, kHorizontalFrame);
    
    //the number of MC interations
    tMC[0] = new TGLabel(fMC, " Number of Monte Carlo iterations");
    fMC->AddFrame(tMC[0], fR);
    transMCIter = new TGNumberEntry(fMC, 100, 9,32, TGNumberFormat::kNESInteger,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,1, 999999);
    fMC->AddFrame(transMCIter, fR);
    
    //the start-stop button
    fStart = new TGTextButton(fGMC, "&Run Monte Carlo Simulation!");
    fStart->SetToolTipText("During each iteration of the Monte Carlo simulation, values for gSF are determined based on randomized choices of the input parameters (witin the limits set by the user). Also, in each iteration specific values for the literature values are chosen within their uncertainties. In each iteration, a best slope alpha is determined, comaring the ShapeIt gSF values with the literature values. In the end, a distribution of the best alpha values is displayed.");
    start = kFALSE;
    
    //add all frames
    fGDialog->AddFrame(fTransform[0], fL);
    fGDialog->AddFrame(fTransform[1], fL);
    fGMC->AddFrame(fMC, fL);
    fGMC->AddFrame(fStart, new TGLayoutHints(kLHintsTop,
                                                4, 3, 3, 4));
    fMainDialog->AddFrame(fGDialog, fL);
    fMainDialog->AddFrame(fGChi2, fL);
    fMainDialog->AddFrame(fGMC, fL);

    //signals
    transBDialog->Connect("ValueSet(Long_t)", "ShapeFrame", caller_obj, "TransGraph()");
    transAlphaDialog->Connect("ValueSet(Long_t)", "ShapeFrame", caller_obj, "TransGraph()");
    alpha[0]->Connect("ValueSet(Long_t)", "ShapeDialogAlpha", this, "HandleAlpha0()");
    alpha[1]->Connect("ValueSet(Long_t)", "ShapeDialogAlpha", this, "HandleAlpha1()");
    fStart->Connect("Clicked()", "ShapeFrame", caller_obj, "MonteCarlo()");
    Chi2Button->Connect("Clicked()", "ShapeFrame", caller_obj, "AlphaChi2()");

    alpha[2]->Connect("ValueSet(Long_t)", "ShapeDialogAlpha", this, "HandleAlpha2()");
    

    fMainDialog->MapSubwindows();
    fMainDialog->Resize();

    // position relative to the parent's window
    fMainDialog->CenterOnParent();

    fMainDialog->SetWindowName("Transformation");

    fMainDialog->MapWindow();
    //gClient->WaitFor(fMain);    // otherwise canvas contextmenu does not work
}

void ShapeDialogAlpha::HandleAlpha2() {
    m_sett->alphaIter = alpha[2]->GetNumber();
}

void ShapeDialogAlpha::HandleAlpha0() {
    if ( alpha[0]->GetNumber() >= alpha[1]->GetNumber() )
        alpha[0]->SetNumber(alpha[0]->GetNumber() - 1);
    m_sett->alphaLimit[0] = alpha[0]->GetNumber();
}

void ShapeDialogAlpha::HandleAlpha1() {
    if ( alpha[1]->GetNumber() <= alpha[0]->GetNumber() )
        alpha[1]->SetNumber(alpha[1]->GetNumber() + 1);
    m_sett->alphaLimit[1] = alpha[1]->GetNumber();
}


void ShapeDialogAlpha::ChangeStartLabel()
{
  // Slot connected to the Clicked() signal.
  // It will toggle labels "Start" and "Stop".

  fStart->SetState(kButtonDown);
  if (!start) {
     fStart->SetText("&Stop MC!");
     start = kTRUE;
      m_sett->doMC = true;
  } else {
     fStart->SetText("&Start MC!");
     start = kFALSE;
      m_sett->doMC = false;
  }
  fStart->SetState(kButtonUp);
}

ShapeDialogAlpha::~ShapeDialogAlpha()
{
    // Delete test dialog widgets.
    
    fMainDialog->DeleteWindow();  // deletes fMain
    //sett doMC to false
    start= false;
    m_sett->doMC = false;
}

double ShapeDialogAlpha::GetAlphaTransform()
{
    return transAlphaDialog->GetNumber();
}

double ShapeDialogAlpha::GetBTransform()
{
    return transBDialog->GetNumber();
}

void ShapeDialogAlpha::CloseWindow()
{
    // Called when window is closed via the window manager.
    
    delete this;
}


