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

#ifndef SHAPEDIALOGALPHA_H
#define SHAPEDIALOGALPHA_H

#include <iostream>

//#include "ShapeFrame.h"

#include <TGFrame.h>
#include <TGNumberEntry.h>
#include <TGWindow.h>
#include <RQ_OBJECT.h>

class ShapeFrame;
class ShapeDialogAlpha {
    
    RQ_OBJECT("ShapeDialogAlpha");
    
private:
    ShapeSetting              *m_sett;
    TGTransientFrame          *fMainDialog;
    TGGroupFrame              *fGDialog;
    TGNumberEntry             *transBDialog;              //transformation B parameter field
    TGNumberEntry             *transAlphaDialog;              //transformation alpha parameter field
    TGCheckButton             *MCButton;                  //checkbutton to toggle MC mode
    TGTextButton              *fStart;                     //start button to start MC simulation
    TGTextButton              *Chi2Button;                  //start search for chi2 minimum of actual gSF values
    TGNumberEntry             *transMCIter;                 //the number of MC iterations; default is 100
    TGLabel                   *tMC[4];                      //all text labels of MC
    bool                      start;                      //status of the start button


public:
    ShapeDialogAlpha(ShapeSetting* t_sett, const TGWindow *p, const TGWindow *main, ShapeFrame *caller_obj, UInt_t w, UInt_t h, double b_init = 1, double alpha_init = 0, UInt_t options = kVerticalFrame);
    
    virtual              ~ShapeDialogAlpha();
    void                 CloseWindow();
    double               GetAlphaTransform();
    double               GetBTransform();
    bool                 GetStartStatus() {return start;}
    void                 ChangeStartLabel();
    void                 HandleMCButton();
    int                  GetNrOfIter() {return transMCIter->GetNumber();}
    TGNumberEntry        *alpha[3];                  //alpha[0,1] are the lower/upper bounds for alpha during chi2 minimum search; alpha[2] is the number of steps in the search
    void                      HandleAlpha0();
    void                      HandleAlpha1();
    void                      HandleAlpha2();
};
#endif
