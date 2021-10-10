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
    TGTransientFrame    *fMainDialog;
    TGGroupFrame        *fGDialog;
    TGCompositeFrame    *fTransform[2];
    TGNumberEntry       *transBDialog;              //transformation B parameter field
    TGNumberEntry       *transAlphaDialog;              //transformation alpha parameter field
   
public:
    ShapeDialogAlpha(const TGWindow *p, const TGWindow *main, ShapeFrame *caller_obj, UInt_t w, UInt_t h, double b_init = 1, double alpha_init = 0, UInt_t options = kVerticalFrame);
    
    virtual ~ShapeDialogAlpha();
    void CloseWindow();
    double GetAlphaTransform();
    double GetBTransform();   
};
#endif
