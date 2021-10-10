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
#ifndef SHAPEINFO_H
#define SHAPEINFO_H

#include <string>
#include <iostream>

#include <TGFrame.h>
#include <TGWindow.h>
#include <TGIcon.h>
#include <TGPicture.h>

class ShapeInfo {
    
private:
    TGTransientFrame  *fMain;
    TGVerticalFrame   *fFrame1;
    TGLayoutHints     *fC;
    TGPicture         *ipic;
    TGIcon            *icon;
    
public:
    ShapeInfo(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h, std::string path);
    virtual ~ShapeInfo();
    
    // slots
    void CloseWindow();
    void DoClose();
};
#endif
