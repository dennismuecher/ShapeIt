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
#include <string>
#include <iostream>

#include <TSystem.h>

#include "ShapeTalys.C"
#include "ShapeFrame.C"
#include "ShapeSetting.C"
#include "ShapeMatrix.C"
#include "ShapeGSF.C"
#include "ShapeRho.C"
#include "ShapeCollector.C"
#include "ShapeAlpha.C"
//#include "ShapeMultiGraph.C"

void ShapeIt() {
  static const std::string path = gSystem->pwd();
	ShapeFrame* shapeIt = new ShapeFrame(gClient->GetRoot(),700,500, path);
    //shapeIt->OpenSettingFile("../Analysis/112Cd/Cd112settings.dat");
    
}
