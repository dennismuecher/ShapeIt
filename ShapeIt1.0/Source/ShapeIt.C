/*************************************************************************
* Copyright (C) 2019-2020, Dennis Muecher.                               *
* All rights reserved.                                                   *
*                                                                        *
* This program is free software: you can redistribute it and/or modify   *
* it under the terms of the GNU General Public License as published by   *
* the Free Software Foundation, either version 3 of the License, or      *
* (at your option) any later version.                                    *
* You should have received a copy of the GNU General Public License      *
* along with this program. If not, see  http://www.gnu.org/licenses/.    *
*************************************************************************/

#include "ShapeFrame.C"
#include "ShapeSetting.C"
#include "ShapeMatrix.C"
#include "ShapeGSF.C"
#include "ShapeChi2.C"
#include "ShapeRho.C"

void ShapeIt() {
    static const string path = gSystem->pwd();
	ShapeFrame* test = new ShapeFrame(gClient->GetRoot(),700,500, path);
}
