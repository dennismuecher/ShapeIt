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

#ifndef SHAPETALYS_H
#define SHAPETALYS_H

#include <iostream>

#include "ShapeSetting.h"


class ShapeTalys {
    
private:
    ShapeSetting *m_sett;
   
public:
    
    ShapeRho(ShapeSetting* t_setting);
    
};
#endif
