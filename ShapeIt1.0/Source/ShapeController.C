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
#include "../Include/ShapeController.h"

// see ShapeFrame::ShapeItBaby() for the pre-refactor version this was lifted from
ShapeCollector* ShapeController::RunAnalysis(ShapeSetting* sett, ShapeMatrix* matrix) {

    sett->nOfBins = sett->SizeToBin();

    std::cout << "Running ShapeIt analysis..." << std::endl;
    
    ShapeCollector* collector = new ShapeCollector(sett, matrix);
    collector->Collect();

    std::cout << "ShapeIt analysis complete." << std::endl;
    
    return collector;
}

// see ShapeFrame::MonteCarlo() for the pre-refactor version this was lifted from
ShapeCollector* ShapeController::RunMonteCarloStep(ShapeSetting* sett, ShapeMatrix* matrix) {

    ShapeCollector* collector = new ShapeCollector(sett, matrix);
    collector->Collect();

    return collector;
}

// see ShapeFrame::AlphaChi2() for the pre-refactor version this was lifted from
ShapeAlpha* ShapeController::FitAlpha(ShapeSetting* sett, ShapeCollector* collector) {

    ShapeAlpha* frameAlpha = new ShapeAlpha(sett, collector);
    frameAlpha->Chi2Loop();

    if (sett->verbose)
        std::cout << "Minimum chi2 value of " << frameAlpha->getMinChi2()
                   << " found for alpha = " << frameAlpha->getMinAlpha() << std::endl;

    return frameAlpha;
}
