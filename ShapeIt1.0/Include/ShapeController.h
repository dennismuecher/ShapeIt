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

#ifndef SHAPECONTROLLER_H
#define SHAPECONTROLLER_H

#include "ShapeSetting.h"
#include "ShapeMatrix.h"
#include "ShapeCollector.h"
#include "ShapeAlpha.h"

// ShapeController holds the "what to compute" logic that used to live inside
// ShapeFrame. It knows nothing about GUI widgets, canvases, or message boxes:
// every method takes the settings/matrix/collector it needs as arguments and
// returns a result object. ShapeFrame is responsible for reading widget
// values into ShapeSetting *before* calling these, and for drawing/displaying
// the result afterwards.
//
// Methods are static for now: there is no per-instance state to carry between
// calls, and this keeps the ShapeFrame integration to a one-line change per
// call site (see ShapeFrame::ShapeItBaby / ShapeFrame::AlphaChi2 for the
// pattern). If shared state between steps becomes necessary later (e.g. to
// avoid re-passing sett/matrix everywhere), this can be made a normal
// instance class without changing the call sites much.
class ShapeController {
public:

    // Runs one full Shape Method analysis (sliding window + bin-size
    // variation as configured in sett) for the given matrix and settings.
    // Equivalent to the old ShapeFrame::ShapeItBaby(), minus the GUI parts
    // (status checks, display mode switch, enabling the fit-width menu
    // entry, and the final ShowGraph() call all stay in ShapeFrame).
    //
    // Caller takes ownership of the returned ShapeCollector.
    static ShapeCollector* RunAnalysis(ShapeSetting* sett, ShapeMatrix* matrix);

    // Runs a single Monte Carlo realization: randomizes bin size / sliding
    // window per sett, collects gSF data, and returns the resulting
    // collector. Equivalent to the per-iteration body of the old
    // ShapeFrame::MonteCarlo() loop, minus the progress display, the
    // start/stop button check, and the histogram fill (still done in
    // ShapeFrame, since GetStartStatus()/UpdateDisplay() are GUI concerns).
    //
    // Caller takes ownership of the returned ShapeCollector.
    static ShapeCollector* RunMonteCarloStep(ShapeSetting* sett, ShapeMatrix* matrix);

    // Fits the slope correction alpha by minimizing chi2 between the gSF
    // data in collector and the literature data referenced in sett.
    // Equivalent to the computational part of the old
    // ShapeFrame::AlphaChi2(); the caller is responsible for drawing
    // frameAlpha->getChi2Graph() / getPaveTextChi2() and reading
    // frameAlpha->getMinAlpha() for its own purposes.
    //
    // Caller takes ownership of the returned ShapeAlpha.
    static ShapeAlpha* FitAlpha(ShapeSetting* sett, ShapeCollector* collector);
};
#endif
