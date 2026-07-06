// WebShapeIt.cxx
//
// First real increment of ShapeIt 2.0 -- NOT a mockup, this calls your actual
// ShapeSetting/ShapeMatrix/ShapeController/ShapeCollector classes unchanged.
// Only the Peaks and Energies panel is wired up; Mode/Options/Integration Bin
// are deliberately left for the next pass, once this path is confirmed working.
//
// Place this file in ShapeIt1.0/Source/ (next to ShapeController.h's siblings)
// alongside webshapeit.html, then run:  root WebShapeIt.cxx
//
// Known limitation, on purpose for now: bin size (sett->exi_size) is hardcoded
// below rather than wired to a UI panel yet, since Integration Bin isn't built
// yet. Mode is always "Integration" (mode=1); Autofit isn't wired yet either.

#include <ROOT/RWebWindow.hxx>
#include "TCanvas.h"
#include "TWebCanvas.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TLine.h"
#include "TBox.h"
#include "TH1.h"
#include "TVirtualPad.h"
#include "Buttons.h"
#include "TTimer.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TList.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>

#include "ShapeSetting.C"
#include "ShapeMatrix.C"
#include "ShapeGSF.C"
#include "ShapeCollector.C"
#include "ShapeAlpha.C"
#include "ShapeController.C"

std::shared_ptr<ROOT::RWebWindow> window;
TCanvas *canvas = nullptr;

ShapeSetting   *sett    = nullptr;
ShapeMatrix    *matrix  = nullptr;
ShapeCollector *gSFColl = nullptr;
std::string currentMatrixPath;
std::string gStartDir;

// mirrors ShapeFrame's displayMode: 1=matrix, 4=summed diagonal projection,
// 5=per-bin projection. Markers (Level 1/2 lines + background boxes) are only
// drawn/draggable in modes 4 and 5, same as the native GUI.
int gDisplayMode = 0;
int gCurrentBin = 1;
TLine *gMarkerLine[4] = { nullptr, nullptr, nullptr, nullptr };
TBox  *gBgBox[4]      = { nullptr, nullptr, nullptr, nullptr };
TH1   *gCurrentHist   = nullptr; // whatever's currently drawn in modes 4/5, used to size markers

// splits "685|900|1497|1732|3500|6700" style messages on '|'
std::vector<double> ParsePipeDoubles(const std::string &s)
{
    std::vector<double> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, '|'))
        out.push_back(std::stod(tok));
    return out;
}

// matches the native GUI's own "mname" convention: just the filename, used
// purely as a histogram title, not for anything functional
std::string BaseName(const std::string &path)
{
    auto pos = path.find_last_of("/\\");
    return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

std::string DirName(const std::string &path)
{
    auto pos = path.find_last_of("/\\");
    return (pos == std::string::npos) ? std::string(".") : path.substr(0, pos);
}

// canvas->Update() on a TWebCanvas blocks the entire ROOT process waiting for
// the browser to confirm it painted the frame -- calling it from inside
// ProcessData() creates a wait for an acknowledgment that can only arrive
// through the same message-processing path currently blocked, which was the
// actual cause of the "everything freezes for a minute" symptom (confirmed via
// `sample`, not guessed: TWebCanvas::WaitWhenCanvasPainted was ~100% of every
// sample taken during the freeze).
//
// TWebCanvas::ForceUpdate() is the documented non-blocking equivalent --
// "force sending data to client - do not wait for reply" -- so use this
// everywhere instead of Modified()+Update().
void PushCanvasUpdate()
{
    auto web_imp = dynamic_cast<TWebCanvas *>(canvas->GetCanvasImp());
    if (web_imp)
        web_imp->ForceUpdate();
}

// ShapeSetting::ReadSettings() stores dataFileName/osloFileName as literal raw
// lines from the settings file with no path resolution at all -- if they're
// relative (e.g. "Raw/Kr88.root"), they're relative to wherever ROOT happened
// to be launched from, not to the settings file's own directory. Resolving
// against the settings file's directory is the robust interpretation
// regardless of launch directory.
std::string ResolveRelativeTo(const std::string &baseDir, const std::string &maybeRelative)
{
    if (maybeRelative.empty() || maybeRelative[0] == '/')
        return maybeRelative; // already absolute
    return baseDir + "/" + maybeRelative;
}

void DumpSettings()
{
    std::cout << "--- sett dump after loading settings file ---\n"
              << "  mode: " << sett->mode << " (1=Integration, 2=Autofit)\n"
              << "  dataFileName: " << sett->dataFileName << "\n"
              << "  osloFileName: " << sett->osloFileName << "\n"
              << "  matrixName: " << sett->matrixName << "\n"
              << "  levEne: " << sett->levEne[0] << " " << sett->levEne[1] << " "
                               << sett->levEne[2] << " " << sett->levEne[3] << "\n"
              << "  exiEne: " << sett->exiEne[0] << " " << sett->exiEne[1] << "\n"
              << "  exi_size: " << sett->exi_size[0] << " " << sett->exi_size[1] << "\n"
              << "  doOslo: " << sett->doOslo << "  doAutoScale: " << sett->doAutoScale
                               << "  doInterpol: " << sett->doInterpol << "\n"
              << "  doSlidingWindow: " << sett->doSlidingWindow
                               << "  doBinVariation: " << sett->doBinVariation
                               << "  doBackground: " << sett->doBackground << "\n"
              << "----------------------------------------------\n";
}

// Recomputes nOfBins (low) from whatever exiEne/exi_size[0] currently are, and
// separately computes the high-bin-size equivalent using the explicit-size
// overload -- mirrors ShapeFrame::DoNumberEntry(), which recalculates both
// live as the user types into Excitation/bin-size fields, well before ever
// running ShapeIt. Only the low value is actually persisted in sett->nOfBins
// (ShapeSetting only has one such field); the high one is purely a live
// display/edit convenience, same as the native GUI's nOfBins[1] widget.
void SendBinSyncValues(unsigned connid, int nLo, int nHi)
{
    std::string msg = "NBINS:" + std::to_string(sett->exi_size[0]) + "|"
                                + std::to_string(sett->exi_size[1]) + "|"
                                + std::to_string(nLo) + "|"
                                + std::to_string(nHi);
    window->Send(connid, msg);
}

// Recomputes both bin counts fresh from the current exi_size/exiEne -- use this
// when bin SIZE (or the excitation range) changed, so the displayed counts
// reflect the new size. Do NOT use this after the user edits a bin COUNT
// directly (see NBINSLO/NBINSHI below) -- recomputing the count back from a
// rounded size can disagree with what was actually typed, which is exactly
// what caused the "always increases by ~2, never decreases" bug: editing
// nOfBins derived a size via BinToSize(), then this function immediately
// recomputed nOfBins again via SizeToBin() on that rounded size, silently
// overwriting the typed value with a different one.
void SendNBins(unsigned connid)
{
    sett->nOfBins = sett->SizeToBin();
    int nHigh = sett->doBinVariation ? sett->SizeToBin(sett->exi_size[1]) : sett->nOfBins;
    SendBinSyncValues(connid, sett->nOfBins, nHigh);
}

// Sends the list of matrix names found in the currently open file, and draws +
// selects the given 1-based index -- mirrors ShapeFrame::MatrixSelector() +
// ShapeFrame::MatrixSelect() combined into one step.
void SendMatrixListAndSelect(unsigned connid, const std::string &matrixPath, int selectIndex)
{
    auto names = matrix->GetMatrixName();

    std::string msg = "MATRIXLIST:" + std::to_string(selectIndex) + "\n";
    for (auto &n : names) msg += n + "\n";
    window->Send(connid, msg);

    matrix->SetMatrix(selectIndex);
    matrix->GetInputMatrix(BaseName(matrixPath))->Draw("colz");
    matrix->Diag();
    canvas->cd();
    PushCanvasUpdate();
}

// Sends the current sett state back to the browser so form fields can be kept
// in sync after loading a settings file (previously loading a settings file
// updated sett server-side but never told the browser, so the Peaks/Options
// panels silently went stale).
void SendSettingsSync(unsigned connid)
{
    std::string msg = "SETTINGS_SYNC:";
    msg += std::to_string(sett->levEne[0]) + "|" + std::to_string(sett->levEne[1]) + "|";
    msg += std::to_string(sett->levEne[2]) + "|" + std::to_string(sett->levEne[3]) + "|";
    msg += std::to_string(sett->exiEne[0]) + "|" + std::to_string(sett->exiEne[1]) + "|";
    msg += std::to_string(sett->levEne_2[0]) + "|" + std::to_string(sett->levEne_2[1]) + "|";
    msg += std::to_string(sett->levEne_2[2]) + "|" + std::to_string(sett->levEne_2[3]) + "|";
    msg += std::to_string(sett->doInterpol ? 1 : 0) + "|";
    msg += std::to_string(sett->doOslo ? 1 : 0) + "|";
    msg += std::to_string(sett->doSlidingWindow ? 1 : 0) + "|";
    msg += std::to_string(sett->doBackground ? 1 : 0) + "|";
    msg += std::to_string(sett->minCounts) + "|";
    msg += std::to_string(sett->gSF_norm) + "|";
    msg += std::to_string(sett->doAutoScale ? 1 : 0) + "|";
    msg += std::to_string(sett->eff_corr) + "|";
    msg += std::to_string(sett->mode) + "|";
    msg += std::to_string(sett->bgEne[0][0]) + "|" + std::to_string(sett->bgEne[0][1]) + "|";
    msg += std::to_string(sett->bgEne[0][2]) + "|" + std::to_string(sett->bgEne[0][3]) + "|";
    msg += std::to_string(sett->bgEne[1][0]) + "|" + std::to_string(sett->bgEne[1][1]) + "|";
    msg += std::to_string(sett->bgEne[1][2]) + "|" + std::to_string(sett->bgEne[1][3]);
    window->Send(connid, msg);
}

// Draws the Level 1/2 marker lines and background region boxes on whatever
// projection is currently displayed -- mirrors ShapeFrame::DrawMarker().
// Only meaningful in display modes 4 (summed projection) and 5 (per-bin
// projection), same as the native GUI.
//
// usePadRange: on the *initial* draw of a histogram, gPad's axis range isn't
// reliably established yet (see note below), so we size markers from the
// histogram's own data range instead. But once the user has actually
// interacted with the plot (zoomed/panned), gPad's range IS current and
// should be used instead, so markers resize/reposition to match -- and
// markers entirely outside the current X range are skipped rather than drawn
// off in space.
void DrawMarkers(bool usePadRange = false)
{
    for (int i = 0; i < 4; i++) {
        if (gMarkerLine[i]) canvas->GetListOfPrimitives()->Remove(gMarkerLine[i]);
        if (gBgBox[i]) canvas->GetListOfPrimitives()->Remove(gBgBox[i]);
    }

    if (gDisplayMode != 4 && gDisplayMode != 5) {
        PushCanvasUpdate();
        return;
    }

    double y1, y2;
    double xmin = -1e18, xmax = 1e18; // effectively "no clipping" unless usePadRange

    if (usePadRange && gPad) {
        y1 = gPad->GetUymin();
        y2 = gPad->GetUymax();
        if (gPad->GetLogy()) {
            y1 = TMath::Power(10, y1);
            y2 = TMath::Power(10, y2);
        }
        xmin = gPad->GetUxmin();
        xmax = gPad->GetUxmax();
    } else {
        // Previously always used gPad->GetUymin()/GetUymax() here too -- but
        // since we no longer call the blocking canvas->Update() (that was the
        // cause of the freeze bug fixed earlier), the pad's own axis-range
        // bookkeeping isn't reliably refreshed immediately after switching to
        // a new histogram. The histogram's own data range is independent of
        // that pad-painting timing and reflects the real data on first draw.
        y1 = 0;
        y2 = gCurrentHist ? gCurrentHist->GetMaximum() * 1.05 : 100;
    }

    for (int i = 0; i < 4; i++) {
        gMarkerLine[i] = new TLine(sett->levEne[i], y1, sett->levEne[i], y2);
        gMarkerLine[i]->SetLineColor(kRed);
        gMarkerLine[i]->SetLineWidth(2);
        if (sett->levEne[i] >= xmin && sett->levEne[i] <= xmax)
            gMarkerLine[i]->Draw();
    }

    gBgBox[0] = new TBox(sett->bgEne[0][0], y1, sett->bgEne[0][1], y2);
    gBgBox[1] = new TBox(sett->bgEne[0][2], y1, sett->bgEne[0][3], y2);
    gBgBox[2] = new TBox(sett->bgEne[1][0], y1, sett->bgEne[1][1], y2);
    gBgBox[3] = new TBox(sett->bgEne[1][2], y1, sett->bgEne[1][3], y2);

    for (int i = 0; i < 4; i++) {
        // Solid fill instead of the alpha-blended kBlue-9/-10 the native GUI
        // uses -- those are quite pale colors even at full opacity, and
        // combined with 45% alpha rendered as barely-visible grey in the
        // browser. kAzure/kOrange are more saturated and should stay visibly
        // colored; SetFillStyle(3013) gives a hatched "see-through" look
        // without depending on alpha compositing support.
        gBgBox[i]->SetFillColor(i < 2 ? kAzure + 1 : kOrange + 1);
        gBgBox[i]->SetFillStyle(3013);
        gBgBox[i]->SetLineColor(i < 2 ? kAzure + 1 : kOrange + 1);
        bool overlaps = gBgBox[i]->GetX2() >= xmin && gBgBox[i]->GetX1() <= xmax;
        if (sett->doBackground && overlaps)
            gBgBox[i]->Draw();
    }

    PushCanvasUpdate();
}

// Tracks the last-seen pad axis range, mirroring ShapeFrame's histY1/histY2
// members -- used to detect zoom/pan so markers can be redrawn to match.
double gLastUxmin = 0, gLastUxmax = 0, gLastUymin = 0, gLastUymax = 0;
bool gHaveLastRange = false;

// Checks whether the pad's visible axis range has changed since last checked,
// and redraws markers to match if so. Used both by HandleCanvasEvent (cheap,
// but confirmed NOT to fire for zoom/pan on a web canvas -- only for actual
// clicks/drags on objects) and by a periodic poll below (the mechanism that
// actually catches zoom/pan, since it doesn't depend on any event firing at
// all -- it just directly reads the pad's current state on a timer).
void CheckRangeChanged()
{
    if (!matrix || (gDisplayMode != 4 && gDisplayMode != 5) || !gPad)
        return;

    double uxmin = gPad->GetUxmin(), uxmax = gPad->GetUxmax();
    double uymin = gPad->GetUymin(), uymax = gPad->GetUymax();

    if (!gHaveLastRange || gLastUxmin != uxmin || gLastUxmax != uxmax
                        || gLastUymin != uymin || gLastUymax != uymax) {
        gLastUxmin = uxmin; gLastUxmax = uxmax;
        gLastUymin = uymin; gLastUymax = uymax;
        gHaveLastRange = true;
        DrawMarkers(true);
    }
}

// Called on every mouse event on the canvas -- mirrors ShapeFrame::HandleMyCanvas().
// This is a free function (not a class method), connected below via
// canvas->Connect() using ROOT's classic "global function slot" syntax, since
// this file has no custom dictionary-registered class to receive the signal.
//
// The actual drag interaction itself needs no code here at all -- JSROOT's
// built-in object editing already round-trips a dragged TLine/TBox's new
// coordinates back to the real C++ object automatically (confirmed early in
// this project: dragging a TLine in the browser shows up as real
// TLine::SetX1()/SetY1() etc. calls on the server, with zero custom code).
// This handler only needs to run *after* a drag finishes, to read the now-
// updated marker positions back into sett and keep the Peaks panel in sync.
// The CheckRangeChanged() call here is a cheap fallback for whatever event
// codes DO fire on this signal -- confirmed via testing that zoom/pan itself
// does NOT reach this handler at all, which is why the periodic poll below
// exists as the actual mechanism.
void HandleCanvasEvent(Int_t event, Int_t /*x*/, Int_t /*y*/, TObject * /*obj*/)
{
    if (!matrix || (gDisplayMode != 4 && gDisplayMode != 5))
        return;

    CheckRangeChanged();

    if (event == kButton1Up || event == kButton2Up || event == kButton3Up) {
        for (int i = 0; i < 4; i++)
            sett->levEne[i] = gMarkerLine[i]->GetX1();

        sett->bgEne[0][0] = gBgBox[0]->GetX1(); sett->bgEne[0][1] = gBgBox[0]->GetX2();
        sett->bgEne[0][2] = gBgBox[1]->GetX1(); sett->bgEne[0][3] = gBgBox[1]->GetX2();
        sett->bgEne[1][0] = gBgBox[2]->GetX1(); sett->bgEne[1][1] = gBgBox[2]->GetX2();
        sett->bgEne[1][2] = gBgBox[3]->GetX1(); sett->bgEne[1][3] = gBgBox[3]->GetX2();

        DrawMarkers();
        SendSettingsSync(0); // connid 0 broadcasts to all connections
    }
}

void RunShapeIt(unsigned connid)
{
    if (!matrix) {
        window->Send(connid, "No matrix loaded yet -- open one first.");
        return;
    }

    std::cout << "About to run with these settings:\n";
    DumpSettings();

    delete gSFColl;
    gSFColl = ShapeController::RunAnalysis(sett, matrix);

    // Testing whether the crash is about the DATA, about error bars specifically
    // (my earlier bare-TGraph test had none -- a real gap in that test), or
    // about TMultiGraph's own wrapping/serialization. This builds fresh copies
    // of each sub-graph's actual type (TGraphErrors or TGraphAsymmErrors) with
    // real error values, drawn directly -- no TMultiGraph involved at all.
    //
    // NEW: also classifies which of the sub-graphs is the literature/Oslo
    // comparison graph (vs. the actual gSF result graph), so it can be drawn
    // as a filled error band ("L3"/"AL3") instead of as points ("P"/"AP").
    // This is purely a change to how this diagnostic scaffolding draws --
    // still no TMultiGraph involved, so it doesn't touch the crash
    // investigation at all.
    TMultiGraph *diagGraph = gSFColl->getMultGraph();
    TList *graphList = diagGraph->GetListOfGraphs();
    std::vector<double> allX, allY;

    struct FreshGraph {
        TGraph *graph;
        bool isLiterature;
    };
    std::vector<FreshGraph> freshGraphs;
    
    // Get the literature graph pointer to reliably identify it
    // Literature data is loaded from osloFileName and only present if doOslo or doMC is true
    TGraph *litGraph = nullptr;
    if ((sett->doOslo || sett->doMC) && !sett->osloFileName.empty()) {
        litGraph = gSFColl->getLitGraph();
    }

    if (graphList) {
        TIter next(graphList);
        TObject *obj;
        int graphIdx = 0;
        while ((obj = next())) {
            TGraph *g = dynamic_cast<TGraph *>(obj);
            if (!g) { graphIdx++; continue; }
            std::cout << "Graph #" << graphIdx << " (" << g->GetName() << ", "
                      << obj->ClassName() << "), " << g->GetN() << " points:\n";

            // Classify by comparing pointer to the known literature graph
            // This is reliable because it checks if this graph IS the literature data
            // loaded from osloFileName, rather than guessing based on name
            bool isLit = (litGraph != nullptr && g == litGraph);

            if (auto *ge = dynamic_cast<TGraphAsymmErrors *>(g)) {
                std::vector<double> x, y, exl, exh, eyl, eyh;
                for (int i = 0; i < ge->GetN(); i++) {
                    x.push_back(ge->GetX()[i]); y.push_back(ge->GetY()[i]);
                    exl.push_back(ge->GetEXlow()[i]); exh.push_back(ge->GetEXhigh()[i]);
                    eyl.push_back(ge->GetEYlow()[i]); eyh.push_back(ge->GetEYhigh()[i]);
                    std::cout << "  [" << i << "] x=" << x.back() << " y=" << y.back()
                              << " eyl=" << eyl.back() << " eyh=" << eyh.back() << "\n";
                    allX.push_back(x.back()); allY.push_back(y.back());
                }
                auto *fresh = new TGraphAsymmErrors((int)x.size(), x.data(), y.data(),
                                                     exl.data(), exh.data(), eyl.data(), eyh.data());
                fresh->Sort();
                freshGraphs.push_back({fresh, isLit});
            }
            else if (auto *ge2 = dynamic_cast<TGraphErrors *>(g)) {
                std::vector<double> x, y, ex, ey;
                for (int i = 0; i < ge2->GetN(); i++) {
                    x.push_back(ge2->GetX()[i]); y.push_back(ge2->GetY()[i]);
                    ex.push_back(ge2->GetEX()[i]); ey.push_back(ge2->GetEY()[i]);
                    std::cout << "  [" << i << "] x=" << x.back() << " y=" << y.back()
                              << " ex=" << ex.back() << " ey=" << ey.back() << "\n";
                    allX.push_back(x.back()); allY.push_back(y.back());
                }
                auto *fresh = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
                fresh->Sort();
                freshGraphs.push_back({fresh, isLit});
            }
            else {
                for (int i = 0; i < g->GetN(); i++) {
                    std::cout << "  [" << i << "] x=" << g->GetX()[i] << " y=" << g->GetY()[i] << "\n";
                    allX.push_back(g->GetX()[i]); allY.push_back(g->GetY()[i]);
                }
            }
            graphIdx++;
        }
    }
    std::cout << "Total points collected: " << allX.size()
              << ", fresh error-bar graphs built: " << freshGraphs.size() << std::endl;

    canvas->cd();

    canvas->Clear();
    for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gBgBox[i] = nullptr; }

    bool firstDrawn = false;
    int colorIdx = 0;
    int colors[] = { kBlue, kRed, kGreen + 2, kMagenta, kOrange + 7, kCyan + 2 };
    for (auto &fg : freshGraphs) {
        TGraph *g = fg.graph;

        if (fg.isLiterature) {
            // Filled error band instead of points. Solid pale fill + hatch
            // style, not alpha -- alpha-blended fills render as near-invisible
            // grey in this browser/JSROOT combo, same issue already worked
            // around for the background boxes in DrawMarkers().
            g->SetFillColor(kBlue - 10);
            g->SetFillStyle(3013);
            g->SetLineColor(kBlue);
            g->SetLineWidth(2);
            g->SetTitle("Literature gSF (error band)");
            g->Draw(firstDrawn ? "L3 SAME" : "AL3");
        } else {
            g->SetMarkerStyle(20);
            g->SetMarkerColor(colors[colorIdx % 6]);
            g->SetLineColor(colors[colorIdx % 6]);
            g->SetTitle("gSF values (fresh error-bar graphs, no TMultiGraph)");
            g->Draw(firstDrawn ? "P SAME" : "AP");
            colorIdx++;
        }
        firstDrawn = true;
    }
    // Fallback: if no error-bar graphs were found (e.g. nothing matched
    // TGraphErrors/TGraphAsymmErrors), fall back to the plain-TGraph test
    // from before, so this still shows something.
    if (!firstDrawn && !allX.empty()) {
        TGraph *simpleGraph = new TGraph((int)allX.size(), allX.data(), allY.data());
        simpleGraph->SetMarkerStyle(20);
        simpleGraph->SetMarkerColor(kBlue);
        simpleGraph->SetTitle("gSF values (bare TGraph test plot, no error-bar graphs found)");
        simpleGraph->Draw("AP");
    }

    // canvas->Clear() just deleted any marker TLine/TBox objects left over
    // from whatever projection view was showing before -- but gDisplayMode
    // was still 4/5, so the 100ms poll timer would soon call DrawMarkers()
    // again, whose first action is Remove()-ing these now-dangling pointers
    // (already nulled above, right after Clear()). That's a use-after-free
    // that plausibly corrupts the same primitive list CreatePadSnapshot walks
    // moments later -- which lines up exactly with where this crash happens.
    // This is display mode 0: "results view, no markers apply".
    gDisplayMode = 0;
    gHaveLastRange = false;

    PushCanvasUpdate();

    window->Send(connid, "Done.");
    SendNBins(connid);
}

// Lists a directory's contents and sends it back as a simple newline-delimited
// message: "DIRLIST:<resolved path>\nD:name1\nD:name2\nF:file1.root\n..."
// D/F prefixes distinguish directories from files. The browser only ever
// displays what the server tells it -- it never touches the real filesystem,
// which sidesteps the browser's file-picker sandbox limitation entirely.
void SendDirListing(unsigned connid, std::string path)
{
    if (path.empty())
        path = gSystem->WorkingDirectory();

    if (gSystem->AccessPathName(path.c_str())) {
        window->Send(connid, "DIRERROR:Path not found: " + path);
        return;
    }

    void *dir = gSystem->OpenDirectory(path.c_str());
    if (!dir) {
        window->Send(connid, "DIRERROR:Could not open directory: " + path);
        return;
    }

    std::vector<std::string> dirs, files;
    const char *entry;
    while ((entry = gSystem->GetDirEntry(dir))) {
        std::string name(entry);
        if (name == "." || name == "..") continue;

        std::string full = path + "/" + name;
        struct stat st;
        if (stat(full.c_str(), &st) != 0) continue;

        if (S_ISDIR(st.st_mode))
            dirs.push_back(name);
        else
            files.push_back(name);
    }
    gSystem->FreeDirectory(dir);

    std::sort(dirs.begin(), dirs.end());
    std::sort(files.begin(), files.end());

    std::string msg = "DIRLIST:" + path + "\n";
    for (auto &d : dirs)  msg += "D:" + d + "\n";
    for (auto &f : files) msg += "F:" + f + "\n";

    window->Send(connid, msg);
}

void ProcessData(unsigned connid, const std::string &arg)
{
    std::cout << "Got message from browser: " << arg << std::endl;

    if (arg.compare(0, 8, "channel:") == 0) {
        int chid = std::stoi(arg.substr(8));
        auto web_imp = dynamic_cast<TWebCanvas *>(canvas->GetCanvasImp());
        if (web_imp) {
            web_imp->ShowWebWindow({ window, connid, chid });
            web_imp->ForceUpdate();
        }
        window->Send(connid, "STARTDIR:" + gStartDir);
    }
    else if (arg.compare(0, 5, "OPEN:") == 0) {
        std::string path = arg.substr(5);

        if (gSystem->AccessPathName(path.c_str())) {
            // AccessPathName returns non-zero (true) when the path does NOT exist --
            // this used to go straight into ShapeMatrix's constructor with no check
            // at all, crashing hard on any typo'd or invalid path.
            window->Send(connid, "File not found: " + path);
            return;
        }

        sett->SetFileName(path);
        currentMatrixPath = path;
        delete matrix;
        matrix = nullptr;
        matrix = new ShapeMatrix(sett);
        SendMatrixListAndSelect(connid, currentMatrixPath, 1);
        SendNBins(connid);
        window->Send(connid, "Matrix opened: " + path);
    }
    else if (arg.compare(0, 13, "SELECTMATRIX:") == 0) {
        if (!matrix) {
            window->Send(connid, "No matrix file open yet.");
            return;
        }
        int idx = std::stoi(arg.substr(13));
        SendMatrixListAndSelect(connid, currentMatrixPath, idx);
        SendNBins(connid);
        window->Send(connid, "Matrix selected.");
    }
    else if (arg.compare(0, 5, "OSLO:") == 0) {
        std::string path = arg.substr(5);

        if (gSystem->AccessPathName(path.c_str())) {
            window->Send(connid, "Literature file not found: " + path);
            return;
        }

        sett->osloFileName = path;
        sett->doOslo = true;
        window->Send(connid, "Literature file set: " + path);
    }
    else if (arg.compare(0, 8, "LISTDIR:") == 0) {
        SendDirListing(connid, arg.substr(8));
    }
    else if (arg.compare(0, 14, "LOAD_SETTINGS:") == 0) {
        std::string path = arg.substr(14);

        if (gSystem->AccessPathName(path.c_str())) {
            window->Send(connid, "Settings file not found: " + path);
            return;
        }

        sett->settFileName = path;
        sett->ReadSettings();

        // dataFileName/osloFileName as read are literal raw lines from the file,
        // unresolved -- if relative, resolve against the settings file's own
        // directory rather than wherever ROOT happened to be launched from.
        std::string settDir = DirName(path);
        sett->dataFileName = ResolveRelativeTo(settDir, sett->dataFileName);
        sett->osloFileName = ResolveRelativeTo(settDir, sett->osloFileName);

        DumpSettings();

        // mirrors ShapeFrame::OpenSettingFile -- the settings file references
        // both a data file and a specific matrix name within it
        if (sett->dataFileName.empty()) {
            window->Send(connid, "Settings loaded, but no data file listed in it.");
        }
        else if (gSystem->AccessPathName(sett->dataFileName.c_str())) {
            window->Send(connid, "Settings loaded, but matrix file not found: " + sett->dataFileName);
        }
        else {
            currentMatrixPath = sett->dataFileName;
            delete matrix;
            matrix = nullptr;
            matrix = new ShapeMatrix(sett);

            auto names = matrix->GetMatrixName();
            int idx = 0;
            for (size_t i = 0; i < names.size(); i++)
                if (names[i] == sett->matrixName) idx = (int)i + 1;

            if (idx > 0) {
                SendMatrixListAndSelect(connid, currentMatrixPath, idx);
                SendNBins(connid);
            }
            else
                window->Send(connid, "Warning: matrix '" + sett->matrixName + "' from settings not found in " + sett->dataFileName);
        }

        SendSettingsSync(connid);
        window->Send(connid, "Settings loaded: " + path);
    }
    else if (arg.compare(0, 14, "SAVE_SETTINGS:") == 0) {
        // Uses sett as it currently stands -- i.e. whatever the last "ShapeIt!"
        // run (or Options change) set it to. Click ShapeIt! at least once before
        // saving so the Peaks panel's current values are actually captured.
        std::string path = arg.substr(14);
        sett->settFileName = path;
        sett->SaveSettings();
        window->Send(connid, "Settings saved: " + path);
    }
    else if (arg.compare(0, 8, "OPTIONS:") == 0) {
        // order: doInterpol|doOslo|doSlidingWindow|doBackground
        // (doBinVariation moved to the Integration Bin panel/BINSIZE message)
        auto v = ParsePipeDoubles(arg.substr(8));
        if (v.size() != 4) {
            window->Send(connid, "Malformed OPTIONS message.");
            return;
        }
        sett->doInterpol      = v[0] != 0.0;
        sett->doOslo          = v[1] != 0.0;
        sett->doSlidingWindow = v[2] != 0.0;
        sett->doBackground    = v[3] != 0.0;
        window->Send(connid, "Options updated.");
    }
    else if (arg.compare(0, 8, "BINSIZE:") == 0) {
        // order: lo|hi|isVariation
        auto v = ParsePipeDoubles(arg.substr(8));
        if (v.size() != 3) {
            window->Send(connid, "Malformed BINSIZE message.");
            return;
        }
        sett->exi_size[0] = v[0];
        sett->doBinVariation = v[2] != 0.0;
        double hi = v[1];
        if (sett->doBinVariation && hi <= sett->exi_size[0] + 50)
            hi = sett->exi_size[0] + 50; // mirrors DoNumberEntry's id==11 clamp
        sett->exi_size[1] = sett->doBinVariation ? hi : sett->exi_size[0];
        SendNBins(connid);
        window->Send(connid, "Bin size updated.");
    }
    else if (arg.compare(0, 8, "NBINSLO:") == 0) {
        // mirrors DoNumberEntry's id==8 case: editing "Nr. of bins" (low)
        // recomputes the corresponding bin size via BinToSize(). Trusts the
        // typed count as-is -- does NOT recompute it back from the derived
        // size, which is what caused the earlier drift bug.
        int nLo = std::stoi(arg.substr(8));
        if (nLo < 1) nLo = 1;
        sett->nOfBins = nLo;
        sett->exi_size[0] = sett->BinToSize();
        if (!sett->doBinVariation)
            sett->exi_size[1] = sett->exi_size[0];
        int nHigh = sett->doBinVariation ? sett->SizeToBin(sett->exi_size[1]) : sett->nOfBins;
        SendBinSyncValues(connid, sett->nOfBins, nHigh);
        window->Send(connid, "Bin size updated.");
    }
    else if (arg.compare(0, 8, "NBINSHI:") == 0) {
        // mirrors DoNumberEntry's id==12 case exactly, including its two
        // correction passes: hi bin size must exceed lo by at least 50 keV,
        // and the resulting bin count must be at least 3. sett->nOfBins (the
        // low count) is untouched here -- echoed back as-is, not recomputed.
        int nHigh = std::stoi(arg.substr(8));
        double size = sett->BinToSize(nHigh);
        if (size <= sett->exi_size[0] + 50) {
            size = sett->exi_size[0] + 50;
            nHigh = sett->SizeToBin(size);
        }
        if (nHigh < 3) {
            nHigh = 3;
            size = sett->BinToSize(nHigh);
        }
        sett->exi_size[1] = size;
        SendBinSyncValues(connid, sett->nOfBins, nHigh);
        window->Send(connid, "Bin size updated.");
    }
    else if (arg.compare(0, 10, "INTPARAMS:") == 0) {
        // order: minCounts|scaling|autoScale|effCorr
        auto v = ParsePipeDoubles(arg.substr(10));
        if (v.size() != 4) {
            window->Send(connid, "Malformed INTPARAMS message.");
            return;
        }
        sett->minCounts   = (int)v[0];
        sett->gSF_norm    = v[1];
        sett->doAutoScale = v[2] != 0.0;
        sett->eff_corr    = v[3];
        window->Send(connid, "Integration bin parameters updated.");
    }
    else if (arg.compare(0, 8, "VERBOSE:") == 0) {
        sett->verbose = std::stoi(arg.substr(8));
        window->Send(connid, "Verbose level set to " + std::to_string(sett->verbose) + ".");
    }
    else if (arg.compare(0, 5, "MODE:") == 0) {
        sett->mode = std::stoi(arg.substr(5)); // 1 = Integration, 2 = Autofit
        window->Send(connid, sett->mode == 2 ? "Mode: Autofit" : "Mode: Integration");
    }
    else if (arg.compare(0, 11, "BGENERGIES:") == 0) {
        // order: bgEne[0][0..3] | bgEne[1][0..3]  (8 values)
        auto v = ParsePipeDoubles(arg.substr(11));
        if (v.size() != 8) {
            window->Send(connid, "Malformed BGENERGIES message.");
            return;
        }
        for (int i = 0; i < 4; i++) sett->bgEne[0][i] = v[i];
        for (int i = 0; i < 4; i++) sett->bgEne[1][i] = v[4 + i];
        DrawMarkers(); // refreshes the draggable boxes to match, if a projection is shown
        window->Send(connid, "Background regions updated.");
    }
    else if (arg == "SHOWMATRIX") {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet -- open one first.");
            return;
        }
        gDisplayMode = 1;
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gBgBox[i] = nullptr; } // Clear() just deleted these -- see RunShapeIt for the full explanation
        matrix->GetInputMatrix(BaseName(currentMatrixPath))->Draw("colz");
        DrawMarkers(); // no-op in mode 1 (removes any leftover markers), matches native
    }
    else if (arg == "SHOWPROJ") {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet -- open one first.");
            return;
        }
        gDisplayMode = 4;
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gBgBox[i] = nullptr; }
        gCurrentHist = matrix->GetDiag(BaseName(currentMatrixPath));
        gCurrentHist->Draw("hist");
        gHaveLastRange = false;
        DrawMarkers();
    }
    else if (arg.compare(0, 12, "SHOWBINPROJ:") == 0) {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet -- open one first.");
            return;
        }
        gCurrentBin = std::stoi(arg.substr(12));
        gDisplayMode = 5;
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gBgBox[i] = nullptr; }
        gCurrentHist = matrix->GetDiagEx(gCurrentBin, BaseName(currentMatrixPath));
        gCurrentHist->Draw();
        gHaveLastRange = false;
        DrawMarkers();
    }
    else if (arg.compare(0, 11, "EXCITATION:") == 0) {
        auto v = ParsePipeDoubles(arg.substr(11));
        if (v.size() == 2) {
            sett->exiEne[0] = v[0];
            sett->exiEne[1] = v[1];
            SendNBins(connid);
        }
    }
    else if (arg.compare(0, 4, "RUN:") == 0) {
        // expected order: lvl1_lo|lvl1_hi|lvl2_lo|lvl2_hi|exc_lo|exc_hi|
        //                  is_doublet1|d1_lo|d1_hi|is_doublet2|d2_lo|d2_hi
        auto v = ParsePipeDoubles(arg.substr(4));
        if (v.size() != 12) {
            window->Send(connid, "Malformed RUN message.");
            return;
        }

        sett->levEne[0] = v[0]; sett->levEne[1] = v[1];
        sett->levEne[2] = v[2]; sett->levEne[3] = v[3];
        sett->exiEne[0] = v[4]; sett->exiEne[1] = v[5];

        bool isDoublet1 = v[6] != 0.0;
        sett->levEne_2[0] = isDoublet1 ? v[7] : 0.0;
        sett->levEne_2[1] = isDoublet1 ? v[8] : 0.0;

        bool isDoublet2 = v[9] != 0.0;
        sett->levEne_2[2] = isDoublet2 ? v[10] : 0.0;
        sett->levEne_2[3] = isDoublet2 ? v[11] : 0.0;

        RunShapeIt(connid);
    }
}

void WebShapeIt()
{
    gEnv->SetValue("WebGui.ConnCredits", "100");
    gStartDir = gSystem->WorkingDirectory();

    sett = new ShapeSetting();
    // Hardcoded placeholder until the Integration Bin panel is wired up.
    sett->exi_size[0] = 400;
    sett->exi_size[1] = 400;
    sett->exiEne[0] = 2500; // matches the native GUI's default excitation-range widget values
    sett->exiEne[1] = 7000;
    sett->mode = 1; // default to Integration mode; toggled via MODE: message
    sett->doBackground = false; // ShapeSetting defaults this to true, but background regions
                                 // are dataset-specific and unset here -- leaving it on with
                                 // zero-width regions is what crashed ShapeCollector::Norm().
                                 // Turn back on once the Options panel lets you set real values.

    // Generic placeholder background regions, matching the native GUI's own constructor
    // defaults exactly -- not physically correct for any specific dataset, but prevents
    // a crash if Background subtraction gets checked via the Options panel before real
    // per-dataset regions are configurable.
    double bg1[4] = {260, 360, 700, 800};
    double bg2[4] = {850, 950, 1350, 1450};
    sett->setBgEne1(bg1);
    sett->setBgEne2(bg2);

    canvas = TWebCanvas::CreateWebCanvas("webshapeit_canvas", "ShapeIt 2.0");

    // Connects to a plain global function (no custom dictionary-registered
    // class needed) -- this is what lets DrawMarkers()/HandleCanvasEvent() know
    // when a dragged marker/box has been released, to read its new position
    // back into sett.
    canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
                     "HandleCanvasEvent(Int_t,Int_t,Int_t,TObject*)");

    // Polls the pad's axis range directly every 100ms and redraws markers if
    // it changed -- this is the actual mechanism for zoom-adaptive markers.
    // Confirmed by testing that zoom/pan on a web canvas doesn't fire
    // ProcessedEvent at all (unlike dragging an object, which does), so
    // there's no event to hook for this -- polling sidesteps that entirely
    // by just directly checking the pad's current state instead of waiting
    // for something to tell us it changed. CheckRangeChanged() only reads
    // gPad's state and calls the non-blocking PushCanvasUpdate() -- it never
    // calls the blocking canvas->Update(), so this doesn't reintroduce the
    // freeze bug fixed earlier.
    static TTimer *zoomPollTimer = new TTimer();
    zoomPollTimer->Connect("Timeout()", 0, 0, "CheckRangeChanged()");
    zoomPollTimer->Start(100, kFALSE);

    window = ROOT::RWebWindow::Create();
    std::string fname = __FILE__;
    auto pos = fname.find("WebShapeIt.cxx");
    std::string dir = (pos != std::string::npos) ? fname.substr(0, pos) : std::string("./");
    window->SetDefaultPage("file:" + dir + "webshapeit.html");
    window->SetDataCallBack(ProcessData);
    window->SetGeometry(1200, 700);
    window->Show();

    std::cout << "\nShapeIt 2.0 prototype running.\n";
    std::cout << "(Note: button clicks may take up to ~1 minute to take effect --\n"
              << "this is a known, unresolved delay from earlier testing, not new.)\n";
}
