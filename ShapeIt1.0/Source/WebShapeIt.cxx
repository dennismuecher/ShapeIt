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
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TList.h"
#include <iostream>
#include <fstream>
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
TLine *gDoubletLine[4] = { nullptr, nullptr, nullptr, nullptr }; // Doublet markers
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

// Parse command-line arguments passed to ROOT
// Also checks environment variable SHAPEIT_SETTINGS as fallback
std::string GetCmdLineArg(const std::string &flag)
{
    TApplication *app = gApplication;
    if (!app) {
        // If no gApplication, check environment variable
        if (flag == "--settings") {
            const char* env = gSystem->Getenv("SHAPEIT_SETTINGS");
            return env ? std::string(env) : "";
        }
        return "";
    }
    
    Int_t argc = app->Argc();
    char **argv = app->Argv();
    
    for (Int_t i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        
        // Support both --settings=path and --settings path formats
        if (arg.find(flag + "=") == 0) {
            return arg.substr(flag.length() + 1);
        }
        else if (arg == flag && i + 1 < argc) {
            return std::string(argv[i + 1]);
        }
    }
    
    // Fallback: check environment variable
    if (flag == "--settings") {
        const char* env = gSystem->Getenv("SHAPEIT_SETTINGS");
        return env ? std::string(env) : "";
    }
    
    return "";
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
    
    // Set pad margins before drawing
    canvas->cd();
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.12);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.10);
    
    TH2* hist = matrix->GetInputMatrix(BaseName(matrixPath));
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("E_{#gamma} (keV)");
    hist->GetYaxis()->SetTitle("E_{x} (keV)");
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetTitleOffset(1.0);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->Draw("col");
    
    // Manually create and position the color palette
    double xmin = hist->GetXaxis()->GetXmin();
    double xmax = hist->GetXaxis()->GetXmax();
    double ymin = hist->GetYaxis()->GetXmin();
    double ymax = hist->GetYaxis()->GetXmax();
    TPaletteAxis *palette = new TPaletteAxis(xmax, ymin, xmax + (xmax-xmin)*0.05, ymax, hist);
    palette->SetX1NDC(0.86);
    palette->SetX2NDC(0.89);
    palette->SetY1NDC(0.10);
    palette->SetY2NDC(0.90);
    palette->Draw();
    
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
    msg += std::to_string(sett->doDoublet[0] ? 1 : 0) + "|";  // Add doublet checkbox states
    msg += std::to_string(sett->doDoublet[1] ? 1 : 0) + "|";
    msg += std::to_string(sett->fixDoubletWidth[0] ? 1 : 0) + "|";  // Add doublet width fix toggles
    msg += std::to_string(sett->fixDoubletWidth[1] ? 1 : 0) + "|";
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
    msg += std::to_string(sett->bgEne[1][2]) + "|" + std::to_string(sett->bgEne[1][3]) + "|";
    msg += std::to_string(sett->displaySingle ? 1 : 0) + "|";
    msg += std::to_string(sett->displayAvg ? 1 : 0) + "|";
    msg += std::to_string(sett->colour ? 1 : 0) + "|";
    msg += std::to_string(sett->doWidthCal ? 1 : 0) + "|";
    msg += std::to_string(sett->doBinVariation ? 1 : 0) + "|";
    msg += std::to_string(sett->exi_size[0]) + "|" + std::to_string(sett->exi_size[1]) + "|";
    msg += std::to_string(sett->verbose) + "|";
    msg += std::to_string(sett->widthCal[0][0]) + "|" + std::to_string(sett->widthCal[0][1]) + "|";
    msg += std::to_string(sett->widthCal[1][0]) + "|" + std::to_string(sett->widthCal[1][1]);
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
        if (gDoubletLine[i]) canvas->GetListOfPrimitives()->Remove(gDoubletLine[i]);
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

    // Draw main peak markers (red)
    for (int i = 0; i < 4; i++) {
        gMarkerLine[i] = new TLine(sett->levEne[i], y1, sett->levEne[i], y2);
        gMarkerLine[i]->SetLineColor(kRed);
        gMarkerLine[i]->SetLineWidth(2);
        if (sett->levEne[i] >= xmin && sett->levEne[i] <= xmax)
            gMarkerLine[i]->Draw();
    }

    // Draw doublet markers (orange) if enabled via checkbox
    // levEne_2[0-1] are for level 1 doublet, levEne_2[2-3] are for level 2 doublet
    // Use doDoublet flag, not zero-detection of energy values
    
    if (sett->doDoublet[0]) {
        for (int i = 0; i < 2; i++) {
            gDoubletLine[i] = new TLine(sett->levEne_2[i], y1, sett->levEne_2[i], y2);
            gDoubletLine[i]->SetLineColor(kOrange);
            gDoubletLine[i]->SetLineWidth(2);
            if (sett->levEne_2[i] >= xmin && sett->levEne_2[i] <= xmax)
                gDoubletLine[i]->Draw();
        }
    }
    
    if (sett->doDoublet[1]) {
        for (int i = 2; i < 4; i++) {
            gDoubletLine[i] = new TLine(sett->levEne_2[i], y1, sett->levEne_2[i], y2);
            gDoubletLine[i]->SetLineColor(kOrange);
            gDoubletLine[i]->SetLineWidth(2);
            if (sett->levEne_2[i] >= xmin && sett->levEne_2[i] <= xmax)
                gDoubletLine[i]->Draw();
        }
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

// Track last marker positions to detect when they're dragged
double gLastLevEne[4] = {0, 0, 0, 0};
double gLastLevEne_2[4] = {0, 0, 0, 0}; // Doublet positions
double gLastBgEne[2][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
bool gHaveLastMarkers = false;

// Track width calibration plot axis ranges (set once when plot is created)
double gWidthCalibXMin = 0, gWidthCalibXMax = 0, gWidthCalibYMin = 0, gWidthCalibYMax = 0;
bool gHaveWidthCalibRanges = false;

// Forward declarations
void RunShapeIt(unsigned connid);
void RunWidthCalibration(unsigned connid);

// Runs a temporary single-iteration autofit analysis purely to generate width
// calibration data, without affecting the user's actual settings for sliding
// window or bin variation. This allows width calibration to be viewed anytime
// (even before pressing ShapeIt) and always uses optimal settings (single bin
// size, no sliding window) to get the maximum number of data points.
void RunWidthCalibration(unsigned connid)
{
    std::cout << "=== Running width calibration analysis ===" << std::endl;
    
    if (!matrix) {
        window->Send(connid, "No matrix loaded yet.");
        return;
    }
    
    // Save ALL current settings that we'll temporarily override
    // This ensures the user's settings are completely untouched
    int savedMode = sett->mode;
    bool savedSlidingWindow = sett->doSlidingWindow;
    bool savedBinVariation = sett->doBinVariation;
    double savedExiSize1 = sett->exi_size[1];  // Save high bin size too
    
    std::cout << "Saving user settings - Mode: " << savedMode 
              << ", SlidingWindow: " << savedSlidingWindow 
              << ", BinVariation: " << savedBinVariation << std::endl;
    
    // Force settings optimal for width calibration
    sett->mode = 2;  // Autofit mode required for width data
    sett->doSlidingWindow = false;  // No sliding window - want single peaks
    sett->doBinVariation = false;   // Single bin size only - maximum data points
    
    std::cout << "Running single-iteration autofit for width calibration..." << std::endl;
    std::cout << "  (Temporarily using: Autofit mode, no sliding window, no bin variation)" << std::endl;
    
    // Run a temporary analysis just to populate the width data in the matrix
    // We discard the ShapeCollector results - only the width fits matter
    ShapeCollector *tempColl = ShapeController::RunAnalysis(sett, matrix);
    delete tempColl;
    
    // Restore ALL original settings - user should see no change
    sett->mode = savedMode;
    sett->doSlidingWindow = savedSlidingWindow;
    sett->doBinVariation = savedBinVariation;
    sett->exi_size[1] = savedExiSize1;  // Restore high bin size
    
    std::cout << "Width calibration complete. User settings restored." << std::endl;
    std::cout << "  Restored - Mode: " << sett->mode 
              << ", SlidingWindow: " << sett->doSlidingWindow 
              << ", BinVariation: " << sett->doBinVariation << std::endl;
    
    // Send the fitted parameters to the UI to populate the fields
    std::string msg = "WIDTH_CALIB_PARAMS:";
    msg += std::to_string(sett->widthCal[0][0]) + "|" + std::to_string(sett->widthCal[0][1]) + "|";
    msg += std::to_string(sett->widthCal[1][0]) + "|" + std::to_string(sett->widthCal[1][1]);
    window->Send(connid, msg);
}

// Clean up autofit display: remove intermediate cyan fits and add clean background lines
void CleanupAutofitDisplay()
{
    if (sett->mode != 2 || !gCurrentHist)
        return;
    
    TList *funcs = gCurrentHist->GetListOfFunctions();
    if (!funcs)
        return;
    
    // Remove all cyan intermediate fits (color 426)
    std::vector<TF1*> toRemove;
    TIter next1(funcs);
    TObject *obj1;
    while ((obj1 = next1())) {
        TF1 *f = dynamic_cast<TF1*>(obj1);
        if (f && f->GetLineColor() == 426) {
            toRemove.push_back(f);
        }
    }
    for (TF1 *f : toRemove) {
        funcs->Remove(f);
    }
    
    // Extract background parameters from final fits and draw background lines
    TIter next2(funcs);
    TObject *obj2;
    while ((obj2 = next2())) {
        TF1 *f = dynamic_cast<TF1*>(obj2);
        if (f && f->GetLineColor() == 6) { // Only process the final fit (magenta)
            std::string fname = f->GetName();
            int level = -1;
            if (fname.find("level1") != std::string::npos) level = 0;
            else if (fname.find("level2") != std::string::npos) level = 1;
            
            if (level >= 0 && level <= 1) {
                // Extract background parameters: par[1]*x + par[2]
                double slope = f->GetParameter(1);
                double intercept = f->GetParameter(2);
                
                // Create a simple linear function for just the background
                std::string bgName = "bg_line_level" + std::to_string(level+1) + "_bin" + std::to_string(gCurrentBin);
                TF1 *bgLine = new TF1(bgName.c_str(), "[0]*x + [1]", 
                                      sett->bgEne[level][0], sett->bgEne[level][3]);
                bgLine->SetParameter(0, slope);
                bgLine->SetParameter(1, intercept);
                bgLine->SetLineColor(kCyan-6);
                bgLine->SetLineWidth(2);
                bgLine->SetLineStyle(2); // Dashed line
                bgLine->SetNpx(500);
                bgLine->Draw("SAME");
            }
        }
    }
}

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

// Checks if markers have been dragged and updates settings if so
void CheckMarkersChanged()
{
    if (!matrix || (gDisplayMode != 4 && gDisplayMode != 5))
        return;
    
    if (!gMarkerLine[0] || !gBgBox[0])
        return;
    
    // Read current marker positions
    double levEne[4];
    double levEne_2[4];
    double bgEne[2][4];
    
    for (int i = 0; i < 4; i++)
        levEne[i] = gMarkerLine[i]->GetX1();
    
    // Read doublet positions ONLY if they were actually drawn (checkbox is ON)
    // If checkbox is OFF (lines not drawn), preserve existing stored values unchanged
    for (int i = 0; i < 4; i++) {
        if (gDoubletLine[i])
            levEne_2[i] = gDoubletLine[i]->GetX1();  // User dragged it, use new position
        else
            levEne_2[i] = sett->levEne_2[i];  // Checkbox off, preserve stored value (don't modify)
    }
    
    bgEne[0][0] = gBgBox[0]->GetX1(); bgEne[0][1] = gBgBox[0]->GetX2();
    bgEne[0][2] = gBgBox[1]->GetX1(); bgEne[0][3] = gBgBox[1]->GetX2();
    bgEne[1][0] = gBgBox[2]->GetX1(); bgEne[1][1] = gBgBox[2]->GetX2();
    bgEne[1][2] = gBgBox[3]->GetX1(); bgEne[1][3] = gBgBox[3]->GetX2();
    
    // Check if anything changed
    bool changed = false;
    if (gHaveLastMarkers) {
        for (int i = 0; i < 4; i++) {
            if (std::abs(levEne[i] - gLastLevEne[i]) > 0.01) {
                changed = true;
                break;
            }
            if (std::abs(levEne_2[i] - gLastLevEne_2[i]) > 0.01) {
                changed = true;
                break;
            }
        }
        if (!changed) {
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 4; j++) {
                    if (std::abs(bgEne[i][j] - gLastBgEne[i][j]) > 0.01) {
                        changed = true;
                        break;
                    }
                }
                if (changed) break;
            }
        }
    }
    
    // Update stored values
    for (int i = 0; i < 4; i++) {
        gLastLevEne[i] = levEne[i];
        gLastLevEne_2[i] = levEne_2[i];
    }
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 4; j++)
            gLastBgEne[i][j] = bgEne[i][j];
    gHaveLastMarkers = true;
    
    // If changed, update settings and UI
    if (changed) {
        std::cout << "Marker position changed, updating settings..." << std::endl;
        
        for (int i = 0; i < 4; i++) {
            sett->levEne[i] = levEne[i];
            sett->levEne_2[i] = levEne_2[i];  // This now preserves the stored value when doublet is off
        }
        
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 4; j++)
                sett->bgEne[i][j] = bgEne[i][j];
        
        // Send only marker positions to UI, not full settings sync
        // This avoids overwriting checkbox states that may have changed in the UI
        // but not yet been sent back to C++
        std::string msg = "MARKER_UPDATE:";
        msg += std::to_string(sett->levEne[0]) + "|" + std::to_string(sett->levEne[1]) + "|";
        msg += std::to_string(sett->levEne[2]) + "|" + std::to_string(sett->levEne[3]) + "|";
        msg += std::to_string(sett->levEne_2[0]) + "|" + std::to_string(sett->levEne_2[1]) + "|";
        msg += std::to_string(sett->levEne_2[2]) + "|" + std::to_string(sett->levEne_2[3]) + "|";
        msg += std::to_string(sett->bgEne[0][0]) + "|" + std::to_string(sett->bgEne[0][1]) + "|";
        msg += std::to_string(sett->bgEne[0][2]) + "|" + std::to_string(sett->bgEne[0][3]) + "|";
        msg += std::to_string(sett->bgEne[1][0]) + "|" + std::to_string(sett->bgEne[1][1]) + "|";
        msg += std::to_string(sett->bgEne[1][2]) + "|" + std::to_string(sett->bgEne[1][3]);
        window->Send(0, msg);
        
        // If in Autofit mode (mode 2), re-fit and redraw the current projection
        // This updates the Gaussian fits on the histogram without recalculating all gSF
        if (sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
            std::cout << "Autofit mode: updating fits for bin " << gCurrentBin << "..." << std::endl;
            
            // Save current axis ranges before redrawing
            double xmin = gPad->GetUxmin();
            double xmax = gPad->GetUxmax();
            double ymin = gPad->GetUymin();
            double ymax = gPad->GetUymax();
            bool isLogy = gPad->GetLogy();
            
            canvas->cd();
            gCurrentHist = matrix->GetDiagEx(gCurrentBin, BaseName(currentMatrixPath));
            
            // Restore axis ranges
            gCurrentHist->GetXaxis()->SetRangeUser(xmin, xmax);
            if (isLogy)
                gCurrentHist->GetYaxis()->SetRangeUser(TMath::Power(10, ymin), TMath::Power(10, ymax));
            else
                gCurrentHist->GetYaxis()->SetRangeUser(ymin, ymax);
            
            gCurrentHist->Draw();
            CleanupAutofitDisplay();
            DrawMarkers(true);
            PushCanvasUpdate();
        }
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

        // Update doublet positions if they exist (only if doublet lines were drawn)
        for (int i = 0; i < 4; i++) {
            if (gDoubletLine[i])
                sett->levEne_2[i] = gDoubletLine[i]->GetX1();
            // else: preserve existing sett->levEne_2[i] value unchanged
        }

        sett->bgEne[0][0] = gBgBox[0]->GetX1(); sett->bgEne[0][1] = gBgBox[0]->GetX2();
        sett->bgEne[0][2] = gBgBox[1]->GetX1(); sett->bgEne[0][3] = gBgBox[1]->GetX2();
        sett->bgEne[1][0] = gBgBox[2]->GetX1(); sett->bgEne[1][1] = gBgBox[2]->GetX2();
        sett->bgEne[1][2] = gBgBox[3]->GetX1(); sett->bgEne[1][3] = gBgBox[3]->GetX2();

        DrawMarkers();
        
        // Send only marker positions to UI, not full settings sync
        // This avoids overwriting checkbox states that may have changed in the UI
        // but not yet been sent back to C++
        std::string msg = "MARKER_UPDATE:";
        msg += std::to_string(sett->levEne[0]) + "|" + std::to_string(sett->levEne[1]) + "|";
        msg += std::to_string(sett->levEne[2]) + "|" + std::to_string(sett->levEne[3]) + "|";
        msg += std::to_string(sett->levEne_2[0]) + "|" + std::to_string(sett->levEne_2[1]) + "|";
        msg += std::to_string(sett->levEne_2[2]) + "|" + std::to_string(sett->levEne_2[3]) + "|";
        msg += std::to_string(sett->bgEne[0][0]) + "|" + std::to_string(sett->bgEne[0][1]) + "|";
        msg += std::to_string(sett->bgEne[0][2]) + "|" + std::to_string(sett->bgEne[0][3]) + "|";
        msg += std::to_string(sett->bgEne[1][0]) + "|" + std::to_string(sett->bgEne[1][1]) + "|";
        msg += std::to_string(sett->bgEne[1][2]) + "|" + std::to_string(sett->bgEne[1][3]);
        window->Send(0, msg);
    }
}

void RunShapeIt(unsigned connid)
{
    std::cout << "=== RunShapeIt called, display mode = " << gDisplayMode << " ===" << std::endl;
    
    if (!matrix) {
        window->Send(connid, "No matrix loaded yet -- open one first.");
        return;
    }

    std::cout << "About to run with these settings:\n";
    DumpSettings();

    std::cout << "Deleting old gSFColl..." << std::endl;
    delete gSFColl;
    std::cout << "Running analysis..." << std::endl;
    gSFColl = ShapeController::RunAnalysis(sett, matrix);
    std::cout << "Analysis complete." << std::endl;

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
        bool isAverage;
    };
    std::vector<FreshGraph> freshGraphs;
    
    // Get the literature graph pointer to reliably identify it
    // Literature data is loaded from osloFileName and only present if doOslo or doMC is true
    TGraph *litGraph = nullptr;
    if ((sett->doOslo || sett->doMC) && !sett->osloFileName.empty()) {
        litGraph = gSFColl->getLitGraph();
    }
    
    // Get the average/smoothed graph pointer to identify it
    TGraph *avgGraph = sett->displayAvg ? gSFColl->getAvgGraph() : nullptr;

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
            bool isAvg = (avgGraph != nullptr && g == avgGraph);

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
                freshGraphs.push_back({fresh, isLit, isAvg});
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
                freshGraphs.push_back({fresh, isLit, isAvg});
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

    std::cout << "About to clear canvas and draw results..." << std::endl;
    std::cout << "Current display mode before clear: " << gDisplayMode << std::endl;
    
    // If coming from width calibration view (mode 7), the TMultiGraph owns the TGraphs
    // we created, and canvas->Clear() will try to delete them. To avoid any potential
    // ownership/deletion issues, manually delete the primitives BEFORE calling Clear().
    if (gDisplayMode == 7) {
        std::cout << "Coming from width calibration view, doing explicit cleanup..." << std::endl;
        canvas->cd();
        TList *prims = canvas->GetListOfPrimitives();
        if (prims) {
            std::cout << "Canvas has " << prims->GetSize() << " primitives before clear" << std::endl;
            // Remove and delete all primitives manually to ensure clean deletion order
            while (prims->GetSize() > 0) {
                TObject *obj = prims->First();
                std::cout << "  Removing: " << obj->ClassName() << " (" << obj->GetName() << ")" << std::endl;
                prims->Remove(obj);
                delete obj;  // Explicitly delete - this will also delete owned graphs if it's a TMultiGraph
            }
            std::cout << "All primitives manually deleted." << std::endl;
        }
        // Now Clear() should have nothing to do
        std::cout << "Calling canvas->Clear() on empty canvas..." << std::endl;
        canvas->Clear();
    } else {
        // Normal case - just clear as usual
        canvas->cd();
        std::cout << "Calling canvas->Clear()..." << std::endl;
        canvas->Clear();
    }
    
    std::cout << "Canvas cleared, nulling marker pointers..." << std::endl;
    for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
    std::cout << "Marker pointers nulled." << std::endl;

    bool firstDrawn = false;
    TGraph *firstGraph = nullptr;  // Track the first graph drawn with "A" option
    int colorIdx = 0;
    // Use colors matching ShapeGSF.C: color 6 (magenta) for Level 1, color 7 (cyan) for Level 2
    // When colour=false, both use color 6
    int color1 = 6;  // kMagenta
    int color2 = sett->colour ? 7 : 6;  // kCyan if colour enabled, otherwise same as Level 1
    
    std::cout << "Drawing " << freshGraphs.size() << " graphs..." << std::endl;
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
            g->Draw(firstDrawn ? "L3 SAME" : "AL3");
            if (!firstDrawn) firstGraph = g;
        } else if (fg.isAverage) {
            // Average/smoothed graph in black
            g->SetMarkerStyle(22);
            g->SetMarkerSize(2);
            g->SetMarkerColor(1);  // black
            g->SetLineColor(1);    // black
            g->Draw(firstDrawn ? "P SAME" : "AP");
            if (!firstDrawn) firstGraph = g;
        } else {
            g->SetMarkerStyle(22);
            g->SetMarkerSize(2);
            // Alternate between Level 1 (color1) and Level 2 (color2) colors
            int color = (colorIdx % 2 == 0) ? color1 : color2;
            g->SetMarkerColor(color);
            g->SetLineColor(color);
            g->Draw(firstDrawn ? "P SAME" : "AP");
            if (!firstDrawn) firstGraph = g;
            colorIdx++;
        }
        
        firstDrawn = true;
    }
    std::cout << "Graphs drawn." << std::endl;
    
    // After all graphs are drawn, set the axis labels and plot title.
    // When a graph is drawn with the "A" option, it owns the histogram that
    // draws the axes. We can access that histogram via GetHistogram().
    if (firstGraph) {
        TH1F *hist = firstGraph->GetHistogram();
        if (hist) {
            std::cout << "Found graph histogram, setting title and axis labels..." << std::endl;
            // Set title with semicolons to separate title;xlabel;ylabel (ROOT convention)
            hist->SetTitle("Gamma Ray Strength Function from Shape Method;E_{#gamma} (keV);f(E_{#gamma}) (MeV^{-3})");
            // Also set axis labels explicitly
            hist->GetXaxis()->SetTitle("E_{#gamma} (keV)");
            hist->GetYaxis()->SetTitle("f(E_{#gamma}) (MeV^{-3})");
            // Force axis titles to be displayed
            hist->GetXaxis()->SetTitleSize(0.04);
            hist->GetYaxis()->SetTitleSize(0.04);
            hist->GetXaxis()->SetTitleOffset(1.0);
            hist->GetYaxis()->SetTitleOffset(1.2);
            // Also set the graph's own title as a fallback
            firstGraph->SetTitle("Gamma Ray Strength Function from Shape Method");
            gPad->Modified();  // Mark pad as modified after changing labels
        } else {
            std::cout << "Warning: graph histogram not found, cannot set axis labels" << std::endl;
        }
    }
    // Fallback: if no error-bar graphs were found (e.g. nothing matched
    // TGraphErrors/TGraphAsymmErrors), fall back to the plain-TGraph test
    // from before, so this still shows something.
    if (!firstGraph && !allX.empty()) {
        TGraph *simpleGraph = new TGraph((int)allX.size(), allX.data(), allY.data());
        simpleGraph->SetMarkerStyle(20);
        simpleGraph->SetMarkerColor(kBlue);
        simpleGraph->SetTitle("Gamma Ray Strength Function from Shape Method");
        simpleGraph->Draw("AP");
        
        // Set axis labels and title for the fallback graph too
        TH1F *hist = simpleGraph->GetHistogram();
        if (hist) {
            std::cout << "Found graph histogram (fallback), setting title and axis labels..." << std::endl;
            hist->SetTitle("Gamma Ray Strength Function from Shape Method;E_{#gamma} (keV);f(E_{#gamma}) (MeV^{-3})");
            hist->GetXaxis()->SetTitle("E_{#gamma} (keV)");
            hist->GetYaxis()->SetTitle("f(E_{#gamma}) (MeV^{-3})");
            hist->GetXaxis()->SetTitleSize(0.04);
            hist->GetYaxis()->SetTitleSize(0.04);
            hist->GetXaxis()->SetTitleOffset(1.0);
            hist->GetYaxis()->SetTitleOffset(1.2);
            gPad->Modified();
        } else {
            std::cout << "Warning: graph histogram not found (fallback)" << std::endl;
        }
    }

    // canvas->Clear() just deleted any marker TLine/TBox objects left over
    // from whatever projection view was showing before -- but gDisplayMode
    // was still 4/5, so the 100ms poll timer would soon call DrawMarkers()
    // again, whose first action is Remove()-ing these now-dangling pointers
    // (already nulled above, right after Clear()). That's a use-after-free
    // that plausibly corrupts the same primitive list CreatePadSnapshot walks
    // moments later -- which lines up exactly with where this crash happens.
    // This is display mode 0: "results view, no markers apply".
    std::cout << "Setting display mode to 0 (results view)..." << std::endl;
    gDisplayMode = 0;
    gHaveLastRange = false;

    std::cout << "Calling PushCanvasUpdate()..." << std::endl;
    PushCanvasUpdate();
    std::cout << "Canvas update pushed." << std::endl;

    std::cout << "Sending completion messages..." << std::endl;
    window->Send(connid, "Done.");
    SendNBins(connid);
    std::cout << "=== RunShapeIt complete ===" << std::endl;
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
    // Buffer stdout and send in batch to avoid queue overflow
    std::ostringstream logBuffer;
    std::streambuf* oldBuf = std::cout.rdbuf(logBuffer.rdbuf());
    
    // Suppress noisy debug messages for channel setup and width calibration updates
    if (arg.compare(0, 8, "channel:") != 0 && arg.compare(0, 25, "UPDATE_WIDTH_CALIB_LINES:") != 0) {
        std::cout << "Got message from browser: " << arg << std::endl;
    }

    if (arg.compare(0, 8, "channel:") == 0) {
        int chid = std::stoi(arg.substr(8));
        auto web_imp = dynamic_cast<TWebCanvas *>(canvas->GetCanvasImp());
        if (web_imp) {
            web_imp->ShowWebWindow({ window, connid, chid });
            web_imp->ForceUpdate();
        }
        window->Send(connid, "STARTDIR:" + gStartDir);
        
        // Send settings file path to frontend (must come before SETTINGS_SYNC)
        if (!sett->settFileName.empty()) {
            window->Send(connid, "SETTINGS_PATH:" + sett->settFileName);
        }
        
        // If matrix was auto-loaded at startup (via environment variable or otherwise),
        // sync it to the UI exactly as LOAD_SETTINGS does
        if (matrix) {
            auto names = matrix->GetMatrixName();
            int idx = 0;
            for (size_t i = 0; i < names.size(); i++)
                if (names[i] == sett->matrixName) idx = (int)i + 1;
            if (idx > 0) {
                gDisplayMode = 1;
                SendMatrixListAndSelect(connid, currentMatrixPath, idx);
                SendNBins(connid);
                // Enable width calibration since matrix is loaded
                window->Send(connid, "WIDTH_CALIB_AVAILABLE:1");
            }
        }
        
        // Sync all settings to UI (this updates all form fields to match loaded settings)
        SendSettingsSync(connid);
        
        // Send status message AFTER everything else is synced (same as LOAD_SETTINGS does)
        // Only send if we actually loaded settings from a file
        if (!sett->settFileName.empty()) {
            window->Send(connid, "Settings loaded: " + sett->settFileName);
        }
    }
    else if (arg.compare(0, 5, "OPEN:") == 0) {
        std::string path = arg.substr(5);

        if (gSystem->AccessPathName(path.c_str())) {
            // AccessPathName returns non-zero (true) when the path does NOT exist
            window->Send(connid, "File not found: " + path);
            return;
        }

        // Check if this is actually a ROOT file before attempting to open it
        // TFile::Open() will print errors to stdout but won't throw exceptions
        TFile *testFile = TFile::Open(path.c_str(), "READ");
        if (!testFile || testFile->IsZombie()) {
            if (testFile) delete testFile;
            window->Send(connid, "ERROR: Not a valid ROOT file or file is corrupted: " + path);
            std::cout << "User attempted to open non-ROOT file as matrix: " << path << std::endl;
            return;
        }
        testFile->Close();
        delete testFile;

        sett->SetFileName(path);
        currentMatrixPath = path;
        delete matrix;
        matrix = nullptr;
        matrix = new ShapeMatrix(sett);
        SendMatrixListAndSelect(connid, currentMatrixPath, 1);
        SendNBins(connid);
        // Enable width calibration now that we have a matrix loaded
        window->Send(connid, "WIDTH_CALIB_AVAILABLE:1");
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
        // Enable width calibration since we have a matrix
        window->Send(connid, "WIDTH_CALIB_AVAILABLE:1");
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
                // Enable width calibration since we have a matrix loaded
                window->Send(connid, "WIDTH_CALIB_AVAILABLE:1");
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
        // order: doInterpol|doOslo|doSlidingWindow|doBackground|doWidthCal
        auto v = ParsePipeDoubles(arg.substr(8));
        if (v.size() != 5) {
            window->Send(connid, "Malformed OPTIONS message.");
            return;
        }
        sett->doInterpol      = v[0] != 0.0;
        sett->doOslo          = v[1] != 0.0;
        sett->doSlidingWindow = v[2] != 0.0;
        sett->doBackground    = v[3] != 0.0;
        sett->doWidthCal      = v[4] != 0.0;
    }
    else if (arg.compare(0, 16, "DISPLAY_OPTIONS:") == 0) {
        // order: displaySingle|displayAvg|colour
        auto v = ParsePipeDoubles(arg.substr(16));
        if (v.size() != 3) {
            window->Send(connid, "Malformed DISPLAY_OPTIONS message.");
            return;
        }
        sett->displaySingle = v[0] != 0.0;
        sett->displayAvg    = v[1] != 0.0;
        sett->colour        = v[2] != 0.0;
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
    }
    else if (arg.compare(0, 8, "VERBOSE:") == 0) {
        sett->verbose = std::stoi(arg.substr(8));
    }
    else if (arg.compare(0, 5, "MODE:") == 0) {
        sett->mode = std::stoi(arg.substr(5)); // 1 = Integration, 2 = Autofit
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
    }
    else if (arg == "UPDATE_MARKERS") {
        // Redraw markers with current settings (triggered by energy changes in UI)
        DrawMarkers(true);
    }
    else if (arg == "SHOW_LEVELS_PANEL") {
        // When Levels panel is clicked, show bin 1 projection ONLY if not already viewing a projection
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet -- open one first.");
            return;
        }
        
        // Only switch to bin 1 projection if not already in a projection view (modes 4 or 5)
        if (gDisplayMode != 4 && gDisplayMode != 5) {
            std::cout << "Levels panel opened - switching from mode " << gDisplayMode << " to bin 1 projection" << std::endl;
            gCurrentBin = 1;
            gDisplayMode = 5;
            canvas->cd();
            canvas->Clear();
            for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
            gCurrentHist = matrix->GetDiagEx(gCurrentBin, BaseName(currentMatrixPath));
            gCurrentHist->Draw();
            gHaveLastRange = false;
            CleanupAutofitDisplay();
            DrawMarkers();
        } else {
            std::cout << "Levels panel opened - already in projection mode " << gDisplayMode << ", keeping current view" << std::endl;
        }
    }
    else if (arg == "SHOWMATRIX") {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet -- open one first.");
            return;
        }
        gDisplayMode = 1;
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
        
        // Set pad margins to accommodate the manually-positioned color palette
        gPad->SetRightMargin(0.15);  // 15% right margin provides space for palette
        gPad->SetLeftMargin(0.12);   // 12% left margin for y-axis label
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.10);
        
        TH2* hist = matrix->GetInputMatrix(BaseName(currentMatrixPath));
        hist->SetStats(0);  // Disable statistics box
        
        // Set axis titles with LaTeX formatting
        hist->GetXaxis()->SetTitle("E_{#gamma} (keV)");
        hist->GetYaxis()->SetTitle("E_{x} (keV)");
        hist->GetXaxis()->SetTitleSize(0.045);
        hist->GetYaxis()->SetTitleSize(0.045);
        hist->GetXaxis()->SetTitleOffset(1.0);
        hist->GetYaxis()->SetTitleOffset(1.1);
        
        // Draw histogram without automatic palette (use "col" not "colz")
        hist->Draw("col");
        
        // Manually create and position the color palette using NDC coordinates
        // First create it with histogram coordinates (required by constructor)
        double xmin = hist->GetXaxis()->GetXmin();
        double xmax = hist->GetXaxis()->GetXmax();
        double ymin = hist->GetYaxis()->GetXmin();
        double ymax = hist->GetYaxis()->GetXmax();
        TPaletteAxis *palette = new TPaletteAxis(xmax, ymin, xmax + (xmax-xmin)*0.05, ymax, hist);
        
        // Override with NDC coordinates to position it within the right margin
        palette->SetX1NDC(0.86);  // Left edge at 86% of canvas width
        palette->SetX2NDC(0.89);  // Right edge at 89% of canvas width
        palette->SetY1NDC(0.10);  // Bottom aligned with pad margin
        palette->SetY2NDC(0.90);  // Top aligned with pad margin
        palette->Draw();
        
        DrawMarkers();
        PushCanvasUpdate();
    }
    else if (arg == "SHOWPROJ") {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet -- open one first.");
            return;
        }
        gDisplayMode = 4;
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
        gCurrentHist = matrix->GetDiag(BaseName(currentMatrixPath));
        gCurrentHist->Draw("hist");
        gHaveLastRange = false;
        DrawMarkers();
    }
    else if (arg == "SHOW_WIDTH_CALIB") {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet.");
            return;
        }
        
        // Save existing calibration parameters before running fresh fit
        double savedWidthCal[2][2];
        savedWidthCal[0][0] = sett->widthCal[0][0];
        savedWidthCal[0][1] = sett->widthCal[0][1];
        savedWidthCal[1][0] = sett->widthCal[1][0];
        savedWidthCal[1][1] = sett->widthCal[1][1];
        
        bool hasExistingCalib = (savedWidthCal[0][0] != 0.0 || savedWidthCal[0][1] != 0.0 ||
                                  savedWidthCal[1][0] != 0.0 || savedWidthCal[1][1] != 0.0);
        
        // ALWAYS run fresh fit to generate graph data points AND fit parameters
        window->Send(connid, "Generating width calibration data...");
        RunWidthCalibration(connid);
        
        // Extract the width data graphs (populated by RunWidthCalibration above)
        TGraph *T1 = matrix->getFitWidthGraph(0);
        TGraph *T2 = matrix->getFitWidthGraph(1);
        
        // Perform linear fits on the width data (mirrors ShapeFrame.C case 7)
        // This populates sett->widthCal with the fit parameters
        if (T1->GetN() > 0) {
            T1->Fit("pol1", "Q");  // Q = quiet mode
            TF1 *fit1 = T1->GetFunction("pol1");
            if (fit1) {
                sett->widthCal[0][0] = fit1->GetParameter(0);
                sett->widthCal[0][1] = fit1->GetParameter(1);
            }
        }
        
        if (T2->GetN() > 0) {
            T2->Fit("pol1", "Q");  // Q = quiet mode
            TF1 *fit2 = T2->GetFunction("pol1");
            if (fit2) {
                sett->widthCal[1][0] = fit2->GetParameter(0);
                sett->widthCal[1][1] = fit2->GetParameter(1);
            }
        }
        
        // Only restore old fit parameters if they were non-zero (i.e., from a settings file with width cal active)
        // If all zeros, keep the fresh fit results instead
        if (hasExistingCalib) {
            sett->widthCal[0][0] = savedWidthCal[0][0];
            sett->widthCal[0][1] = savedWidthCal[0][1];
            sett->widthCal[1][0] = savedWidthCal[1][0];
            sett->widthCal[1][1] = savedWidthCal[1][1];
        }
        // else: keep the fresh fit parameters from the fits above
        
        canvas->cd();
        canvas->Clear();
        // Null out marker pointers since Clear() deleted them
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
        gDisplayMode = 7;  // Width calibration display mode
        
        // Create fit functions using parameters from sett->widthCal
        // These will either be from the file or from the fresh fit we just ran
        TF1 *fit1 = nullptr;
        TF1 *fit2 = nullptr;
        
        if (T1->GetN() > 0) {
            fit1 = new TF1("fit1", "[0] + [1]*x", 0, 10000);
            fit1->SetParameter(0, sett->widthCal[0][0]);
            fit1->SetParameter(1, sett->widthCal[0][1]);
            fit1->SetLineColor(kRed);
            fit1->SetLineWidth(2);
        }
        
        if (T2->GetN() > 0) {
            fit2 = new TF1("fit2", "[0] + [1]*x", 0, 10000);
            fit2->SetParameter(0, sett->widthCal[1][0]);
            fit2->SetParameter(1, sett->widthCal[1][1]);
            fit2->SetLineColor(kBlue);
            fit2->SetLineWidth(2);
        }
        
        // Find the combined range of both graphs for proper axis scaling
        double xMin = 1e9, xMax = -1e9, yMin = 1e9, yMax = -1e9;
        
        if (T1->GetN() > 0) {
            double xmin1, xmax1, ymin1, ymax1;
            T1->ComputeRange(xmin1, ymin1, xmax1, ymax1);
            xMin = std::min(xMin, xmin1);
            xMax = std::max(xMax, xmax1);
            yMin = std::min(yMin, ymin1);
            yMax = std::max(yMax, ymax1);
        }
        
        if (T2->GetN() > 0) {
            double xmin2, xmax2, ymin2, ymax2;
            T2->ComputeRange(xmin2, ymin2, xmax2, ymax2);
            xMin = std::min(xMin, xmin2);
            xMax = std::max(xMax, xmax2);
            yMin = std::min(yMin, ymin2);
            yMax = std::max(yMax, ymax2);
        }
        
        // Add 5% padding to ranges
        double xRange = xMax - xMin;
        double yRange = yMax - yMin;
        xMin -= 0.05 * xRange;
        xMax += 0.05 * xRange;
        yMin -= 0.05 * yRange;
        yMax += 0.05 * yRange;
        
        // Save these ranges for future updates
        gWidthCalibXMin = xMin;
        gWidthCalibXMax = xMax;
        gWidthCalibYMin = yMin;
        gWidthCalibYMax = yMax;
        gHaveWidthCalibRanges = true;
        
        // Draw graphs directly without TMultiGraph to avoid any potential Update() calls
        // Style the graphs
        T1->SetMarkerStyle(20);
        T1->SetMarkerColor(kRed);
        T1->SetLineColor(kRed);
        T1->SetMarkerSize(1);
        
        T2->SetMarkerStyle(21);
        T2->SetMarkerColor(kBlue);
        T2->SetLineColor(kBlue);
        T2->SetMarkerSize(1);
        
        // Draw first graph with axes - use explicit range
        T1->Draw("AP");
        T1->SetTitle("Peak Widths from Autofit");
        T1->GetXaxis()->SetTitle("E_{#gamma} (keV)");
        T1->GetYaxis()->SetTitle("Width (keV)");
        
        // Set the axis ranges to include both graphs
        T1->GetXaxis()->SetLimits(xMin, xMax);
        T1->GetHistogram()->SetMinimum(yMin);
        T1->GetHistogram()->SetMaximum(yMax);
        
        // Draw second graph on top
        if (T2->GetN() > 0) {
            T2->Draw("P SAME");
        }
        
        // Draw fit lines manually - extend fit range to cover the full axis range
        if (fit1) {
            fit1->SetRange(xMin, xMax);
            fit1->Draw("SAME");
        }
        if (fit2) {
            fit2->SetRange(xMin, xMax);
            fit2->Draw("SAME");
        }
        
        // Add legend (smaller size)
        TLegend *leg = new TLegend(0.75, 0.80, 0.90, 0.90);
        leg->SetFillColor(0);
        if (T1->GetN() > 0) leg->AddEntry(T1, "level 1", "lp");
        if (T2->GetN() > 0) leg->AddEntry(T2, "level 2", "lp");
        leg->Draw();
        
        // Enable the width calibration checkbox by telling frontend it's available
        window->Send(connid, "WIDTH_CALIB_AVAILABLE:1");
        
        // Send the final calibration parameters to UI (either restored from file or fresh fit)
        std::string msg = "WIDTH_CALIB_PARAMS:";
        msg += std::to_string(sett->widthCal[0][0]) + "|" + std::to_string(sett->widthCal[0][1]) + "|";
        msg += std::to_string(sett->widthCal[1][0]) + "|" + std::to_string(sett->widthCal[1][1]);
        window->Send(connid, msg);
        
        // DON'T delete T1/T2 here - they're now drawn on the canvas and will be
        // cleaned up automatically when the canvas is cleared.
        
        PushCanvasUpdate();
    }
    else if (arg == "SHOW_SETTINGS_FILE") {
        if (sett->settFileName.empty()) {
            window->Send(connid, "No settings file loaded.");
            return;
        }
        
        std::ifstream file(sett->settFileName.c_str());
        if (!file.good()) {
            window->Send(connid, "Could not read settings file: " + sett->settFileName);
            return;
        }
        
        std::stringstream buffer;
        buffer << file.rdbuf();
        file.close();
        
        std::string content = buffer.str();
        window->Send(connid, "SETTINGS_FILE_CONTENT:" + content);
    }
    else if (arg == "SHOW_LEVEL_DENSITY") {
        if (sett->settFileName.empty()) {
            window->Send(connid, "No settings file loaded.");
            return;
        }
        
        // Check if level density file exists in settings
        if (sett->rhoFileName.empty()) {
            window->Send(connid, "No level density file specified in settings.");
            return;
        }
        
        // Resolve relative path against settings file directory
        std::string settDir = DirName(sett->settFileName);
        std::string rhoPath = ResolveRelativeTo(settDir, sett->rhoFileName);
        
        // Check if file exists
        if (gSystem->AccessPathName(rhoPath.c_str())) {
            window->Send(connid, "Level density file not found: " + rhoPath);
            return;
        }
        
        // Read level density data from file (format: Ex rho rho_error)
        TGraphErrors *grLevelDensity = new TGraphErrors(rhoPath.c_str(), "%lg %lg %lg");
        
        if (grLevelDensity->GetN() == 0) {
            delete grLevelDensity;
            window->Send(connid, "No level density data in file: " + rhoPath);
            return;
        }
        
        // Apply scaling factor if set
        if (sett->rhoScale != 1.0) {
            for (int i = 0; i < grLevelDensity->GetN(); i++) {
                grLevelDensity->SetPoint(i, grLevelDensity->GetX()[i], 
                                         sett->rhoScale * grLevelDensity->GetY()[i]);
            }
        }
        
        grLevelDensity->SetMarkerStyle(20);
        grLevelDensity->SetMarkerSize(0.8);
        grLevelDensity->SetMarkerColor(kBlue);
        grLevelDensity->SetLineColor(kBlue);
        grLevelDensity->SetTitle("Level Density");
        
        gDisplayMode = 0; // Not a mode with markers
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
        gCurrentHist = nullptr;
        
        // Draw level density graph first
        grLevelDensity->Draw("APE");
        grLevelDensity->GetXaxis()->SetTitle("Excitation Energy (keV)");
        grLevelDensity->GetYaxis()->SetTitle("Level Density (1/keV)");
        
        // Read and overlay discrete levels histogram if available
        TH1F *discreteLevel = nullptr;
        if (!sett->discreteLevelFile.empty()) {
            std::string discPath = ResolveRelativeTo(settDir, sett->discreteLevelFile);
            
            if (!gSystem->AccessPathName(discPath.c_str())) {
                std::ifstream discFile(discPath.c_str());
                if (discFile.good()) {
                    std::vector<double> ene, discLev;
                    double e, disc;
                    
                    // Read discrete level data (format: energy level_density)
                    while (discFile >> e >> disc) {
                        ene.push_back(e);
                        discLev.push_back(disc);
                    }
                    discFile.close();
                    
                    if (!ene.empty()) {
                        // Find maximum energy
                        double discreteMax = *std::max_element(ene.begin(), ene.end());
                        
                        // Create histogram with bin size from settings
                        int nBins = (int)((1000.0 * discreteMax) / sett->discreteBins);
                        discreteLevel = new TH1F("discreteLevel", "discrete levels", nBins, 0, discreteMax);
                        
                        // Fill histogram
                        for (size_t i = 0; i < ene.size(); i++) {
                            discreteLevel->Fill(ene[i], discLev[i]);
                        }
                        
                        // Style to match ShapeRhoCollector
                        discreteLevel->SetFillColorAlpha(kAzure-9, 0.4);
                        discreteLevel->SetFillStyle(3002);
                        discreteLevel->SetLineColorAlpha(kBlack, 0.6);
                        
                        // Draw on same canvas
                        discreteLevel->Draw("same hist");
                        
                        window->Send(connid, "Displaying level density (" + std::to_string(grLevelDensity->GetN()) + 
                                     " points) with discrete levels (" + std::to_string(ene.size()) + " levels)");
                    }
                }
            }
        }
        
        if (!discreteLevel) {
            window->Send(connid, "Displaying level density: " + std::to_string(grLevelDensity->GetN()) + " points");
        }
        
        PushCanvasUpdate();
    }
    else if (arg == "SHOW_GSF_RESULTS") {
        std::cout << "=== SHOW_GSF_RESULTS handler called ===" << std::endl;
        
        if (!gSFColl) {
            std::cout << "ERROR: gSFColl is null" << std::endl;
            window->Send(connid, "GSF_RESULTS:ERROR:No gSF results available. Run ShapeIt first.");
            return;
        }
        
        std::cout << "gSFColl exists, checking display options..." << std::endl;
        std::cout << "  displayAvg: " << sett->displayAvg << std::endl;
        std::cout << "  displaySingle: " << sett->displaySingle << std::endl;
        std::cout << "  GetNSmooth(): " << gSFColl->GetNSmooth() << std::endl;
        std::cout << "  GetN(): " << gSFColl->GetN() << std::endl;
        
        // Check if we have anything to display
        bool hasAverage = sett->displayAvg && (gSFColl->GetNSmooth() > 0);
        bool hasIndividual = sett->displaySingle && (gSFColl->GetN() > 0);
        
        if (!hasAverage && !hasIndividual) {
            std::cout << "ERROR: No display options enabled" << std::endl;
            window->Send(connid, "GSF_RESULTS:ERROR:No gSF display options selected. Enable 'Display Average' or 'Display Individual' in Options panel.");
            return;
        }
        
        // Build gSF results text by directly extracting data from the collector graphs
        std::ostringstream resultText;
        
        // Show smoothed/average graph if enabled
        if (sett->displayAvg && hasAverage) {
            std::cout << "Extracting average graph data..." << std::endl;
            TGraphAsymmErrors *smoothGraph = gSFColl->getAvgGraph();
            if (smoothGraph && smoothGraph->GetN() > 0) {
                std::cout << "  Average graph has " << smoothGraph->GetN() << " points" << std::endl;
                resultText << "gSF values for smoothed graph:\n";
                resultText << "energy    gSF   error gSF\n";
                
                for (int i = 0; i < smoothGraph->GetN(); i++) {
                    double e = smoothGraph->GetX()[i];
                    double g = smoothGraph->GetY()[i];
                    double dgHigh = smoothGraph->GetEYhigh()[i];
                    double dgLow = smoothGraph->GetEYlow()[i];
                    
                    if (std::abs(dgHigh - dgLow) < 1e-10)
                        resultText << e << " " << g << " " << dgHigh << "\n";
                    else
                        resultText << e << " " << g << " + " << dgHigh << " - " << dgLow << "\n";
                }
                resultText << "\n";
            } else {
                std::cout << "WARNING: Average graph is null or empty" << std::endl;
            }
        }
        
        // Show individual merged data points if enabled (sorted by energy)
        if (sett->displaySingle && hasIndividual) {
            std::cout << "Extracting individual merged data points..." << std::endl;
            TGraphErrors *mergedGraph = gSFColl->getMergedGraph();
            if (mergedGraph && mergedGraph->GetN() > 0) {
                std::cout << "  Merged graph has " << mergedGraph->GetN() << " points" << std::endl;
                resultText << "Individual gSF data points (all iterations merged, sorted by energy):\n";
                resultText << "energy    gSF   error\n";
                
                for (int i = 0; i < mergedGraph->GetN(); i++) {
                    double e = mergedGraph->GetX()[i];
                    double g = mergedGraph->GetY()[i];
                    double dg = mergedGraph->GetEY()[i];
                    
                    resultText << e << " " << g << " " << dg << "\n";
                }
                resultText << "\n";
            } else {
                std::cout << "WARNING: Merged graph is null or empty" << std::endl;
            }
        }
        
        // Show literature data if loaded
        if (sett->doOslo) {
            std::cout << "Extracting literature graph data..." << std::endl;
            TGraphErrors *litGraph = gSFColl->getLitGraph();
            if (litGraph && litGraph->GetN() > 0) {
                std::cout << "  Literature graph has " << litGraph->GetN() << " points" << std::endl;
                resultText << "Literature gSF values:\n";
                resultText << "energy    gSF    error\n";
                
                for (int i = 0; i < litGraph->GetN(); i++) {
                    double e = litGraph->GetX()[i];
                    double g = litGraph->GetY()[i];
                    double dg = litGraph->GetEY()[i];
                    
                    resultText << e << " " << g << " " << dg << "\n";
                }
            } else {
                std::cout << "  Literature graph is null or empty" << std::endl;
            }
        }
        
        std::string finalText = resultText.str();
        std::cout << "Generated " << finalText.length() << " characters of output" << std::endl;
        std::cout << "First 200 chars: " << finalText.substr(0, std::min((size_t)200, finalText.length())) << std::endl;
        
        window->Send(connid, "GSF_RESULTS:SUCCESS:" + finalText);
        std::cout << "=== SHOW_GSF_RESULTS complete ===" << std::endl;
    }
    else if (arg.compare(0, 12, "SHOWBINPROJ:") == 0) {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet -- open one first.");
            return;
        }
        
        // Save only X-axis zoom state before switching bins
        // Let Y-axis auto-scale for each new bin (different bins have different count ranges)
        bool hadRange = gHaveLastRange;
        double savedXmin = gLastUxmin;
        double savedXmax = gLastUxmax;
        
        gCurrentBin = std::stoi(arg.substr(12));
        gDisplayMode = 5;
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
        gCurrentHist = matrix->GetDiagEx(gCurrentBin, BaseName(currentMatrixPath));
        gCurrentHist->Draw();
        
        // Restore X-axis zoom only after drawing new histogram
        if (hadRange && gPad) {
            gCurrentHist->GetXaxis()->SetRangeUser(savedXmin, savedXmax);
            // Y-axis is left at its auto-scaled default
            gPad->Modified();
            
            // Update tracking variables with current state
            gLastUxmin = savedXmin;
            gLastUxmax = savedXmax;
            gLastUymin = gPad->GetUymin();
            gLastUymax = gPad->GetUymax();
            gHaveLastRange = true;
        } else {
            gHaveLastRange = false;
        }
        
        CleanupAutofitDisplay();
        DrawMarkers(hadRange);
    }
    else if (arg.compare(0, 11, "EXCITATION:") == 0) {
        auto v = ParsePipeDoubles(arg.substr(11));
        if (v.size() == 2) {
            sett->exiEne[0] = v[0];
            sett->exiEne[1] = v[1];
            SendNBins(connid);
        }
    }
    else if (arg.compare(0, 14, "LEVELENERGIES:") == 0) {
        std::cout << "*** LEVELENERGIES HANDLER CALLED ***" << std::endl;
        auto v = ParsePipeDoubles(arg.substr(14));
        std::cout << "*** Parsed " << v.size() << " values ***" << std::endl;
        if (v.size() != 12) {
            std::cout << "*** ERROR: Expected 12 values, got " << v.size() << " ***" << std::endl;
            return;
        }
        std::cout << "*** Setting levEne[0-3] to: " << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << " ***" << std::endl;
        
        // Store previous doublet checkbox states to detect changes
        bool hadDoublet1 = sett->doDoublet[0];
        bool hadDoublet2 = sett->doDoublet[1];
        bool hadFixWidth1 = sett->fixDoubletWidth[0];
        bool hadFixWidth2 = sett->fixDoubletWidth[1];
        
        // Update main peak energies
        sett->levEne[0] = v[0]; 
        sett->levEne[1] = v[1];
        sett->levEne[2] = v[2]; 
        sett->levEne[3] = v[3];
        
        // Update doublet checkbox states
        sett->doDoublet[0] = v[4] != 0.0;
        sett->doDoublet[1] = v[7] != 0.0;
        
        // Update doublet width fix toggles
        sett->fixDoubletWidth[0] = v[10] != 0.0;
        sett->fixDoubletWidth[1] = v[11] != 0.0;
        
        // ALWAYS update doublet energy values regardless of checkbox state
        sett->levEne_2[0] = v[5];
        sett->levEne_2[1] = v[6];
        sett->levEne_2[2] = v[8];
        sett->levEne_2[3] = v[9];
        
        std::cout << "*** Level energies updated successfully ***" << std::endl;
        std::cout << "*** Doublet 1: " << (sett->doDoublet[0] ? "ENABLED" : "DISABLED") 
                  << ", values: " << sett->levEne_2[0] << ", " << sett->levEne_2[1] 
                  << ", fix width: " << (sett->fixDoubletWidth[0] ? "YES" : "NO") << " ***" << std::endl;
        std::cout << "*** Doublet 2: " << (sett->doDoublet[1] ? "ENABLED" : "DISABLED") 
                  << ", values: " << sett->levEne_2[2] << ", " << sett->levEne_2[3] 
                  << ", fix width: " << (sett->fixDoubletWidth[1] ? "YES" : "NO") << " ***" << std::endl;
        
        // Detect if doublet state or width fix state changed
        bool doubletStateChanged = (hadDoublet1 != sett->doDoublet[0]) || (hadDoublet2 != sett->doDoublet[1]);
        bool widthFixChanged = (hadFixWidth1 != sett->fixDoubletWidth[0]) || (hadFixWidth2 != sett->fixDoubletWidth[1]);
        
        // If in Autofit mode and viewing a bin projection, re-fit when doublet checkbox or width fix changes
        if ((doubletStateChanged || widthFixChanged) && sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
            std::cout << "Doublet settings changed in Autofit mode: re-fitting bin " << gCurrentBin << "..." << std::endl;
            
            // Save current axis ranges before redrawing
            double xmin = gPad->GetUxmin();
            double xmax = gPad->GetUxmax();
            double ymin = gPad->GetUymin();
            double ymax = gPad->GetUymax();
            bool isLogy = gPad->GetLogy();
            
            canvas->cd();
            gCurrentHist = matrix->GetDiagEx(gCurrentBin, BaseName(currentMatrixPath));
            
            // Restore axis ranges
            gCurrentHist->GetXaxis()->SetRangeUser(xmin, xmax);
            if (isLogy)
                gCurrentHist->GetYaxis()->SetRangeUser(TMath::Power(10, ymin), TMath::Power(10, ymax));
            else
                gCurrentHist->GetYaxis()->SetRangeUser(ymin, ymax);
            
            gCurrentHist->Draw();
            CleanupAutofitDisplay();
            DrawMarkers(true);
            PushCanvasUpdate();
        } else {
            // Just redraw markers (handles Integration mode or when not viewing projection)
            DrawMarkers(true);
        }
    }
    else if (arg.compare(0, 19, "WIDTH_CALIB_PARAMS:") == 0) {
        // order: l1_offset|l1_slope|l2_offset|l2_slope
        auto v = ParsePipeDoubles(arg.substr(19));
        if (v.size() != 4) {
            window->Send(connid, "Malformed WIDTH_CALIB_PARAMS message.");
            return;
        }
        sett->widthCal[0][0] = v[0];  // Level 1 offset
        sett->widthCal[0][1] = v[1];  // Level 1 slope
        sett->widthCal[1][0] = v[2];  // Level 2 offset
        sett->widthCal[1][1] = v[3];  // Level 2 slope
    }
    else if (arg.compare(0, 25, "UPDATE_WIDTH_CALIB_LINES:") == 0) {
        // Just update the fit line parameters without re-running analysis
        auto v = ParsePipeDoubles(arg.substr(25));
        if (v.size() != 4) {
            window->Send(connid, "Malformed UPDATE_WIDTH_CALIB_LINES message.");
            return;
        }
        sett->widthCal[0][0] = v[0];
        sett->widthCal[0][1] = v[1];
        sett->widthCal[1][0] = v[2];
        sett->widthCal[1][1] = v[3];
        
        // Only redraw if we're currently viewing width calibration (mode 7)
        if (gDisplayMode == 7 && matrix) {
            // Use the saved ranges from when the plot was first created
            if (!gHaveWidthCalibRanges) {
                return;
            }
            
            TGraph *T1 = matrix->getFitWidthGraph(0);
            TGraph *T2 = matrix->getFitWidthGraph(1);
            
            // Clear canvas and redraw everything from scratch with correct ranges
            canvas->cd();
            canvas->Clear();
            
            // Style the graphs
            T1->SetMarkerStyle(20);
            T1->SetMarkerColor(kRed);
            T1->SetLineColor(kRed);
            T1->SetMarkerSize(1);
            
            T2->SetMarkerStyle(21);
            T2->SetMarkerColor(kBlue);
            T2->SetLineColor(kBlue);
            T2->SetMarkerSize(1);
            
            // Draw first graph with axes - force the saved range
            T1->Draw("AP");
            T1->SetTitle("Peak Widths from Autofit");
            T1->GetXaxis()->SetTitle("E_{#gamma} (keV)");
            T1->GetYaxis()->SetTitle("Width (keV)");
            T1->GetXaxis()->SetLimits(gWidthCalibXMin, gWidthCalibXMax);
            T1->GetHistogram()->SetMinimum(gWidthCalibYMin);
            T1->GetHistogram()->SetMaximum(gWidthCalibYMax);
            
            // Draw second graph
            if (T2->GetN() > 0) {
                T2->Draw("P SAME");
            }
            
            // Now create and draw fit functions with the extended range
            if (T1->GetN() > 0) {
                TF1 *fit1 = new TF1("fit1", "[0] + [1]*x", gWidthCalibXMin, gWidthCalibXMax);
                fit1->SetParameter(0, sett->widthCal[0][0]);
                fit1->SetParameter(1, sett->widthCal[0][1]);
                fit1->SetLineColor(kRed);
                fit1->SetLineWidth(2);
                fit1->Draw("SAME");
            }
            
            if (T2->GetN() > 0) {
                TF1 *fit2 = new TF1("fit2", "[0] + [1]*x", gWidthCalibXMin, gWidthCalibXMax);
                fit2->SetParameter(0, sett->widthCal[1][0]);
                fit2->SetParameter(1, sett->widthCal[1][1]);
                fit2->SetLineColor(kBlue);
                fit2->SetLineWidth(2);
                fit2->Draw("SAME");
            }
            
            // Add legend
            TLegend *leg = new TLegend(0.75, 0.80, 0.90, 0.90);
            leg->SetFillColor(0);
            if (T1->GetN() > 0) leg->AddEntry(T1, "level 1", "lp");
            if (T2->GetN() > 0) leg->AddEntry(T2, "level 2", "lp");
            leg->Draw();
            
            PushCanvasUpdate();
        }
    }
    else if (arg.compare(0, 17, "SAVE_WIDTH_CALIB:") == 0) {
        std::string path = arg.substr(17);
        sett->settFileName = path;
        sett->SaveSettings();
        window->Send(connid, "Width calibration saved to: " + path);
    }
    else if (arg.compare(0, 4, "RUN:") == 0) {
        // expected order: lvl1_lo|lvl1_hi|lvl2_lo|lvl2_hi|exc_lo|exc_hi|
        //                  is_doublet1|d1_lo|d1_hi|is_doublet2|d2_lo|d2_hi|fix_width1|fix_width2
        auto v = ParsePipeDoubles(arg.substr(4));
        if (v.size() != 14) {
            window->Send(connid, "Malformed RUN message.");
            return;
        }

        sett->levEne[0] = v[0]; sett->levEne[1] = v[1];
        sett->levEne[2] = v[2]; sett->levEne[3] = v[3];
        sett->exiEne[0] = v[4]; sett->exiEne[1] = v[5];

        // Set doublet checkbox states
        sett->doDoublet[0] = v[6] != 0.0;
        sett->doDoublet[1] = v[9] != 0.0;
        
        // ALWAYS store doublet energy values regardless of checkbox state
        sett->levEne_2[0] = v[7];
        sett->levEne_2[1] = v[8];
        sett->levEne_2[2] = v[10];
        sett->levEne_2[3] = v[11];
        
        // Set doublet width fix toggles
        sett->fixDoubletWidth[0] = v[12] != 0.0;
        sett->fixDoubletWidth[1] = v[13] != 0.0;

        // Keep stdout redirected to logBuffer so verbose output is captured
        std::cout << "*** About to call RunShapeIt() ***" << std::endl;
        
        RunShapeIt(connid);
        
        std::cout << "*** RunShapeIt() returned ***" << std::endl;
        
        // Restore stdout and send all captured output including verbose logs
        std::cout.rdbuf(oldBuf);
        std::string logs = logBuffer.str();
        if (!logs.empty() && window) {
            window->Send(connid, "LOGBATCH:" + logs);
        }
        
        // Return early since we've already handled stdout restoration
        return;
    }
    
    // Restore stdout and send batched log
    std::cout.rdbuf(oldBuf);
    std::string logs = logBuffer.str();
    if (!logs.empty() && window) {
        window->Send(connid, "LOGBATCH:" + logs);
    }
}

void WebShapeIt()
{
    // Check for help flag first
    if (GetCmdLineArg("--help") == "--help" || GetCmdLineArg("-h") == "-h") {
        std::cout << "\nWebShapeIt 2.0 - Usage:\n"
                  << "  root -l WebShapeIt.cxx\n"
                  << "  SHAPEIT_SETTINGS=<path> root -l WebShapeIt.cxx\n\n"
                  << "Options:\n"
                  << "  SHAPEIT_SETTINGS   Environment variable to load settings file at startup\n"
                  << "  --help, -h         Show this help message\n\n"
                  << "Examples:\n"
                  << "  # Load specific settings file:\n"
                  << "  SHAPEIT_SETTINGS=../Analysis/88Kr/test.dat root -l WebShapeIt.cxx\n\n"
                  << "  # Start with empty settings:\n"
                  << "  root -l WebShapeIt.cxx\n"
                  << std::endl;
        return;  // Exit without starting the GUI
    }
    
    gEnv->SetValue("WebGui.ConnCredits", "100");
    gStartDir = gSystem->WorkingDirectory();

    // Configure stat box size (default is ~0.3x0.2 NDC units, reduce by factor of 2)
    gStyle->SetStatW(0.15);  // width: 0.3 -> 0.15
    gStyle->SetStatH(0.10);  // height: 0.2 -> 0.10

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

    // Check for command-line --settings argument (via environment variable)
    std::string cmdLineSettings = GetCmdLineArg("--settings");
    
    if (!cmdLineSettings.empty()) {
        std::ifstream testFile(cmdLineSettings.c_str());
        if (testFile.good()) {
            testFile.close();
            
            // Suppress verbose ReadSettings() output
            std::ostringstream devnull;
            std::streambuf* oldBuf = std::cout.rdbuf(devnull.rdbuf());
            
            sett->settFileName = cmdLineSettings;
            sett->ReadSettings();
            
            // Restore stdout
            std::cout.rdbuf(oldBuf);
            
            // Resolve relative paths in settings file
            std::string settDir = DirName(cmdLineSettings);
            sett->dataFileName = ResolveRelativeTo(settDir, sett->dataFileName);
            sett->osloFileName = ResolveRelativeTo(settDir, sett->osloFileName);
            
            std::cout << "Loaded settings from: " << cmdLineSettings << std::endl;
            
            // Auto-load the matrix file referenced in settings (same as LOAD_SETTINGS handler)
            if (!sett->dataFileName.empty() && !gSystem->AccessPathName(sett->dataFileName.c_str())) {
                currentMatrixPath = sett->dataFileName;
                matrix = new ShapeMatrix(sett);
                
                // Find the correct matrix by name (same as LOAD_SETTINGS handler)
                auto names = matrix->GetMatrixName();
                int idx = 0;
                for (size_t i = 0; i < names.size(); i++)
                    if (names[i] == sett->matrixName) idx = (int)i + 1;
                
                if (idx > 0) {
                    // Set the matrix (this will be drawn when UI connects)
                    matrix->SetMatrix(idx);
                } else {
                    std::cout << "Warning: matrix '" << sett->matrixName 
                              << "' not found in " << sett->dataFileName << std::endl;
                }
            } else if (!sett->dataFileName.empty()) {
                std::cout << "Warning: matrix file not found: " << sett->dataFileName << std::endl;
            }
        } else {
            std::cout << "Warning: settings file not found: " << cmdLineSettings << std::endl;
        }
    } else {
        // No environment variable set - start with empty settings (no output needed)
    }

    canvas = TWebCanvas::CreateWebCanvas("webshapeit_canvas", "ShapeIt 2.0");
    
    // Enable crosshair and coordinate display - shows x,y values as mouse moves
    // kCrosshair = 1: both vertical and horizontal lines
    canvas->SetCrosshair(1);
    
    // Force the status bar to be shown (displays coordinates)
    canvas->ToggleEventStatus();
    canvas->SetBit(TCanvas::kShowEventStatus);

    // Connects to a plain global function (no custom dictionary-registered
    // class needed) -- this is what lets DrawMarkers()/HandleCanvasEvent() know
    // when a dragged marker/box has been released, to read its new position
    // back into sett.
    canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
                     "HandleCanvasEvent(Int_t,Int_t,Int_t,TObject*)");

    // Polls the pad's axis range and marker positions every 200ms
    // - CheckRangeChanged() redraws markers when zoom/pan changes
    // - CheckMarkersChanged() updates settings when markers are dragged
    static TTimer *pollTimer = new TTimer();
    pollTimer->Connect("Timeout()", 0, 0, "CheckRangeChanged()");
    pollTimer->Connect("Timeout()", 0, 0, "CheckMarkersChanged()");
    pollTimer->Start(200, kFALSE);

    window = ROOT::RWebWindow::Create();
    window->SetMaxQueueLength(100);  // Increase from default 10 to handle verbose output
    std::string fname = __FILE__;
    auto pos = fname.find("WebShapeIt.cxx");
    std::string dir = (pos != std::string::npos) ? fname.substr(0, pos) : std::string("./");
    window->SetDefaultPage("file:" + dir + "webshapeit.html");
    window->SetDataCallBack(ProcessData);
    window->SetGeometry(1200, 700);
    window->Show();

    std::cout << "\nShapeIt 2.0 prototype running.\n";
}
