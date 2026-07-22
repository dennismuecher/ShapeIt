// WebShapeIt.cxx
//
// BACKUP TIMESTAMP: 2026-07-13 - Before string comparison refactor
// TO RESTORE: Copy from git or search for "BACKUP_TIMESTAMP" in version control
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
#include "TRandom3.h"
#include "TMarker.h"
#include "TPaveText.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <climits>  // for PATH_MAX
#include <cstdlib>  // for realpath

#include "ShapeSetting.C"
#include "ShapeMatrix.C"
#include "ShapeGSF.C"
#include "ShapeRho.C"
#include "ShapeCollector.C"
#include "ShapeAlpha.C"
#include "ShapeController.C"

// ============================================================================
// String comparison helpers - eliminates manual length counting bugs
// ============================================================================

// Check if a string starts with a given prefix
inline bool starts_with(const std::string &str, const std::string &prefix) {
    return str.size() >= prefix.size() && 
           str.compare(0, prefix.size(), prefix) == 0;
}

// Extract payload after a prefix (returns empty string if prefix not found)
inline std::string after_prefix(const std::string &str, const std::string &prefix) {
    if (starts_with(str, prefix)) {
        return str.substr(prefix.size());
    }
    return "";
}

// Check prefix and extract payload in one step
inline bool extract_payload(const std::string &str, const std::string &prefix, 
                            std::string &payload) {
    if (starts_with(str, prefix)) {
        payload = str.substr(prefix.size());
        return true;
    }
    return false;
}

// ============================================================================

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
    std::cout << "[ParsePipeDoubles] Input: '" << s << "'" << std::endl;
    std::cout.flush();
    
    std::vector<double> out;
    std::stringstream ss(s);
    std::string tok;
    int idx = 0;
    while (std::getline(ss, tok, '|')) {
        std::cout << "[ParsePipeDoubles] Token " << idx << ": '" << tok << "'" << std::endl;
        std::cout.flush();
        
        double val = std::stod(tok);
        out.push_back(val);
        
        std::cout << "[ParsePipeDoubles] Converted to: " << val << std::endl;
        std::cout.flush();
        idx++;
    }
    
    std::cout << "[ParsePipeDoubles] Total parsed: " << out.size() << " values" << std::endl;
    std::cout.flush();
    
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
    
    // Combine paths
    std::string combined = baseDir + "/" + maybeRelative;
    
    // Use TSystem::GetWorkingDirectory() trick: temporarily prepend the absolute baseDir
    // and let ROOT normalize the path, which will properly resolve .. components
    TString normalized = gSystem->GetDirName(combined.c_str());
    normalized = gSystem->GetDirName(normalized.Data());
    normalized += "/";
    normalized += gSystem->BaseName(combined.c_str());
    
    // Actually, better: use realpath-style resolution
    // Convert to absolute path and normalize
    char resolved[PATH_MAX];
    if (realpath(combined.c_str(), resolved) != nullptr) {
        return std::string(resolved);
    }
    
    // If realpath fails (file doesn't exist yet), do manual normalization
    // This handles .. and . in the path
    std::vector<std::string> parts;
    std::istringstream ss(combined);
    std::string part;
    
    while (std::getline(ss, part, '/')) {
        if (part.empty() || part == ".") {
            continue; // skip empty parts and current directory
        } else if (part == "..") {
            if (!parts.empty() && parts.back() != "..") {
                parts.pop_back(); // go up one directory
            }
        } else {
            parts.push_back(part);
        }
    }
    
    // Reconstruct the path
    std::string result = combined[0] == '/' ? "/" : "";
    for (size_t i = 0; i < parts.size(); ++i) {
        if (i > 0) result += "/";
        result += parts[i];
    }
    
    return result;
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
              << "  peakPos: " << sett->peakPos[0] << " " << sett->peakPos[1] << "\n"
              << "  doubletPeakPos: " << sett->doubletPeakPos[0] << " " << sett->doubletPeakPos[1] << "\n"
              << "  tripletPeakPos: " << sett->tripletPeakPos[0] << " " << sett->tripletPeakPos[1] << "\n"
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
    
    // Clear canvas and set up fresh, exactly like SHOWMATRIX does
    canvas->cd();
    canvas->Clear();
    
    // Set pad margins to accommodate the manually-positioned color palette
    gPad->SetRightMargin(0.15);  // 15% right margin provides space for palette
    gPad->SetLeftMargin(0.12);   // 12% left margin for y-axis label
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.10);
    
    // Enable logarithmic z-axis scale
    gPad->SetLogz(1);
    
    TH2* hist = matrix->GetInputMatrix(BaseName(matrixPath));
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
    
    matrix->Diag();
    gPad->Modified();
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
    msg += std::to_string(sett->doDoublet[0] ? 1 : 0) + "|";  // Add doublet checkbox states
    msg += std::to_string(sett->doDoublet[1] ? 1 : 0) + "|";
    msg += std::to_string(sett->fixDoubletWidth[0] ? 1 : 0) + "|";  // Add doublet width fix toggles
    msg += std::to_string(sett->fixDoubletWidth[1] ? 1 : 0) + "|";
    msg += std::to_string(sett->doTriplet[0] ? 1 : 0) + "|";  // Add triplet checkbox states
    msg += std::to_string(sett->doTriplet[1] ? 1 : 0) + "|";
    msg += std::to_string(sett->fixTripletWidth[0] ? 1 : 0) + "|";  // Add triplet width fix toggles
    msg += std::to_string(sett->fixTripletWidth[1] ? 1 : 0) + "|";
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
    msg += std::to_string(sett->widthCal[1][0]) + "|" + std::to_string(sett->widthCal[1][1]) + "|";
    msg += std::to_string(sett->lit_alpha) + "|" + std::to_string(sett->lit_norm) + "|";
    msg += std::to_string(sett->fixPeakPos[0] ? 1 : 0) + "|" + std::to_string(sett->peakPos[0]) + "|";
    msg += std::to_string(sett->fixPeakPos[1] ? 1 : 0) + "|" + std::to_string(sett->peakPos[1]) + "|";
    msg += std::to_string(sett->fixDoubletPeakPos[0] ? 1 : 0) + "|" + std::to_string(sett->doubletPeakPos[0]) + "|";
    msg += std::to_string(sett->fixDoubletPeakPos[1] ? 1 : 0) + "|" + std::to_string(sett->doubletPeakPos[1]) + "|";
    msg += std::to_string(sett->fixTripletPeakPos[0] ? 1 : 0) + "|" + std::to_string(sett->tripletPeakPos[0]) + "|";
    msg += std::to_string(sett->fixTripletPeakPos[1] ? 1 : 0) + "|" + std::to_string(sett->tripletPeakPos[1]);
    

    
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

// Global markers for level 1 and level 2 peak positions
TMarker *gLevel1PeakMarker = nullptr;
TMarker *gLevel2PeakMarker = nullptr;
TMarker *gDoublet1PeakMarker = nullptr;
TMarker *gDoublet2PeakMarker = nullptr;
TMarker *gTriplet1PeakMarker = nullptr;
TMarker *gTriplet2PeakMarker = nullptr;

void DrawMarkers(bool usePadRange = false)
{
    for (int i = 0; i < 4; i++) {
        if (gMarkerLine[i]) canvas->GetListOfPrimitives()->Remove(gMarkerLine[i]);
        if (gDoubletLine[i]) canvas->GetListOfPrimitives()->Remove(gDoubletLine[i]);
        if (gBgBox[i]) canvas->GetListOfPrimitives()->Remove(gBgBox[i]);
    }
    
    // Remove the peak markers
    if (gLevel1PeakMarker) canvas->GetListOfPrimitives()->Remove(gLevel1PeakMarker);
    if (gLevel2PeakMarker) canvas->GetListOfPrimitives()->Remove(gLevel2PeakMarker);
    if (gDoublet1PeakMarker) canvas->GetListOfPrimitives()->Remove(gDoublet1PeakMarker);
    if (gDoublet2PeakMarker) canvas->GetListOfPrimitives()->Remove(gDoublet2PeakMarker);
    if (gTriplet1PeakMarker) canvas->GetListOfPrimitives()->Remove(gTriplet1PeakMarker);
    if (gTriplet2PeakMarker) canvas->GetListOfPrimitives()->Remove(gTriplet2PeakMarker);

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
        // Use 1.15 (15%) padding to ensure star markers at 110% amplitude are visible
        y2 = gCurrentHist ? gCurrentHist->GetMaximum() * 1.15 : 100;
    }
    
    // Helper function to get peak amplitude at a given X position
    auto getPeakAmplitude = [&](double xPos) -> double {
        if (!gCurrentHist) return (y1 + y2) / 2.0;  // fallback to mid-point
        int bin = gCurrentHist->FindBin(xPos);
        double amplitude = gCurrentHist->GetBinContent(bin);
        return amplitude > 0 ? amplitude : (y1 + y2) / 2.0;  // fallback if bin is empty
    };

    // Draw ONLY the main fit region markers (Level 1 and Level 2 boundaries)
    // NO MORE doublet/triplet level markers - only fit regions
    for (int i = 0; i < 4; i++) {
        gMarkerLine[i] = new TLine(sett->levEne[i], y1, sett->levEne[i], y2);
        if (i < 2) {
            gMarkerLine[i]->SetLineColor(kRed);  // Level 1 boundaries
        } else {
            gMarkerLine[i]->SetLineColor(kOrange);  // Level 2 boundaries
        }
        gMarkerLine[i]->SetLineWidth(2);
        if (sett->levEne[i] >= xmin && sett->levEne[i] <= xmax)
            gMarkerLine[i]->Draw();
    }
    
    // Draw star markers for ALL active peak positions
    // Star Y-position is at 110% of peak amplitude
    // Star X-position uses the input field values directly
    
    // Peak 1 Level 1 (always active - always show)
    // X-position: Use peakPos directly from settings file (sett->peakPos[0])
    double peak1Pos = sett->peakPos[0];
    double peak1Amp = getPeakAmplitude(peak1Pos);
    gLevel1PeakMarker = new TMarker(peak1Pos, 1.1 * peak1Amp, 29);  // 29 = star, Y at 110% amplitude
    gLevel1PeakMarker->SetMarkerColor(kRed);
    gLevel1PeakMarker->SetMarkerSize(2.0);
    if (peak1Pos >= xmin && peak1Pos <= xmax)
        gLevel1PeakMarker->Draw();
    
    // Peak 1 Level 2 (always active - always show)
    // X-position: Use peakPos directly from settings file (sett->peakPos[1])
    double peak2Pos = sett->peakPos[1];
    double peak2Amp = getPeakAmplitude(peak2Pos);
    gLevel2PeakMarker = new TMarker(peak2Pos, 1.1 * peak2Amp, 29);  // Y at 110% amplitude
    gLevel2PeakMarker->SetMarkerColor(kOrange);
    gLevel2PeakMarker->SetMarkerSize(2.0);
    if (peak2Pos >= xmin && peak2Pos <= xmax)
        gLevel2PeakMarker->Draw();
    
    // Peak 2 Level 1 (doublet - show star if active)
    // X-position: Use doubletPeakPos directly from settings file (sett->doubletPeakPos[0])
    if (sett->doDoublet[0]) {
        double doublet1Pos = sett->doubletPeakPos[0];
        double doublet1Amp = getPeakAmplitude(doublet1Pos);
        gDoublet1PeakMarker = new TMarker(doublet1Pos, 1.1 * doublet1Amp, 29);  // Y at 110% amplitude
        gDoublet1PeakMarker->SetMarkerColor(kRed);
        gDoublet1PeakMarker->SetMarkerSize(2.0);
        if (doublet1Pos >= xmin && doublet1Pos <= xmax)
            gDoublet1PeakMarker->Draw();
    }
    
    // Peak 2 Level 2 (doublet - show star if active)
    // X-position: Use doubletPeakPos directly from settings file (sett->doubletPeakPos[1])
    if (sett->doDoublet[1]) {
        double doublet2Pos = sett->doubletPeakPos[1];
        double doublet2Amp = getPeakAmplitude(doublet2Pos);
        gDoublet2PeakMarker = new TMarker(doublet2Pos, 1.1 * doublet2Amp, 29);  // Y at 110% amplitude
        gDoublet2PeakMarker->SetMarkerColor(kOrange);
        gDoublet2PeakMarker->SetMarkerSize(2.0);
        if (doublet2Pos >= xmin && doublet2Pos <= xmax)
            gDoublet2PeakMarker->Draw();
    }
    
    // Peak 3 Level 1 (triplet - show star if active)
    // X-position: Use tripletPeakPos directly from settings file (sett->tripletPeakPos[0])
    if (sett->doTriplet[0]) {
        double triplet1Pos = sett->tripletPeakPos[0];
        double triplet1Amp = getPeakAmplitude(triplet1Pos);
        gTriplet1PeakMarker = new TMarker(triplet1Pos, 1.1 * triplet1Amp, 29);  // Y at 110% amplitude
        gTriplet1PeakMarker->SetMarkerColor(kRed);
        gTriplet1PeakMarker->SetMarkerSize(2.0);
        if (triplet1Pos >= xmin && triplet1Pos <= xmax)
            gTriplet1PeakMarker->Draw();
    }
    
    // Peak 3 Level 2 (triplet - show star if active)
    // X-position: Use tripletPeakPos directly from settings file (sett->tripletPeakPos[1])
    if (sett->doTriplet[1]) {
        double triplet2Pos = sett->tripletPeakPos[1];
        double triplet2Amp = getPeakAmplitude(triplet2Pos);
        gTriplet2PeakMarker = new TMarker(triplet2Pos, 1.1 * triplet2Amp, 29);  // Y at 110% amplitude
        gTriplet2PeakMarker->SetMarkerColor(kOrange);
        gTriplet2PeakMarker->SetMarkerSize(2.0);
        if (triplet2Pos >= xmin && triplet2Pos <= xmax)
            gTriplet2PeakMarker->Draw();
    }

    // Background boxes
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
int gSkipRangeChecks = 0;  // Skip this many polling cycles (for context menu operations)

// Track last marker positions to detect when they're dragged
double gLastLevEne[4] = {0, 0, 0, 0};
double gLastBgEne[2][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
bool gHaveLastMarkers = false;

// Track width calibration plot axis ranges (set once when plot is created)
double gWidthCalibXMin = 0, gWidthCalibXMax = 0, gWidthCalibYMin = 0, gWidthCalibYMax = 0;
bool gHaveWidthCalibRanges = false;

// Monte Carlo state for timer-based iteration
struct MCState {
    int currentIter;
    int totalIters;
    double exiMin, exiMax, savedExiLow, savedExiHigh;
    TGraph *graph;
    std::vector<double> alphas, chi2s, exis;
    TRandom3 rng;
    unsigned connid;
    bool stopRequested;  // Flag to indicate user requested early stop
};

MCState *gMCState = nullptr;

// Forward declarations
void RunShapeIt(unsigned connid);
void RunWidthCalibration(unsigned connid);
void RunMonteCarlo(unsigned connid, int nIterations, double exiLowMin, double exiLowMax);
void MCTimerCallback();

// Timer callback for Monte Carlo iterations - runs ONE iteration at a time with live display
void MCTimerCallback()
{
    // Check if user requested stop OR all iterations complete
    if (!gMCState || gMCState->currentIter >= gMCState->totalIters || gMCState->stopRequested) {
        if (gMCState) {
            // Check if we have any results at all
            if (gMCState->alphas.empty()) {
                std::cout << "Monte Carlo stopped with no results." << std::endl;
                window->Send(gMCState->connid, "MC_STOPPED:No iterations completed before stop.");
                
                // Restore settings
                sett->exiEne[0] = gMCState->savedExiLow;
                sett->exiEne[1] = gMCState->savedExiHigh;
                sett->doMC = false;
                
                delete gMCState->graph;
                delete gMCState;
                gMCState = nullptr;
                return;
            }
            
            std::string stopReason = gMCState->stopRequested ? "stopped by user" : "completed all iterations";
            std::cout << "Monte Carlo " << stopReason << " after " << gMCState->currentIter << " iterations, creating final display..." << std::endl;
            
            // Final display with best fit marker
            canvas->cd();
            canvas->Clear();
            canvas->Divide(2, 2);
            
            auto minIt = std::min_element(gMCState->chi2s.begin(), gMCState->chi2s.end());
            int minIdx = std::distance(gMCState->chi2s.begin(), minIt);
            double bestAlpha = gMCState->alphas[minIdx];
            double bestChi2 = *minIt;
            double bestExLow = gMCState->exis[minIdx];
            
            std::cout << "Best fit: alpha=" << bestAlpha << ", chi2=" << bestChi2 
                      << ", Ex_low=" << bestExLow << std::endl;
            
            // Top left: chi2 vs delta alpha with best point marked
            canvas->cd(1);
            TGraph *gFinal = new TGraph(gMCState->alphas.size(), &gMCState->alphas[0], &gMCState->chi2s[0]);
            gFinal->SetTitle("Chi^{2} vs #Delta#alpha (Final);#Delta#alpha [MeV^{-1}];#chi^{2}");
            gFinal->SetMarkerStyle(20);
            gFinal->SetMarkerSize(0.8);
            gFinal->SetMarkerColor(kBlue);
            gFinal->Draw("AP");
            
            TMarker *bestPoint = new TMarker(bestAlpha, bestChi2, 29);
            bestPoint->SetMarkerColor(kRed);
            bestPoint->SetMarkerSize(2);
            bestPoint->Draw();
            
            // Top right: Ex_low vs delta alpha with best point marked
            canvas->cd(2);
            TGraph *gExFinal = new TGraph(gMCState->exis.size(), &gMCState->exis[0], &gMCState->alphas[0]);
            gExFinal->SetTitle("Ex_{low} vs #Delta#alpha;Ex_{low} [keV];#Delta#alpha [MeV^{-1}]");
            gExFinal->SetMarkerStyle(20);
            gExFinal->SetMarkerSize(0.8);
            gExFinal->SetMarkerColor(kRed);
            gExFinal->Draw("AP");
            
            TMarker *bestPointEx = new TMarker(bestExLow, bestAlpha, 29);
            bestPointEx->SetMarkerColor(kBlue);
            bestPointEx->SetMarkerSize(2);
            bestPointEx->Draw();
            
            // Bottom left: delta alpha histogram with best line
            canvas->cd(3);
            double aMin = *std::min_element(gMCState->alphas.begin(), gMCState->alphas.end());
            double aMax = *std::max_element(gMCState->alphas.begin(), gMCState->alphas.end());
            TH1D *hFinal = new TH1D("hAlphaFinal", 
                                     Form("#Delta#alpha Distribution (Best: %.3f);#Delta#alpha [MeV^{-1}];Counts", bestAlpha),
                                     50, aMin, aMax);
            for (double a : gMCState->alphas) hFinal->Fill(a);
            hFinal->SetLineColor(kBlue);
            hFinal->SetLineWidth(2);
            hFinal->Draw();
            
            TLine *bestLine = new TLine(bestAlpha, 0, bestAlpha, hFinal->GetMaximum());
            bestLine->SetLineColor(kRed);
            bestLine->SetLineWidth(2);
            bestLine->Draw();
            
            // Bottom right: Ex_low histogram with best line
            canvas->cd(4);
            TH1D *hExFinal = new TH1D("hExFinal",
                                       Form("Ex_{low} Distribution (Best: %.0f);Ex_{low} [keV];Counts", bestExLow),
                                       50, gMCState->exiMin, gMCState->exiMax);
            for (double ex : gMCState->exis) hExFinal->Fill(ex);
            hExFinal->SetLineColor(kRed);
            hExFinal->SetLineWidth(2);
            hExFinal->Draw();
            
            TLine *bestLineEx = new TLine(bestExLow, 0, bestExLow, hExFinal->GetMaximum());
            bestLineEx->SetLineColor(kBlue);
            bestLineEx->SetLineWidth(2);
            bestLineEx->Draw();
            
            canvas->cd();
            gPad->Modified();
            PushCanvasUpdate();
            
            // Restore settings
            sett->exiEne[0] = gMCState->savedExiLow;
            sett->exiEne[1] = gMCState->savedExiHigh;
            sett->doMC = false;  // Disable MC mode after completion
            
            std::string msg = "MC_RESULT:" + std::to_string(bestAlpha) + "|" + std::to_string(bestChi2) 
                              + "|" + std::to_string(bestExLow) + "|" 
                              + std::to_string(gMCState->currentIter) + "|"
                              + std::to_string(gMCState->totalIters);
            window->Send(gMCState->connid, msg);
            
            if (gMCState->stopRequested) {
                std::cout << "Monte Carlo stopped by user after " << gMCState->currentIter << " of " 
                          << gMCState->totalIters << " iterations. Best alpha: " << bestAlpha 
                          << " (chi2: " << bestChi2 << ", Ex_low: " << bestExLow << " keV)" << std::endl;
            } else {
                std::cout << "Monte Carlo complete. Best alpha: " << bestAlpha 
                          << " (chi2: " << bestChi2 << ", Ex_low: " << bestExLow << " keV)" << std::endl;
            }
            std::cout << "MC mode disabled (doMC = false)" << std::endl;
            
            delete gMCState->graph;
            delete gMCState;
            gMCState = nullptr;
        }
        return;
    }
    
    // Run ONE iteration
    // Randomize lower excitation boundary
    double exiLow = gMCState->rng.Uniform(gMCState->exiMin, gMCState->exiMax);
    sett->exiEne[0] = exiLow;
    sett->exiEne[1] = gMCState->savedExiHigh;
    
    // Run analysis (suppress verbose output)
    int savedVerbose = sett->verbose;
    sett->verbose = 0;
    ShapeCollector *tempColl = ShapeController::RunAnalysis(sett, matrix);
    
    if (!tempColl) {
        std::cout << "  WARNING: iteration " << (gMCState->currentIter+1) << " returned null" << std::endl;
        sett->verbose = savedVerbose;
        gMCState->currentIter++;
        return;
    }
    
    // Fit alpha to find best value for this iteration
    ShapeAlpha *alphaFit = ShapeController::FitAlpha(sett, tempColl);
    double alpha = alphaFit->getMinAlpha();
    double chi2 = alphaFit->getMinChi2();
    
    delete alphaFit;
    sett->verbose = savedVerbose;
    
    gMCState->alphas.push_back(alpha);
    gMCState->chi2s.push_back(chi2);
    gMCState->exis.push_back(exiLow);
    
    std::cout << "Iteration " << (gMCState->currentIter+1) << "/" << gMCState->totalIters 
              << ": Ex_low=" << exiLow << ", alpha=" << alpha << ", chi2=" << chi2 << std::endl;
    
    // Apply alpha transformation to get properly scaled gSF data
    tempColl->Transform(sett->lit_norm, alpha);
    
    // Get THIS iteration's data AFTER alpha transformation (don't accumulate)
    TGraphErrors *mergedGraph = tempColl->getMergedGraph();
    TGraph *litGraph = tempColl->getLitGraph();
    TGraphErrors *litCopy = nullptr;
    
    if (litGraph && litGraph->GetN() > 0) {
        TGraphErrors *litOrig = dynamic_cast<TGraphErrors*>(litGraph);
        if (litOrig) litCopy = new TGraphErrors(*litOrig);
    }
    
    // CREATE FRESH GRAPH FOR THIS ITERATION ONLY (don't accumulate into gMCState->graph)
    TGraph *freshGraph = nullptr;
    if (mergedGraph && mergedGraph->GetN() > 0) {
        std::vector<double> x(mergedGraph->GetN()), y(mergedGraph->GetN());
        for (int j = 0; j < mergedGraph->GetN(); j++) {
            x[j] = mergedGraph->GetX()[j];
            y[j] = mergedGraph->GetY()[j];
        }
        freshGraph = new TGraph((int)x.size(), x.data(), y.data());
    }
    
    delete tempColl;
    
    // Redraw: lit first, then THIS iteration's points on top
    canvas->cd();
    gPad->SetLogy(1);  // Enable log scale on Y axis
    
    if (litCopy) {
        litCopy->SetLineColor(kRed);
        litCopy->SetLineWidth(2);
        litCopy->SetFillColor(kRed-10);
        litCopy->SetFillStyle(3013);
        litCopy->Draw("AL3");
        litCopy->SetTitle(Form("MC Iteration %d/%d (#Delta#alpha=%.3f);E_{#gamma} (keV);f(E_{#gamma}) (MeV^{-3})", 
                               gMCState->currentIter+1, gMCState->totalIters, alpha));
    }
    
    if (freshGraph) {
        freshGraph->SetMarkerStyle(20);
        freshGraph->SetMarkerSize(0.5);
        freshGraph->SetMarkerColor(kBlue);
        freshGraph->SetLineColor(kBlue);
        
        if (litCopy) {
            freshGraph->Draw("P SAME");
        } else {
            freshGraph->Draw("AP");
            freshGraph->SetTitle(Form("MC Iteration %d/%d (#Delta#alpha=%.3f);E_{#gamma} (keV);f(E_{#gamma}) (MeV^{-3})", 
                                       gMCState->currentIter+1, gMCState->totalIters, alpha));
        }
    }
    
    gMCState->currentIter++;
    
    gPad->Modified();
    PushCanvasUpdate();
}

// Monte Carlo to find best alpha: runs multiple iterations with varying lower
// excitation energy boundary, accumulates chi2 and alpha values, displays
// live progress in main canvas. Uses timer-based batching for live updates.
void RunMonteCarlo(unsigned connid, int nIterations, double exiLowMin, double exiLowMax)
{
    std::cout << "\n=== RunMonteCarlo ENTERED ===" << std::endl;
    std::cout << "Parameters: nIterations=" << nIterations 
              << ", exiLowMin=" << exiLowMin << ", exiLowMax=" << exiLowMax << std::endl;
    
    if (!matrix || !sett->doOslo) {
        std::cout << "ERROR: Prerequisites not met!" << std::endl;
        std::cout << "  matrix = " << (matrix ? "loaded" : "NULL") << std::endl;
        std::cout << "  doOslo = " << sett->doOslo << std::endl;
        window->Send(connid, "ERROR: Monte Carlo requires matrix and Oslo literature data loaded.");
        return;
    }
    
    // Enable MC mode for Gaussian randomization of literature and experimental data
    sett->doMC = true;
    std::cout << "MC mode enabled (doMC = true) - data will be randomized" << std::endl;
    
    std::cout << "Prerequisites OK, starting MC run..." << std::endl;
    
    // Clear canvas and set up progress display
    canvas->cd();
    canvas->Clear();
    gDisplayMode = 8;
    for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
    
    // Initialize MC state
    gMCState = new MCState{
        0,                    // currentIter
        nIterations,          // totalIters
        exiLowMin,            // exiMin
        exiLowMax,            // exiMax
        sett->exiEne[0],      // savedExiLow
        sett->exiEne[1],      // savedExiHigh
        new TGraph(),         // graph
        {},                   // alphas
        {},                   // chi2s
        {},                   // exis
        TRandom3(0),          // rng
        connid,               // connid
        false                 // stopRequested
    };
    
    gMCState->graph->SetTitle("MC Progress: Accumulating gSF Results;E_{#gamma} (keV);f(E_{#gamma}) (MeV^{-3})");
    gMCState->graph->SetMarkerStyle(20);
    gMCState->graph->SetMarkerSize(0.5);
    gMCState->graph->SetMarkerColor(kBlue);
    gMCState->graph->Draw("AP");
    PushCanvasUpdate();
    
    // Tell browser that MC has started so button can change to "Stop"
    window->Send(connid, "MC_STARTED");
    
    // Start timer to run iterations one at a time with live display
    static TTimer *mcTimer = new TTimer();
    mcTimer->Connect("Timeout()", 0, 0, "MCTimerCallback()");
    mcTimer->Start(100, kFALSE); // Every 100ms, run 1 iteration
    
    std::cout << "MC timer started - will run " << nIterations << " iterations with live display" << std::endl;
}

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

    // If we're skipping checks (context menu operation in progress), decrement and return
    if (gSkipRangeChecks > 0) {
        gSkipRangeChecks--;
        return;
    }

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
    double bgEne[2][4];
    
    for (int i = 0; i < 4; i++)
        levEne[i] = gMarkerLine[i]->GetX1();
    
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
    }
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 4; j++)
            gLastBgEne[i][j] = bgEne[i][j];
    gHaveLastMarkers = true;
    
    // If changed, update settings and UI
    if (changed) {
        for (int i = 0; i < 4; i++) {
            sett->levEne[i] = levEne[i];
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
        msg += std::to_string(sett->bgEne[0][0]) + "|" + std::to_string(sett->bgEne[0][1]) + "|";
        msg += std::to_string(sett->bgEne[0][2]) + "|" + std::to_string(sett->bgEne[0][3]) + "|";
        msg += std::to_string(sett->bgEne[1][0]) + "|" + std::to_string(sett->bgEne[1][1]) + "|";
        msg += std::to_string(sett->bgEne[1][2]) + "|" + std::to_string(sett->bgEne[1][3]);
        window->Send(0, msg);
        
        // If in Autofit mode (mode 2), re-fit and redraw the current projection
        // This updates the Gaussian fits on the histogram without recalculating all gSF
        if (sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
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

    // On right-click (context menu), pause polling for 3 cycles (~1.5 seconds)
    // This gives ROOT's context menu operations (like unzoom) time to complete
    // without our marker redrawing interfering
    if (event == kButton3Down) {
        gSkipRangeChecks = 3;  // Skip next 3 polling cycles
        return;
    }

    CheckRangeChanged();

    if (event == kButton1Up || event == kButton2Up || event == kButton3Up) {
        for (int i = 0; i < 4; i++)
            sett->levEne[i] = gMarkerLine[i]->GetX1();

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
    
    // Add chi2 and alpha text box if Oslo data is displayed (matches ShapeFrame::getPaveTextgSF())
    if (sett->doOslo) {
        TPaveText *t = new TPaveText(0.8, 0.85, 0.95, 0.95, "brNDC");
        t->SetTextSize(0.025);
        t->SetTextAlign(13);
        t->SetFillColor(10);
        t->SetTextColor(61);
        t->AddText(Form("slope #alpha: %4.2f", sett->lit_alpha));
        t->AddText(Form("#chi^{2} value: %4.2f", gSFColl->getChi2()));
        t->Draw();
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
    // TEMPORARILY DISABLED: stdout buffering to help debug Monte Carlo issue
    // Just print directly to terminal like a normal C++ program
    
    // Suppress noisy debug messages for channel setup and width calibration updates
    // (Also suppress common messages like SHOWPROJ, EXIT to reduce console noise)
    if (!starts_with(arg, "channel:") && 
        !starts_with(arg, "UPDATE_WIDTH_CALIB_LINES:") &&
        !starts_with(arg, "SHOWPROJ") &&
        !starts_with(arg, "EXIT")) {
        std::cout << "Got message from browser: " << arg << std::endl;
    }

    if (starts_with(arg, "channel:")) {
        int chid = std::stoi(after_prefix(arg, "channel:"));
        auto web_imp = dynamic_cast<TWebCanvas *>(canvas->GetCanvasImp());
        if (web_imp) {
            web_imp->ShowWebWindow({ window, connid, chid });
            // DON'T ForceUpdate yet - wait until after we've drawn the matrix with log scale
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
                // Send efficiency correction status
                std::string effiMsg = "EFFI_INFO:" + std::string(sett->doEffi ? "1" : "0") + "|" + sett->effiFileName;
                window->Send(connid, effiMsg);
            }
        } else {
            // No matrix loaded - push an empty canvas
            if (web_imp) web_imp->ForceUpdate();
            // Still send efficiency correction status
            std::string effiMsg = "EFFI_INFO:" + std::string(sett->doEffi ? "1" : "0") + "|" + sett->effiFileName;
            window->Send(connid, effiMsg);
        }
        
        // Sync all settings to UI (this updates all form fields to match loaded settings)
        SendSettingsSync(connid);
        
        // Send status message AFTER everything else is synced (same as LOAD_SETTINGS does)
        // Only send if we actually loaded settings from a file
        if (!sett->settFileName.empty()) {
            window->Send(connid, "Settings loaded: " + sett->settFileName);
        }
    }
    else if (starts_with(arg, "OPEN:")) {
        std::string path = after_prefix(arg, "OPEN:");
        
        // Use ROOT's ExpandPathName to resolve relative paths and canonicalize
        // ExpandPathName returns a new char* with the expanded path
        char* expandedPath = gSystem->ExpandPathName(path.c_str());
        if (expandedPath) {
            path = expandedPath;
            delete[] expandedPath;
        }

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
    else if (starts_with(arg, "SELECTMATRIX:")) {
        if (!matrix) {
            window->Send(connid, "No matrix file open yet.");
            return;
        }
        int idx = std::stoi(after_prefix(arg, "SELECTMATRIX:"));
        SendMatrixListAndSelect(connid, currentMatrixPath, idx);
        SendNBins(connid);
        // Enable width calibration since we have a matrix
        window->Send(connid, "WIDTH_CALIB_AVAILABLE:1");
    }
    else if (starts_with(arg, "OSLO:")) {
        std::string path = after_prefix(arg, "OSLO:");
        
        // Use ROOT's ExpandPathName to resolve relative paths and canonicalize
        char* expandedPath = gSystem->ExpandPathName(path.c_str());
        if (expandedPath) {
            path = expandedPath;
            delete[] expandedPath;
        }

        if (gSystem->AccessPathName(path.c_str())) {
            window->Send(connid, "Literature file not found: " + path);
            return;
        }

        sett->osloFileName = path;
        sett->doOslo = true;
        window->Send(connid, "gSF literature file set: " + path);
    }
    else if (starts_with(arg, "RHO:")) {
        std::string path = after_prefix(arg, "RHO:");
        
        // Use ROOT's ExpandPathName to resolve relative paths and canonicalize
        // This handles '..' segments and converts to absolute paths
        char* expandedPath = gSystem->ExpandPathName(path.c_str());
        if (expandedPath) {
            path = expandedPath;
            delete[] expandedPath;
        }

        if (gSystem->AccessPathName(path.c_str())) {
            window->Send(connid, "NLD file not found: " + path);
            return;
        }

        sett->rhoFileName = path;
        
        // Create ShapeRho object to read the file and check for unit conversion
        ShapeRho* tempRho = new ShapeRho(sett);
        
        // Check if MeV to keV conversion happened and notify browser
        if (tempRho && tempRho->wasConvertedFromMeV) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << tempRho->originalMaxEnergy;
            std::string warning = "UNIT_WARNING:" + oss.str();
            window->Send(connid, warning);
        }
        
        delete tempRho;  // Clean up temporary object
        
        window->Send(connid, "NLD literature file set: " + path);
    }
    else if (starts_with(arg, "DISCRETE:")) {
        std::string path = after_prefix(arg, "DISCRETE:");
        
        // Use ROOT's ExpandPathName to resolve relative paths and canonicalize
        char* expandedPath = gSystem->ExpandPathName(path.c_str());
        if (expandedPath) {
            path = expandedPath;
            delete[] expandedPath;
        }

        if (gSystem->AccessPathName(path.c_str())) {
            window->Send(connid, "Discrete levels file not found: " + path);
            return;
        }

        sett->discreteLevelFile = path;
        window->Send(connid, "Discrete levels file set: " + path);
    }
    else if (starts_with(arg, "EFFI:")) {
        std::string path = after_prefix(arg, "EFFI:");
        
        // Use ROOT's ExpandPathName to resolve relative paths and canonicalize
        char* expandedPath = gSystem->ExpandPathName(path.c_str());
        if (expandedPath) {
            path = expandedPath;
            delete[] expandedPath;
        }

        if (gSystem->AccessPathName(path.c_str())) {
            window->Send(connid, "Efficiency correction file not found: " + path);
            return;
        }

        sett->effiFileName = path;
        sett->doEffi = true;
        sett->readEffi();  // Read the file immediately
        
        // Send updated status back to frontend
        std::string msg = "EFFI_INFO:1|" + sett->effiFileName;
        window->Send(connid, msg);
        window->Send(connid, "Efficiency correction file loaded: " + path);
    }
    else if (starts_with(arg, "DOEFFI:")) {
        int val = std::stoi(after_prefix(arg, "DOEFFI:"));
        sett->doEffi = (val == 1);
        
        // Send updated status back to frontend
        std::string msg = "EFFI_INFO:" + std::string(sett->doEffi ? "1" : "0") + "|" + sett->effiFileName;
        window->Send(connid, msg);
        window->Send(connid, sett->doEffi ? "Energy-dependent efficiency correction enabled" : "Constant efficiency correction enabled");
    }
    else if (starts_with(arg, "LISTDIR:")) {
        SendDirListing(connid, after_prefix(arg, "LISTDIR:"));
    }
    else if (starts_with(arg, "LOAD_SETTINGS:")) {
        std::string path = after_prefix(arg, "LOAD_SETTINGS:");

        if (gSystem->AccessPathName(path.c_str())) {
            window->Send(connid, "Settings file not found: " + path);
            return;
        }

        sett->settFileName = path;
        sett->ReadSettings();

        // dataFileName/osloFileName/rhoFileName/discreteLevelFile as read are literal raw lines from the file,
        // unresolved -- if relative, resolve against the settings file's own
        // directory rather than wherever ROOT happened to be launched from.
        std::string settDir = DirName(path);
        sett->dataFileName = ResolveRelativeTo(settDir, sett->dataFileName);
        sett->osloFileName = ResolveRelativeTo(settDir, sett->osloFileName);
        sett->rhoFileName = ResolveRelativeTo(settDir, sett->rhoFileName);
        sett->discreteLevelFile = ResolveRelativeTo(settDir, sett->discreteLevelFile);
        sett->effiFileName = ResolveRelativeTo(settDir, sett->effiFileName);

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
                // Send efficiency correction status
                std::string effiMsg = "EFFI_INFO:" + std::string(sett->doEffi ? "1" : "0") + "|" + sett->effiFileName;
                window->Send(connid, effiMsg);
            }
            else
                window->Send(connid, "Warning: matrix '" + sett->matrixName + "' from settings not found in " + sett->dataFileName);
        }

        SendSettingsSync(connid);
        window->Send(connid, "Settings loaded: " + path);
    }
    else if (starts_with(arg, "SAVE_SETTINGS:")) {
        // Uses sett as it currently stands -- i.e. whatever the last "ShapeIt!"
        // run (or Options change) set it to. Click ShapeIt! at least once before
        // saving so the Peaks panel's current values are actually captured.
        std::string path = after_prefix(arg, "SAVE_SETTINGS:");
        sett->settFileName = path;
        sett->SaveSettings();
        window->Send(connid, "Settings saved: " + path);
    }
    else if (starts_with(arg, "OPTIONS:")) {
        // order: doInterpol|doOslo|doSlidingWindow|doBackground|doWidthCal
        auto v = ParsePipeDoubles(after_prefix(arg, "OPTIONS:"));
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
    else if (starts_with(arg, "DISPLAY_OPTIONS:")) {
        // order: displaySingle|displayAvg|colour
        auto v = ParsePipeDoubles(after_prefix(arg, "DISPLAY_OPTIONS:"));
        if (v.size() != 3) {
            window->Send(connid, "Malformed DISPLAY_OPTIONS message.");
            return;
        }
        sett->displaySingle = v[0] != 0.0;
        sett->displayAvg    = v[1] != 0.0;
        sett->colour        = v[2] != 0.0;
    }
    else if (starts_with(arg, "BINSIZE:")) {
        // order: lo|hi|isVariation
        auto v = ParsePipeDoubles(after_prefix(arg, "BINSIZE:"));
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
    else if (starts_with(arg, "NBINSLO:")) {
        // mirrors DoNumberEntry's id==8 case: editing "Nr. of bins" (low)
        // recomputes the corresponding bin size via BinToSize(). Trusts the
        // typed count as-is -- does NOT recompute it back from the derived
        // size, which is what caused the earlier drift bug.
        int nLo = std::stoi(after_prefix(arg, "NBINSLO:"));
        if (nLo < 1) nLo = 1;
        sett->nOfBins = nLo;
        sett->exi_size[0] = sett->BinToSize();
        if (!sett->doBinVariation)
            sett->exi_size[1] = sett->exi_size[0];
        int nHigh = sett->doBinVariation ? sett->SizeToBin(sett->exi_size[1]) : sett->nOfBins;
        SendBinSyncValues(connid, sett->nOfBins, nHigh);
    }
    else if (starts_with(arg, "NBINSHI:")) {
        // mirrors DoNumberEntry's id==12 case exactly, including its two
        // correction passes: hi bin size must exceed lo by at least 50 keV,
        // and the resulting bin count must be at least 3. sett->nOfBins (the
        // low count) is untouched here -- echoed back as-is, not recomputed.
        int nHigh = std::stoi(after_prefix(arg, "NBINSHI:"));
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
    else if (starts_with(arg, "INTPARAMS:")) {
        // order: minCounts|scaling|autoScale|effCorr
        auto v = ParsePipeDoubles(after_prefix(arg, "INTPARAMS:"));
        if (v.size() != 4) {
            window->Send(connid, "Malformed INTPARAMS message.");
            return;
        }
        sett->minCounts   = (int)v[0];
        sett->gSF_norm    = v[1];
        sett->doAutoScale = v[2] != 0.0;
        sett->eff_corr    = v[3];
    }
    else if (starts_with(arg, "VERBOSE:")) {
        sett->verbose = std::stoi(after_prefix(arg, "VERBOSE:"));
    }
    else if (starts_with(arg, "MODE:")) {
        sett->mode = std::stoi(after_prefix(arg, "MODE:")); // 1 = Integration, 2 = Autofit
    }
    else if (starts_with(arg, "BGENERGIES:")) {
        // order: bgEne[0][0..3] | bgEne[1][0..3]  (8 values)
        auto v = ParsePipeDoubles(after_prefix(arg, "BGENERGIES:"));
        if (v.size() != 8) {
            window->Send(connid, "Malformed BGENERGIES message.");
            return;
        }
        for (int i = 0; i < 4; i++) sett->bgEne[0][i] = v[i];
        for (int i = 0; i < 4; i++) sett->bgEne[1][i] = v[4 + i];
        
        // If in Autofit mode and viewing a bin projection, refit with new background regions
        if (sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
            std::cout << "Background energies changed in Autofit mode: re-fitting bin " << gCurrentBin << "..." << std::endl;
            
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
            DrawMarkers(); // refreshes the draggable boxes to match, if a projection is shown
        }
    }
    else if (arg == "UPDATE_MARKERS") {
        // Redraw markers with current settings (triggered by energy changes in UI)
        DrawMarkers(true);
    }
    else if (arg == "SHOW_LEVELS_PANEL") {
        // When Levels panel is clicked, don't change the current view
        // The panel is just for editing parameters, not for forcing a specific display
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet -- open one first.");
            return;
        }
        
        // Don't switch views - just acknowledge the panel was opened
        std::cout << "Levels panel opened - keeping current view (mode " << gDisplayMode << ")" << std::endl;
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
        
        // Enable logarithmic z-axis scale
        gPad->SetLogz(1);
        
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
        
        gPad->Modified();
        PushCanvasUpdate();
    }
    else if (arg == "SHOWPROJ") {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet -- open one first.");
            return;
        }
        gCurrentBin = 1;
        gDisplayMode = 5;
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
        gCurrentHist = matrix->GetDiagEx(gCurrentBin, BaseName(currentMatrixPath));
        if (!gCurrentHist) {
            window->Send(connid, "ERROR: Failed to get projection histogram");
            return;
        }
        gCurrentHist->Draw();
        gHaveLastRange = false;
        CleanupAutofitDisplay();
        DrawMarkers();
        PushCanvasUpdate();
        window->Send(connid, "Showing bin 1 projection");
    }
    else if (arg == "SHOW_WIDTH_CALIB") {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet.");
            return;
        }
        
        // Check if we have existing calibration parameters from settings file
        bool hasExistingCalib = (sett->widthCal[0][0] != 0.0 || sett->widthCal[0][1] != 0.0 ||
                                  sett->widthCal[1][0] != 0.0 || sett->widthCal[1][1] != 0.0);
        
        // Run calibration to generate graph data points (needed for display)
        // but DON'T fit the data yet - just show the points
        window->Send(connid, "Generating width calibration data...");
        RunWidthCalibration(connid);
        
        // Extract the width data graphs (populated by RunWidthCalibration above)
        TGraph *T1 = matrix->getFitWidthGraph(0);
        TGraph *T2 = matrix->getFitWidthGraph(1);
        
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
        
        canvas->cd();
        canvas->Clear();
        // Null out marker pointers since Clear() deleted them
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
        gDisplayMode = 7;  // Width calibration display mode
        
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
        
        // Only draw fit lines if we have existing calibration parameters from settings file
        if (hasExistingCalib) {
            TF1 *fit1 = nullptr;
            TF1 *fit2 = nullptr;
            
            if (T1->GetN() > 0 && (sett->widthCal[0][0] != 0.0 || sett->widthCal[0][1] != 0.0)) {
                fit1 = new TF1("fit1", "[0] + [1]*x", xMin, xMax);
                fit1->SetParameter(0, sett->widthCal[0][0]);
                fit1->SetParameter(1, sett->widthCal[0][1]);
                fit1->SetLineColor(kRed);
                fit1->SetLineWidth(2);
                fit1->Draw("SAME");
            }
            
            if (T2->GetN() > 0 && (sett->widthCal[1][0] != 0.0 || sett->widthCal[1][1] != 0.0)) {
                fit2 = new TF1("fit2", "[0] + [1]*x", xMin, xMax);
                fit2->SetParameter(0, sett->widthCal[1][0]);
                fit2->SetParameter(1, sett->widthCal[1][1]);
                fit2->SetLineColor(kBlue);
                fit2->SetLineWidth(2);
                fit2->Draw("SAME");
            }
        }
        
        // Add legend (smaller size)
        TLegend *leg = new TLegend(0.75, 0.80, 0.90, 0.90);
        leg->SetFillColor(0);
        if (T1->GetN() > 0) leg->AddEntry(T1, "level 1", "lp");
        if (T2->GetN() > 0) leg->AddEntry(T2, "level 2", "lp");
        leg->Draw();
        
        // Enable the width calibration checkbox by telling frontend it's available
        window->Send(connid, "WIDTH_CALIB_AVAILABLE:1");
        
        // Send the current calibration parameters to UI (from settings file)
        std::string msg = "WIDTH_CALIB_PARAMS:";
        msg += std::to_string(sett->widthCal[0][0]) + "|" + std::to_string(sett->widthCal[0][1]) + "|";
        msg += std::to_string(sett->widthCal[1][0]) + "|" + std::to_string(sett->widthCal[1][1]);
        window->Send(connid, msg);
        
        // DON'T delete T1/T2 here - they're now drawn on the canvas and will be
        // cleaned up automatically when the canvas is cleared.
        
        PushCanvasUpdate();
    }
    else if (arg == "REFRESH_WIDTH_CALIB_FITS") {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet.");
            return;
        }
        
        if (gDisplayMode != 7) {
            window->Send(connid, "Not currently viewing width calibration.");
            return;
        }
        
        window->Send(connid, "Refreshing width calibration fits...");
        
        // Extract the width data graphs (already populated by previous SHOW_WIDTH_CALIB)
        TGraph *T1 = matrix->getFitWidthGraph(0);
        TGraph *T2 = matrix->getFitWidthGraph(1);
        
        // Perform linear fits on the width data to get NEW parameters
        if (T1->GetN() > 0) {
            TF1 *tempFit1 = new TF1("tempfit1", "pol1", gWidthCalibXMin, gWidthCalibXMax);
            T1->Fit(tempFit1, "QN");  // Q = quiet, N = don't store function in graph
            sett->widthCal[0][0] = tempFit1->GetParameter(0);
            sett->widthCal[0][1] = tempFit1->GetParameter(1);
            delete tempFit1;
        }
        
        if (T2->GetN() > 0) {
            TF1 *tempFit2 = new TF1("tempfit2", "pol1", gWidthCalibXMin, gWidthCalibXMax);
            T2->Fit(tempFit2, "QN");  // Q = quiet, N = don't store function in graph
            sett->widthCal[1][0] = tempFit2->GetParameter(0);
            sett->widthCal[1][1] = tempFit2->GetParameter(1);
            delete tempFit2;
        }
        
        // Send the new parameters to UI
        std::string msg = "WIDTH_CALIB_PARAMS:";
        msg += std::to_string(sett->widthCal[0][0]) + "|" + std::to_string(sett->widthCal[0][1]) + "|";
        msg += std::to_string(sett->widthCal[1][0]) + "|" + std::to_string(sett->widthCal[1][1]);
        window->Send(connid, msg);
        
        // Redraw the plot with the new fit lines
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
        
        // Draw first graph with axes
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
        
        // Draw fit lines with new parameters
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
        window->Send(connid, "Width calibration fits refreshed.");
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
        
        std::string rhoPath = sett->rhoFileName;
        
        // Check if file exists
        if (gSystem->AccessPathName(rhoPath.c_str())) {
            window->Send(connid, "Level density file not found: " + rhoPath);
            return;
        }
        
        // Read level density data through ShapeRho to get proper MeV→keV conversion
        ShapeRho* tempRho = new ShapeRho(sett);
        TGraphErrors *grLevelDensity = tempRho->rhoGraph;
        
        if (!grLevelDensity || grLevelDensity->GetN() == 0) {
            delete tempRho;
            window->Send(connid, "No level density data in file: " + rhoPath);
            return;
        }
        
        delete tempRho;
        
        grLevelDensity->SetMarkerStyle(20);
        grLevelDensity->SetMarkerSize(0.8);
        grLevelDensity->SetMarkerColor(kBlue);
        grLevelDensity->SetLineColor(kBlue);
        grLevelDensity->SetTitle("Level Density");
        
        gDisplayMode = 0;
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
        gCurrentHist = nullptr;
        
        grLevelDensity->Draw("APE");
        grLevelDensity->GetXaxis()->SetTitle("Excitation Energy (keV)");
        grLevelDensity->GetYaxis()->SetTitle("Level Density (1/keV)");
        
        window->Send(connid, "Displaying level density: " + std::to_string(grLevelDensity->GetN()) + " points");
        
        PushCanvasUpdate();
    }
    else if (arg == "SHOW_DISCRETE_LEVELS") {
        // Check if discrete level file is set
        if (sett->discreteLevelFile.empty()) {
            window->Send(connid, "No discrete level file specified in settings. Load one via Settings > Load discrete levels data...");
            return;
        }
        
        // Check if file exists
        if (gSystem->AccessPathName(sett->discreteLevelFile.c_str())) {
            window->Send(connid, "Discrete level file not found: " + sett->discreteLevelFile);
            return;
        }
        
        // Read the discrete levels data file (2 column: Ex [keV], number_of_states)
        std::ifstream infile(sett->discreteLevelFile.c_str());
        if (!infile.is_open()) {
            window->Send(connid, "ERROR: Could not open discrete level file: " + sett->discreteLevelFile);
            return;
        }
        
        std::vector<double> energies;
        std::vector<double> states;
        
        double ex, nstates;
        while (infile >> ex >> nstates) {
            energies.push_back(ex);
            states.push_back(nstates);
        }
        infile.close();
        
        if (energies.empty()) {
            window->Send(connid, "ERROR: No data found in discrete level file");
            return;
        }
        
        // Determine bin width from data spacing
        double binWidth = sett->discreteBins;  // Default from settings
        if (energies.size() > 1) {
            // Use spacing between first two points as bin width
            binWidth = energies[1] - energies[0];
        }
        
        // Create histogram
        double eMin = energies.front() - binWidth/2.0;
        double eMax = energies.back() + binWidth/2.0;
        int nBins = (int)((eMax - eMin) / binWidth + 0.5);
        
        TH1D* hDiscrete = new TH1D("hDiscrete", "Discrete Level Density;E_{x} [keV];#rho [keV^{-1}]", 
                                    nBins, eMin, eMax);
        
        // Fill histogram
        for (size_t i = 0; i < energies.size(); i++) {
            // Convert from number_of_states to level density (divide by bin width in keV)
            double rho = states[i] / binWidth;
            hDiscrete->Fill(energies[i], rho);
        }
        
        // Draw on canvas
        gDisplayMode = 0;
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
        gCurrentHist = nullptr;
        
        hDiscrete->SetLineColor(kBlue);
        hDiscrete->SetLineWidth(2);
        hDiscrete->Draw("HIST");
        canvas->SetLogy();
        
        std::ostringstream msg;
        msg << "Displayed discrete levels from: " << sett->discreteLevelFile 
            << " (" << energies.size() << " bins, " << binWidth << " keV bin width)";
        window->Send(connid, msg.str());
        
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
    else if (starts_with(arg, "SHOWBINPROJ:")) {
        if (!matrix) {
            window->Send(connid, "No matrix loaded yet -- open one first.");
            return;
        }
        
        // Save only X-axis zoom state before switching bins
        // Let Y-axis auto-scale for each new bin (different bins have different count ranges)
        bool hadRange = gHaveLastRange;
        double savedXmin = gLastUxmin;
        double savedXmax = gLastUxmax;
        
        gCurrentBin = std::stoi(after_prefix(arg, "SHOWBINPROJ:"));
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
    else if (starts_with(arg, "EXCITATION:")) {
        auto v = ParsePipeDoubles(after_prefix(arg, "EXCITATION:"));
        if (v.size() == 2) {
            sett->exiEne[0] = v[0];
            sett->exiEne[1] = v[1];
            SendNBins(connid);
        }
    }
    else if (starts_with(arg, "LEVELENERGIES:")) {
        std::cout << "*** LEVELENERGIES HANDLER CALLED ***" << std::endl;
        auto v = ParsePipeDoubles(after_prefix(arg, "LEVELENERGIES:"));
        std::cout << "*** Parsed " << v.size() << " values ***" << std::endl;
        if (v.size() != 12) {
            std::cout << "*** ERROR: Expected 12 values, got " << v.size() << " ***" << std::endl;
            return;
        }
        std::cout << "*** Setting levEne[0-3] to: " << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << " ***" << std::endl;
        
        // Store previous multiplet checkbox states to detect changes
        bool hadDoublet1 = sett->doDoublet[0];
        bool hadDoublet2 = sett->doDoublet[1];
        bool hadTriplet1 = sett->doTriplet[0];
        bool hadTriplet2 = sett->doTriplet[1];
        bool hadFixDoubletWidth1 = sett->fixDoubletWidth[0];
        bool hadFixDoubletWidth2 = sett->fixDoubletWidth[1];
        bool hadFixTripletWidth1 = sett->fixTripletWidth[0];
        bool hadFixTripletWidth2 = sett->fixTripletWidth[1];
        
        // Update main peak energies (indices 0-3)
        sett->levEne[0] = v[0]; 
        sett->levEne[1] = v[1];
        sett->levEne[2] = v[2]; 
        sett->levEne[3] = v[3];
        
        // Update doublet checkbox states (indices 4, 5)
        sett->doDoublet[0] = v[4] != 0.0;
        sett->doDoublet[1] = v[5] != 0.0;
        
        // Update doublet width fix toggles (indices 6-7)
        sett->fixDoubletWidth[0] = v[6] != 0.0;
        sett->fixDoubletWidth[1] = v[7] != 0.0;
        
        // Update triplet checkbox states (indices 8, 9)
        sett->doTriplet[0] = v[8] != 0.0;
        sett->doTriplet[1] = v[9] != 0.0;
        
        // Update triplet width fix toggles (indices 10-11)
        sett->fixTripletWidth[0] = v[10] != 0.0;
        sett->fixTripletWidth[1] = v[11] != 0.0;
        
        std::cout << "*** Level energies updated successfully ***" << std::endl;
        std::cout << "*** Doublet 1: " << (sett->doDoublet[0] ? "ENABLED" : "DISABLED") 
                  << ", fix width: " << (sett->fixDoubletWidth[0] ? "YES" : "NO") << " ***" << std::endl;
        std::cout << "*** Doublet 2: " << (sett->doDoublet[1] ? "ENABLED" : "DISABLED") 
                  << ", fix width: " << (sett->fixDoubletWidth[1] ? "YES" : "NO") << " ***" << std::endl;
        std::cout << "*** Triplet 1: " << (sett->doTriplet[0] ? "ENABLED" : "DISABLED") 
                  << ", fix width: " << (sett->fixTripletWidth[0] ? "YES" : "NO") << " ***" << std::endl;
        std::cout << "*** Triplet 2: " << (sett->doTriplet[1] ? "ENABLED" : "DISABLED") 
                  << ", fix width: " << (sett->fixTripletWidth[1] ? "YES" : "NO") << " ***" << std::endl;
        
        // Detect if multiplet state or width fix state changed
        bool multipletStateChanged = (hadDoublet1 != sett->doDoublet[0]) || (hadDoublet2 != sett->doDoublet[1]) ||
                                     (hadTriplet1 != sett->doTriplet[0]) || (hadTriplet2 != sett->doTriplet[1]);
        bool widthFixChanged = (hadFixDoubletWidth1 != sett->fixDoubletWidth[0]) || 
                               (hadFixDoubletWidth2 != sett->fixDoubletWidth[1]) ||
                               (hadFixTripletWidth1 != sett->fixTripletWidth[0]) || 
                               (hadFixTripletWidth2 != sett->fixTripletWidth[1]);
        
        // If in Autofit mode and viewing a bin projection, ALWAYS re-fit when any level energy changes
        if (sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
            std::cout << "Multiplet settings changed in Autofit mode: re-fitting bin " << gCurrentBin << "..." << std::endl;
            
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
    else if (starts_with(arg, "PEAKPOS:")) {
        // Format: PEAKPOS:level|enable|position
        // level: 0=level1, 1=level2
        // enable: 0=off, 1=on
        // position: peak position in keV
        auto v = ParsePipeDoubles(after_prefix(arg, "PEAKPOS:"));
        if (v.size() != 3) {
            window->Send(connid, "Malformed PEAKPOS message.");
            return;
        }
        
        int level = (int)v[0];
        bool enable = v[1] != 0.0;
        double position = v[2];
        
        if (level < 0 || level > 1) {
            window->Send(connid, "Invalid level in PEAKPOS message.");
            return;
        }
        
        sett->fixPeakPos[level] = enable;
        sett->peakPos[level] = position;
        
        std::cout << "Peak position for level " << (level+1) << ": " 
                  << (enable ? "ENABLED" : "DISABLED") 
                  << " at " << position << " keV" << std::endl;
        
        // If in Autofit mode and viewing a bin projection, ALWAYS re-fit when peak position changes
        if (sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
            std::cout << "Peak position changed in Autofit mode: re-fitting bin " << gCurrentBin << "..." << std::endl;
            
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
            // Just redraw markers
            DrawMarkers(true);
        }
    }
    else if (starts_with(arg, "DOUBLETPEAKPOS:")) {
        // Format: DOUBLETPEAKPOS:level|enable|position
        // level: 0=level1_doublet, 1=level2_doublet
        // enable: 0=off, 1=on
        // position: peak position in keV
        auto v = ParsePipeDoubles(after_prefix(arg, "DOUBLETPEAKPOS:"));
        if (v.size() != 3) {
            window->Send(connid, "Malformed DOUBLETPEAKPOS message.");
            return;
        }
        
        int level = (int)v[0];
        bool enable = v[1] != 0.0;
        double position = v[2];
        
        if (level < 0 || level > 1) {
            window->Send(connid, "Invalid level in DOUBLETPEAKPOS message.");
            return;
        }
        
        sett->fixDoubletPeakPos[level] = enable;
        sett->doubletPeakPos[level] = position;
        
        std::cout << "Doublet peak position for level " << (level+1) << ": " 
                  << (enable ? "ENABLED" : "DISABLED") 
                  << " at " << position << " keV" << std::endl;
        
        // If in Autofit mode and viewing a bin projection, ALWAYS re-fit when doublet peak position changes
        if (sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
            std::cout << "Doublet peak position changed in Autofit mode: re-fitting bin " << gCurrentBin << "..." << std::endl;
            
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
            // Just redraw markers
            DrawMarkers(true);
        }
    }
    else if (starts_with(arg, "TRIPLETPEAKPOS:")) {
        // Format: TRIPLETPEAKPOS:level|enable|position
        // level: 0=level1_triplet, 1=level2_triplet
        // enable: 0=off, 1=on
        // position: peak position in keV
        auto v = ParsePipeDoubles(after_prefix(arg, "TRIPLETPEAKPOS:"));
        if (v.size() != 3) {
            window->Send(connid, "Malformed TRIPLETPEAKPOS message.");
            return;
        }
        
        int level = (int)v[0];
        bool enable = v[1] != 0.0;
        double position = v[2];
        
        if (level < 0 || level > 1) {
            window->Send(connid, "Invalid level in TRIPLETPEAKPOS message.");
            return;
        }
        
        sett->fixTripletPeakPos[level] = enable;
        sett->tripletPeakPos[level] = position;
        
        std::cout << "Triplet peak position for level " << (level+1) << ": " 
                  << (enable ? "ENABLED" : "DISABLED") 
                  << " at " << position << " keV" << std::endl;
        
        // If in Autofit mode and viewing a bin projection, ALWAYS re-fit when triplet peak position changes
        if (sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
            std::cout << "Triplet peak position changed in Autofit mode: re-fitting bin " << gCurrentBin << "..." << std::endl;
            
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
            // Just redraw markers
            DrawMarkers(true);
        }
    }
    else if (starts_with(arg, "NEW_PEAK_CONFIG:")) {
        // Format: NEW_PEAK_CONFIG:enablePeak2L1|fixPeak2L1|peakPos2L1|enablePeak3L1|fixPeak3L1|peakPos3L1|keepWidths1|
        //                          enablePeak2L2|fixPeak2L2|peakPos2L2|enablePeak3L2|fixPeak3L2|peakPos3L2|keepWidths2
        // This maps to:
        //   - doDoublet[0-1], fixDoubletPeakPos[0-1], doubletPeakPos[0-1]
        //   - doTriplet[0-1], fixTripletPeakPos[0-1], tripletPeakPos[0-1]
        //   - fixDoubletWidth[0-1] (via keepWidths flags)
        auto v = ParsePipeDoubles(after_prefix(arg, "NEW_PEAK_CONFIG:"));
        if (v.size() != 14) {
            std::cout << "ERROR: NEW_PEAK_CONFIG expected 14 values, got " << v.size() << std::endl;
            window->Send(connid, "Malformed NEW_PEAK_CONFIG message.");
            return;
        }
        
        std::cout << "Processing NEW_PEAK_CONFIG message..." << std::endl;
        
        // Level 1 Peak 2 (doublet)
        sett->doDoublet[0] = v[0] != 0.0;
        sett->fixDoubletPeakPos[0] = v[1] != 0.0;
        sett->doubletPeakPos[0] = v[2];
        
        // Level 1 Peak 3 (triplet)
        sett->doTriplet[0] = v[3] != 0.0;
        sett->fixTripletPeakPos[0] = v[4] != 0.0;
        sett->tripletPeakPos[0] = v[5];
        
        // Level 1 keep widths equal
        sett->fixDoubletWidth[0] = v[6] != 0.0;
        sett->fixTripletWidth[0] = v[6] != 0.0;  // Both controlled by same checkbox
        
        // Level 2 Peak 2 (doublet)
        sett->doDoublet[1] = v[7] != 0.0;
        sett->fixDoubletPeakPos[1] = v[8] != 0.0;
        sett->doubletPeakPos[1] = v[9];
        
        // Level 2 Peak 3 (triplet)
        sett->doTriplet[1] = v[10] != 0.0;
        sett->fixTripletPeakPos[1] = v[11] != 0.0;
        sett->tripletPeakPos[1] = v[12];
        
        // Level 2 keep widths equal
        sett->fixDoubletWidth[1] = v[13] != 0.0;
        sett->fixTripletWidth[1] = v[13] != 0.0;  // Both controlled by same checkbox
        
        std::cout << "NEW_PEAK_CONFIG applied:" << std::endl;
        std::cout << "  Level 1 Peak 2: " << (sett->doDoublet[0] ? "ON" : "OFF") 
                  << " at " << sett->doubletPeakPos[0] << " keV" 
                  << (sett->fixDoubletPeakPos[0] ? " (FIXED)" : " (FREE)") << std::endl;
        std::cout << "  Level 1 Peak 3: " << (sett->doTriplet[0] ? "ON" : "OFF") 
                  << " at " << sett->tripletPeakPos[0] << " keV"
                  << (sett->fixTripletPeakPos[0] ? " (FIXED)" : " (FREE)") << std::endl;
        std::cout << "  Level 2 Peak 2: " << (sett->doDoublet[1] ? "ON" : "OFF") 
                  << " at " << sett->doubletPeakPos[1] << " keV"
                  << (sett->fixDoubletPeakPos[1] ? " (FIXED)" : " (FREE)") << std::endl;
        std::cout << "  Level 2 Peak 3: " << (sett->doTriplet[1] ? "ON" : "OFF") 
                  << " at " << sett->tripletPeakPos[1] << " keV"
                  << (sett->fixTripletPeakPos[1] ? " (FIXED)" : " (FREE)") << std::endl;
        
        // If in Autofit mode and viewing a bin projection, refit with new peak configuration
        if (sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
            std::cout << "NEW_PEAK_CONFIG in Autofit mode: re-fitting bin " << gCurrentBin << "..." << std::endl;
            
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
            // Just redraw markers if not in the right mode/view
            DrawMarkers(true);
        }
    }
    else if (starts_with(arg, "WIDTH_CALIB_PARAMS:")) {
        // order: l1_offset|l1_slope|l2_offset|l2_slope
        auto v = ParsePipeDoubles(after_prefix(arg, "WIDTH_CALIB_PARAMS:"));
        if (v.size() != 4) {
            window->Send(connid, "Malformed WIDTH_CALIB_PARAMS message.");
            return;
        }
        sett->widthCal[0][0] = v[0];  // Level 1 offset
        sett->widthCal[0][1] = v[1];  // Level 1 slope
        sett->widthCal[1][0] = v[2];  // Level 2 offset
        sett->widthCal[1][1] = v[3];  // Level 2 slope
    }
    else if (starts_with(arg, "UPDATE_WIDTH_CALIB_LINES:")) {
        // Just update the fit line parameters without re-running analysis
        auto v = ParsePipeDoubles(after_prefix(arg, "UPDATE_WIDTH_CALIB_LINES:"));
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
    else if (starts_with(arg, "SAVE_WIDTH_CALIB:")) {
        std::string path = after_prefix(arg, "SAVE_WIDTH_CALIB:");
        sett->settFileName = path;
        sett->SaveSettings();
        window->Send(connid, "Width calibration saved to: " + path);
    }
    else if (starts_with(arg, "ALPHA_TRANSFORM:")) {
        // Expected format: alpha|lit_norm
        auto v = ParsePipeDoubles(after_prefix(arg, "ALPHA_TRANSFORM:"));
        if (v.size() != 2) {
            window->Send(connid, "Malformed ALPHA_TRANSFORM message.");
            return;
        }
        
        sett->lit_alpha = v[0];
        sett->lit_norm = v[1];
        
        // If we have a collector (after ShapeIt has been run), apply the transformation
        if (gSFColl) {
            std::cout << "Applying alpha transformation: alpha=" << sett->lit_alpha 
                      << ", norm=" << sett->lit_norm << std::endl;
            gSFColl->Transform(sett->lit_norm, sett->lit_alpha);
            
            // Force complete redraw by clearing and rebuilding the graph display
            // This is the same logic as in RunShapeIt but without re-running the analysis
            canvas->cd();
            canvas->Clear();
            for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }
            
            // Extract and redraw the transformed graphs
            TMultiGraph *diagGraph = gSFColl->getMultGraph();
            TList *graphList = diagGraph->GetListOfGraphs();
            
            TGraph *litGraph = nullptr;
            if ((sett->doOslo || sett->doMC) && !sett->osloFileName.empty()) {
                litGraph = gSFColl->getLitGraph();
            }
            TGraph *avgGraph = sett->displayAvg ? gSFColl->getAvgGraph() : nullptr;
            
            bool firstDrawn = false;
            TGraph *firstGraph = nullptr;
            int colorIdx = 0;
            int color1 = 6;
            int color2 = sett->colour ? 7 : 6;
            
            if (graphList) {
                TIter next(graphList);
                TObject *obj;
                while ((obj = next())) {
                    TGraph *g = dynamic_cast<TGraph *>(obj);
                    if (!g) continue;
                    
                    bool isLit = (litGraph != nullptr && g == litGraph);
                    bool isAvg = (avgGraph != nullptr && g == avgGraph);
                    
                    if (isLit) {
                        g->SetFillColor(kBlue - 10);
                        g->SetFillStyle(3013);
                        g->SetLineColor(kBlue);
                        g->SetLineWidth(2);
                        g->Draw(firstDrawn ? "L3 SAME" : "AL3");
                        if (!firstDrawn) firstGraph = g;
                    } else if (isAvg) {
                        g->SetMarkerStyle(22);
                        g->SetMarkerSize(2);
                        g->SetMarkerColor(1);
                        g->SetLineColor(1);
                        g->Draw(firstDrawn ? "P SAME" : "AP");
                        if (!firstDrawn) firstGraph = g;
                    } else {
                        g->SetMarkerStyle(22);
                        g->SetMarkerSize(2);
                        int color = (colorIdx % 2 == 0) ? color1 : color2;
                        g->SetMarkerColor(color);
                        g->SetLineColor(color);
                        g->Draw(firstDrawn ? "P SAME" : "AP");
                        if (!firstDrawn) firstGraph = g;
                        colorIdx++;
                    }
                    firstDrawn = true;
                }
            }
            
            // Set axis labels and reset Y-axis range for new data
            if (firstGraph) {
                // Recompute axis ranges from ALL graphs to ensure everything is visible
                double xmin = 1e99, xmax = -1e99, ymin = 1e99, ymax = -1e99;
                
                if (graphList) {
                    TIter next(graphList);
                    TObject *obj;
                    while ((obj = next())) {
                        TGraph *g = dynamic_cast<TGraph *>(obj);
                        if (!g || g->GetN() == 0) continue;
                        
                        double gxmin, gxmax, gymin, gymax;
                        g->ComputeRange(gxmin, gymin, gxmax, gymax);
                        xmin = std::min(xmin, gxmin);
                        xmax = std::max(xmax, gxmax);
                        ymin = std::min(ymin, gymin);
                        ymax = std::max(ymax, gymax);
                    }
                }
                
                TH1F *hist = firstGraph->GetHistogram();
                if (hist) {
                    hist->SetTitle("Gamma Ray Strength Function from Shape Method;E_{#gamma} (keV);f(E_{#gamma}) (MeV^{-3})");
                    hist->GetXaxis()->SetTitle("E_{#gamma} (keV)");
                    hist->GetYaxis()->SetTitle("f(E_{#gamma}) (MeV^{-3})");
                    hist->GetXaxis()->SetTitleSize(0.04);
                    hist->GetYaxis()->SetTitleSize(0.04);
                    hist->GetXaxis()->SetTitleOffset(1.0);
                    hist->GetYaxis()->SetTitleOffset(1.2);
                    
                    // Set Y-axis range based on actual data with proper padding
                    // For lower bound, use larger of absolute padding or zero (avoid negative on log scale)
                    double yrange = ymax - ymin;
                    double yminPadded = ymin - 0.1*yrange;  // 10% padding below
                    if (yminPadded < 0 || ymin <= 0) yminPadded = ymin * 0.5;  // For log scale, use 50% of min
                    hist->SetMinimum(yminPadded);
                    hist->SetMaximum(ymax * 1.1);  // 10% above max
                    
                    // Don't change log scale - preserve whatever the user had set
                    gPad->Modified();
                }
            }
            
            // Add chi2 and alpha text box if Oslo data is displayed
            if (sett->doOslo) {
                TPaveText *t = new TPaveText(0.8, 0.85, 0.95, 0.95, "brNDC");
                t->SetTextSize(0.025);
                t->SetTextAlign(13);
                t->SetFillColor(10);
                t->SetTextColor(61);
                t->AddText(Form("slope #alpha: %4.2f", sett->lit_alpha));
                t->AddText(Form("#chi^{2} value: %4.2f", gSFColl->getChi2()));
                t->Draw();
            }
            
            gDisplayMode = 0;
            gHaveLastRange = false;
            PushCanvasUpdate();
        }
    }
    else if (arg == "REDRAW_GRAPH") {
        std::cout << "DEBUG: REDRAW_GRAPH handler reached" << std::endl;
        std::cout.flush();
        
        // Redraw the gSF graph with current transformation (mirrors TransGraph() from native code)
        if (!gSFColl) {
            window->Send(connid, "Run ShapeIt first before applying slope correction.");
            return;
        }
        
        std::cout << "Redrawing graph with current transformation..." << std::endl;
        
        // Get fresh graphs from collector (already transformed)
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
        TGraph *litGraph = nullptr;
        if ((sett->doOslo || sett->doMC) && !sett->osloFileName.empty()) {
            litGraph = gSFColl->getLitGraph();
        }
        
        // Get the average/smoothed graph pointer to identify it
        TGraph *avgGraph = sett->displayAvg ? gSFColl->getAvgGraph() : nullptr;

        if (graphList) {
            TIter next(graphList);
            TObject *obj;
            while ((obj = next())) {
                TGraph *g = dynamic_cast<TGraph *>(obj);
                if (!g) continue;

                bool isLit = (litGraph != nullptr && g == litGraph);
                bool isAvg = (avgGraph != nullptr && g == avgGraph);

                if (auto *ge = dynamic_cast<TGraphAsymmErrors *>(g)) {
                    std::vector<double> x, y, exl, exh, eyl, eyh;
                    for (int i = 0; i < ge->GetN(); i++) {
                        x.push_back(ge->GetX()[i]); y.push_back(ge->GetY()[i]);
                        exl.push_back(ge->GetEXlow()[i]); exh.push_back(ge->GetEXhigh()[i]);
                        eyl.push_back(ge->GetEYlow()[i]); eyh.push_back(ge->GetEYhigh()[i]);
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
                        allX.push_back(x.back()); allY.push_back(y.back());
                    }
                    auto *fresh = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
                    fresh->Sort();
                    freshGraphs.push_back({fresh, isLit, isAvg});
                }
                else {
                    for (int i = 0; i < g->GetN(); i++) {
                        allX.push_back(g->GetX()[i]); allY.push_back(g->GetY()[i]);
                    }
                }
            }
        }

        // Clear and redraw
        canvas->cd();
        canvas->Clear();
        for (int i = 0; i < 4; i++) { gMarkerLine[i] = nullptr; gDoubletLine[i] = nullptr; gBgBox[i] = nullptr; }

        bool firstDrawn = false;
        TGraph *firstGraph = nullptr;
        int colorIdx = 0;
        int color1 = 6;  // kMagenta
        int color2 = sett->colour ? 7 : 6;  // kCyan if colour enabled
        
        for (auto &fg : freshGraphs) {
            TGraph *g = fg.graph;

            if (fg.isLiterature) {
                g->SetFillColor(kBlue - 10);
                g->SetFillStyle(3013);
                g->SetLineColor(kBlue);
                g->SetLineWidth(2);
                g->Draw(firstDrawn ? "L3 SAME" : "AL3");
                if (!firstDrawn) firstGraph = g;
            } else if (fg.isAverage) {
                g->SetMarkerStyle(22);
                g->SetMarkerSize(2);
                g->SetMarkerColor(1);
                g->SetLineColor(1);
                g->Draw(firstDrawn ? "P SAME" : "AP");
                if (!firstDrawn) firstGraph = g;
            } else {
                g->SetMarkerStyle(22);
                g->SetMarkerSize(2);
                int color = (colorIdx % 2 == 0) ? color1 : color2;
                g->SetMarkerColor(color);
                g->SetLineColor(color);
                g->Draw(firstDrawn ? "P SAME" : "AP");
                if (!firstDrawn) firstGraph = g;
                colorIdx++;
            }
            
            firstDrawn = true;
        }
        
        // Set axis labels
        if (firstGraph) {
            TH1F *hist = firstGraph->GetHistogram();
            if (hist) {
                hist->SetTitle("Gamma Ray Strength Function from Shape Method;E_{#gamma} (keV);f(E_{#gamma}) (MeV^{-3})");
                hist->GetXaxis()->SetTitle("E_{#gamma} (keV)");
                hist->GetYaxis()->SetTitle("f(E_{#gamma}) (MeV^{-3})");
                hist->GetXaxis()->SetTitleSize(0.04);
                hist->GetYaxis()->SetTitleSize(0.04);
                hist->GetXaxis()->SetTitleOffset(1.0);
                hist->GetYaxis()->SetTitleOffset(1.2);
                firstGraph->SetTitle("Gamma Ray Strength Function from Shape Method");
                gPad->Modified();
            }
        }
        
        // Add chi2 and alpha text box if Oslo data is displayed
        if (sett->doOslo) {
            TPaveText *t = new TPaveText(0.8, 0.85, 0.95, 0.95, "brNDC");
            t->SetTextSize(0.025);
            t->SetTextAlign(13);
            t->SetFillColor(10);
            t->SetTextColor(61);
            t->AddText(Form("slope #alpha: %4.2f", sett->lit_alpha));
            t->AddText(Form("#chi^{2} value: %4.2f", gSFColl->getChi2()));
            t->Draw();
        }
        
        gDisplayMode = 0;  // Results view
        gHaveLastRange = false;
        
        PushCanvasUpdate();
        window->Send(connid, "Graph redrawn with alpha=" + std::to_string(sett->lit_alpha));
    }
    else if (starts_with(arg, "RUN_MC_ALPHA:")) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "RUN_MC_ALPHA HANDLER ENTERED" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout.flush();  // Force immediate output
        
        std::cout << "Raw message: " << arg << std::endl;
        std::cout.flush();
        
        std::cout << "About to parse doubles..." << std::endl;
        std::cout.flush();
        
        auto v = ParsePipeDoubles(after_prefix(arg, "RUN_MC_ALPHA:"));
        
        std::cout << "Parsing complete. Got " << v.size() << " values" << std::endl;
        std::cout.flush();
        
        if (v.size() != 3) {
            std::cout << "ERROR: Expected 3 values, got " << v.size() << std::endl;
            std::cout.flush();
            window->Send(connid, "Malformed RUN_MC_ALPHA message.");
            return;
        }
        
        std::cout << "Extracting parameters..." << std::endl;
        std::cout.flush();
        
        int nIter = (int)v[0];
        double exiMin = v[1];
        double exiMax = v[2];
        
        std::cout << "Parameters extracted:" << std::endl;
        std::cout << "  nIter = " << nIter << std::endl;
        std::cout << "  exiMin = " << exiMin << std::endl;
        std::cout << "  exiMax = " << exiMax << std::endl;
        std::cout.flush();
        
        std::cout << "About to call RunMonteCarlo()..." << std::endl;
        std::cout.flush();
        
        RunMonteCarlo(connid, nIter, exiMin, exiMax);
        
        std::cout << "RunMonteCarlo() returned successfully" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout.flush();
    }
    else if (arg == "STOP_MC") {
        std::cout << "STOP_MC command received" << std::endl;
        if (gMCState) {
            std::cout << "Setting stopRequested flag on MC state" << std::endl;
            gMCState->stopRequested = true;
            window->Send(connid, "MC_STOPPING");
        } else {
            std::cout << "No active MC simulation to stop" << std::endl;
            window->Send(connid, "No active Monte Carlo simulation.");
        }
    }
    else if (starts_with(arg, "RUN:")) {
        // expected order: lvl1_lo|lvl1_hi|lvl2_lo|lvl2_hi|exc_lo|exc_hi|
        //                  is_doublet1|d1_lo|d1_hi|is_doublet2|d2_lo|d2_hi|fix_doublet_width1|fix_doublet_width2|
        //                  is_triplet1|t1_lo|t1_hi|is_triplet2|t2_lo|t2_hi|fix_triplet_width1|fix_triplet_width2
        auto v = ParsePipeDoubles(after_prefix(arg, "RUN:"));
        if (v.size() != 22) {
            window->Send(connid, "Malformed RUN message (expected 22 values with triplet support).");
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
        
        // Set triplet checkbox states
        sett->doTriplet[0] = v[14] != 0.0;
        sett->doTriplet[1] = v[17] != 0.0;
        
        // ALWAYS store triplet energy values regardless of checkbox state
        sett->levEne_3[0] = v[15];
        sett->levEne_3[1] = v[16];
        sett->levEne_3[2] = v[18];
        sett->levEne_3[3] = v[19];
        
        // Set triplet width fix toggles
        sett->fixTripletWidth[0] = v[20] != 0.0;
        sett->fixTripletWidth[1] = v[21] != 0.0;

        std::cout << "*** About to call RunShapeIt() ***" << std::endl;
        std::cout << "*** Doublet 1: " << (sett->doDoublet[0] ? "ENABLED" : "DISABLED") << " ***" << std::endl;
        std::cout << "*** Doublet 2: " << (sett->doDoublet[1] ? "ENABLED" : "DISABLED") << " ***" << std::endl;
        std::cout << "*** Triplet 1: " << (sett->doTriplet[0] ? "ENABLED" : "DISABLED") << " ***" << std::endl;
        std::cout << "*** Triplet 2: " << (sett->doTriplet[1] ? "ENABLED" : "DISABLED") << " ***" << std::endl;
        
        RunShapeIt(connid);
        
        std::cout << "*** RunShapeIt() returned ***" << std::endl;
    }
    else {
        // Message not recognized - silently ignore
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
            sett->rhoFileName = ResolveRelativeTo(settDir, sett->rhoFileName);
            sett->discreteLevelFile = ResolveRelativeTo(settDir, sett->discreteLevelFile);
            
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

    // Polls the pad's axis range and marker positions every 500ms
    // - CheckRangeChanged() redraws markers when zoom/pan changes
    // - CheckMarkersChanged() updates settings when markers are dragged
    // Longer interval (500ms vs 200ms) reduces interference with ROOT's own
    // context menu operations (like unzoom), which need time to complete
    static TTimer *pollTimer = new TTimer();
    pollTimer->Connect("Timeout()", 0, 0, "CheckRangeChanged()");
    pollTimer->Connect("Timeout()", 0, 0, "CheckMarkersChanged()");
    pollTimer->Start(500, kFALSE);

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
