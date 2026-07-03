# API reference

ShapeIt is organized as a set of ROOT-interpreted C++ classes, one per file
under `ShapeIt1.0/Source`:

| Class / macro | File | Role |
|---|---|---|
| `ShapeFrame` | `ShapeFrame.C` | Main GUI window: menus, control panel, event handling |
| `ShapeSetting` | `ShapeSetting.C` | Holds and (de)serializes all analysis settings |
| `ShapeMatrix` | `ShapeMatrix.C` | Loads the ROOT file, holds the input `TH2`, builds the diagonalized Ex-binned projections |
| `ShapeFitFunction` | `ShapeFitFunction.C` | The Gaussian(s)-on-quadratic-background peak fit function |
| `ShapeCollector` | `ShapeCollector.C` | Drives extraction across all excitation-energy bins for one run |
| `ShapeGSF` | `ShapeGSF.C` | Builds the Level 1 / Level 2 gSF graphs and performs the sewing/merge |
| `ShapeRho` / `ShapeRhoCollector` | `ShapeRho.C`, `ShapeRhoCollector.C` | Level-density (ρ) extraction and collection |
| `ShapeAlpha` / `ShapeDialogAlpha` | `ShapeAlpha.C`, `ShapeDialogAlpha.C` | Monte Carlo transformation/uncertainty dialog |
| `ShapeTalys` | `ShapeTalys.C` | Comparison against TALYS-based level-density/gSF models |
| `ShapeInfo` | `ShapeInfo.C` | The "About" info panel |

!!! info "Auto-generated reference"
    This table is a manually curated map of the codebase. If you'd like a
    fully auto-generated, function-level API reference (e.g. via Doxygen),
    I'm happy to set that up as a follow-up — it works well alongside this
    narrative guide.
