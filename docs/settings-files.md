# Settings files

ShapeIt can save your entire analysis configuration (matrix, energy windows,
mode, options, calibrations) to a plain-text settings file, and reload it
later — useful for reproducing an analysis or sharing a configuration with
a collaborator.

## Saving and loading

Under the **Settings** menu:

- **Load Settings...** — load a previously saved `.dat` settings file
- **Save Settings** — overwrite the currently loaded settings file
- **Save Settings as...** — save to a new file
- **Print Settings to terminal** — dump the current, active settings to the
  ROOT terminal without writing a file (handy for a quick sanity check)

## File format

The settings file is a simple whitespace-separated text format, one entry
per line: a keyword followed by its value(s). The main keys are:

| Key | Meaning |
|---|---|
| `MeV:` | Unit scale factor (1 = matrix already in keV, 1000 = matrix in MeV) |
| `mode` | Integration (1) or Autofit (2) |
| `level1`, `level2` | Window (low, high) on the Ex − Eγ axis bracketing final level $L_1$ / $L_2$ |
| `level1_2`, `level2_2` | Optional second window per level, for **doublet** peak fitting — see the walkthrough in [Defining the two levels](user-guide/defining-levels.md#fitting-a-peak-as-a-doublet). Leave at `0 0` for a single peak. |
| `bg_level1`, `bg_level2` | Four background-window edges (same Ex − Eγ axis) for each level |
| `excitation` | Ex range (low, high) analyzed |
| `excitationMCLimits` | Range for the randomized lower Ex edge in the Monte Carlo procedure (see [Monte Carlo uncertainty](user-guide/monte-carlo-uncertainty.md)) |
| `excitation_bin_1`, `excitation_bin_2` | Bin size(s) [keV] |
| `nOfBins` | Number of excitation-energy bins |
| `minCounts` | Minimum peak counts to accept a bin |
| `eff_corr` | Efficiency-correction factor |
| `alphaLimit` | Trial slope-correction ($\alpha$) range for the χ² search and Monte Carlo simulation |
| `alphaIter` | Number of steps in the χ² search over that $\alpha$ range |
| `doInterpol` | Sewing interpolation on/off |
| `doBackground` | Background subtraction on/off |
| `doSlidingWindow` | Sliding-window variation on/off |
| `doBinVariation` | Bin-size variation on/off |
| `doWidthCal` | Use width calibration on/off |
| `doOslo` | Overlay literature/expectation values on/off |
| `doAutoScale` | Automatic gSF normalization on/off |
| `gSF_norm`, `lit_norm` | Normalization values |
| `widthCal` | Stored peak-width calibration parameters |
| `rhoFileName`, `rhoScale` | Level-density comparison file and scale |
| `discreteLevelFile`, `discreteBins` | Discrete level scheme file, for the level-density panel |

!!! tip
    You don't need to hand-edit these files — **Save Settings as...** writes
    a complete, valid file from your current GUI state. Hand-editing is
    mainly useful for scripting a batch of related analyses.

## Related loaders and their file formats

Under **Settings**, you can separately load three kinds of external data
file. Each is a plain-text file with one data point per line and
whitespace-separated columns — but the **column meaning differs for each**,
which is a common source of confusion. The exact formats are below.

### Load Literature Values gSF...

The comparison γSF (e.g. from an Oslo-method analysis of the same data),
overlaid when **Display expectation** is enabled and used as the reference
for the slope determination.

Three columns per line:

```
E_gamma    gSF_high    gSF_low
```

| Column | Meaning |
|---|---|
| `E_gamma` | γ-ray energy, in keV |
| `gSF_high` | Upper edge of the γSF value (value + its uncertainty) |
| `gSF_low` | Lower edge of the γSF value (value − its uncertainty) |

!!! warning "This file gives an error band, not value + error"
    ShapeIt does **not** read a central value and a symmetric error. It reads
    the **high and low edges** of each point, then internally takes the
    central value as $(\text{high}+\text{low})/2$ and the 1σ uncertainty as
    $(\text{high}-\text{low})/2$. If you have value ± error instead, convert
    to `value+error` and `value−error` before saving.

### Load Literature Values rho...

The comparison level density (e.g. an Oslo-method ρ for the same nucleus),
overlaid on the level-density panel.

Three columns per line:

```
energy    rho    rho_error
```

| Column | Meaning |
|---|---|
| `energy` | Excitation energy, in **MeV** |
| `rho` | Level density, in MeV⁻¹ |
| `rho_error` | Uncertainty on the level density, in MeV⁻¹ |

!!! note
    Unlike the γSF and efficiency files (which use keV), the level-density
    file's energy column is in **MeV**, matching the level-density plot's
    axis. The `rho` column is additionally multiplied by the `rhoScale`
    setting (default 1) when loaded.

### Load Efficiency Calibration...

This is the **energy-dependent** efficiency correction — yes, it's the
energy-dependent counterpart to the scalar **Eff. Corr.** field. When
loaded, it's interpolated and applied to *both* levels at their respective
γ-ray energies (returning a factor of 1 outside the energy range covered by
the file), and it combines with the scalar **Eff. Corr.** (which always
additionally multiplies Level 2). See
[Choosing mode and binning](user-guide/mode-and-binning.md#eff-corr).

Two columns per line:

```
E_gamma    factor
```

| Column | Meaning |
|---|---|
| `E_gamma` | γ-ray energy, in keV |
| `factor` | Multiplicative scaling factor applied to the peak intensity at that energy |

You can visualize a loaded curve with **Display → Plot Efficiency Data from
File**.

!!! tip "How the file paths are stored"
    Once loaded, the paths to these files are saved *inside* your settings
    file (`osloFileName`, `rhoFileName`, `effiFileName`), so re-loading a
    settings file automatically re-loads the associated data files — you
    don't have to re-select them each time.
