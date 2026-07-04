# 3. Choosing mode and binning

## Mode: Integration vs Autofit

In the **Mode** box, choose how each peak's intensity is extracted in every
excitation-energy bin:

- **Integration** (default) — sums the raw counts in the peak window
  (a simple rectangular sum over the bins).
- **Autofit** — fits each peak with a Gaussian on a background (see
  [The algorithm](../algorithm.md#peak-fitting) for the functional form),
  useful when peaks sit on a steeply varying background or overlap with
  nearby contaminant lines. The fitted Gaussian's integral is used as the
  peak's intensity, in place of the raw sum.

### How background subtraction works in each mode

When **Background subtraction** is enabled, ShapeIt performs the
Autofit-style Gaussian-on-background fit internally to obtain a background
estimate under each peak, *regardless* of which Mode is selected. What
differs between the two modes is how the **signal** is treated:

- In **Autofit** mode, the peak's reported intensity is the fitted
  Gaussian's integral, with the fitted background's integral subtracted.
- In **Integration** mode, the peak's reported intensity is the raw
  rectangular sum over the peak window, rescaled by the fraction of it that
  the background fit identifies as signal rather than background, i.e. by
  $(\text{raw sum} - \text{fitted background})/\text{raw sum}$.

Both modes therefore use the same underlying background fit when background
subtraction is on; Integration mode applies it as a correction factor to the
raw counts rather than replacing the counts with the fitted peak shape. If
**Background subtraction** is off, Integration mode uses the raw rectangular
sum with no correction.

## Binning the excitation-energy axis

In the **Integration bin** box:

| Field | Meaning |
|---|---|
| Bin size [keV] | Width of each excitation-energy slice |
| Nr. of Bins | Number of slices to analyze |
| Min. Counts | Minimum peak counts required to accept a bin (bins below this are skipped) — see below |
| gSF scaling / Auto | Overall normalization of the extracted γSF — see below |
| Eff. Corr. | Correction factor applied to Level 2's intensity — see below |

By default only a single bin size and a single bin count are used (the
second "Nr. of Bins" / "Bin size" fields are disabled). Enabling
**Bin size variation** (see [next step](options-and-sewing.md)) activates
the second pair of fields, letting ShapeIt repeat the extraction over a
range of bin sizes/bin counts as a systematic-uncertainty check.

## Min. Counts

If either peak's intensity in a bin falls below this threshold, the entire
bin is skipped — it doesn't contribute a γSF pair. This protects the
strength ratio (see [The algorithm](../algorithm.md)) against bins where a
peak isn't statistically significant above background, which would
otherwise produce a wildly unreliable ratio.

## gSF scaling and Auto

The **gSF scaling** field sets an overall multiplicative normalization
applied to your final extracted γSF curve — the *shape* of the sewn curve is
fixed by the sewing algorithm, but its absolute vertical placement on the
plot is still free until you set this.

- **Auto checked** (and literature γSF data loaded via
  **Settings → Load Literature Values gSF...**): ShapeIt determines this
  factor automatically, by a least-squares fit of your extracted γSF onto
  the literature curve over their overlapping energy range. The field
  becomes read-only and shows the resulting value.
- **Auto unchecked**, or no literature data loaded: you set the value by
  hand. Internally, without a literature anchor, ShapeIt instead normalizes
  multiple sliding-window/bin-size realizations to each other (using the
  very first realization as the reference) so they overlay consistently —
  see [Running ShapeIt and reading results](running-and-results.md#reading-the-averaged-plot).

## Eff. Corr.

Level 1 and Level 2 aren't always detected with the same relative
efficiency — for example if they have different spins, or due to detector
geometry effects. ShapeIt has two, independent ways to correct for this,
which combine when both are used:

- **Eff. Corr.** (the number field here) is a single, **energy-independent**
  factor that always multiplies **Level 2's** intensity before the ratio
  $R$ (see [The algorithm](../algorithm.md)) is formed.
- **Settings → Load Efficiency Calibration...** loads an optional,
  **energy-dependent** efficiency curve from a file, interpolated and
  applied to *both* levels at their respective γ-ray energies (falling back
  to a factor of 1 outside the energies covered by the file). You can
  visualize a loaded curve via **Display → Plot Efficiency Data from File**.

The recommended way to determine the scalar **Eff. Corr.** value: simulate
your experiment (e.g. with GEANT4) using a known input γSF, run that
simulated data through the same Shape-method analysis, and adjust
**Eff. Corr.** until ShapeIt's extracted γSF reproduces the input γSF you
simulated with.

Next: [Options and sewing](options-and-sewing.md).
