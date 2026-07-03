# Input data format

ShapeIt reads a **two-dimensional ROOT histogram (`TH2`)** from a `.root`
file. Any `TH2`-derived histogram present in the file can be selected from
the matrix dropdown once the file is loaded.

## Axis convention

| Axis | Quantity | Unit |
|---|---|---|
| x | γ-ray energy, Eγ | keV |
| y | Excitation energy, Ex | keV |

!!! warning "Units matter"
    ShapeIt assumes the matrix is binned in **keV** by default. If your
    matrix is binned in MeV, use **Settings → Transformation...** to set the
    MeV scale factor (this multiplies your axis values by 1000 internally)
    before doing anything else — energy windows you type into the "Energies"
    panel are always interpreted in keV.

## What the matrix should contain

The matrix should be a **first-generation (primary) γ-ray matrix**
$P(E_\gamma, E_x)$ — for each excitation-energy bin (a row of constant Ex),
the projection onto Eγ shows the *primary* γ-ray transitions directly
depopulating that state. This is the same first-generation matrix used in
Oslo-method analyses.

Transitions feeding a specific low-lying final level $L_j$ (fixed
$E_x - E_\gamma = E_{Lj}$) form a **diagonal** in this matrix — see
[Defining the two levels](user-guide/defining-levels.md) for why this
matters and how ShapeIt uses it. You'll need at least two such diagonals,
corresponding to two known, well-resolved final levels, that are populated
with enough statistics across the excitation-energy range you want to
analyze.

## Multiple matrices in one file

A single ROOT file can contain several `TH2` histograms (e.g. different
gates, different runs, or Ex-squared/Ex-cubed variants used internally by
ShapeIt). All of them appear in the **Input Matrix** dropdown after loading
the file; pick the one you want to analyze.

Once your matrix is ready, continue to
[Loading a matrix](user-guide/loading-a-matrix.md).
