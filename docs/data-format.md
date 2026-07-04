# Input data format

ShapeIt reads a **two-dimensional ROOT histogram (`TH2`)** from a `.root`
file. Any `TH2`-derived histogram present in the file can be selected from
the matrix dropdown once the file is loaded.

## Axis convention

| Axis | Quantity | Unit |
|---|---|---|
| x | γ-ray energy, Eγ | keV |
| y | Excitation energy, Ex | keV |

!!! warning "Units: keV is strongly recommended"
    ShapeIt works most naturally with a matrix binned in **keV**, and all the
    energy windows you type into the **Energies** panel are interpreted in
    keV. There is no GUI control to switch units — the only way to tell
    ShapeIt your matrix is in MeV is to set the `MeV:` key in the
    [settings file](settings-files.md#file-format) by hand
    (`MeV: 1` for a keV matrix, `MeV: 1000` for a MeV matrix), then re-read
    it with **Settings → Load Settings**. ShapeIt does **not** notice
    settings-file edits made in the background, so the reload is required for
    the change to take effect. Because this path is easy to forget, the
    simplest approach is to prepare your input matrix in keV in the first
    place.

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
