# Installation

## Requirements

- [ROOT](https://root.cern/) — developed and tested against ROOT 6.34
- A C++ compiler compatible with your ROOT installation (ShapeIt is compiled
  on the fly by ROOT's interpreter/JIT, so no separate build step is needed)

!!! note
    ShapeIt is distributed under the GNU General Public License, version 3
    or later. Bug reports and contributions are welcome via the
    [issue tracker](https://github.com/dennismuecher/ShapeIt/issues).

## Getting the code

```bash
git clone https://github.com/dennismuecher/ShapeIt.git
cd ShapeIt
git checkout ShapeItLatest
```

## Launching ShapeIt

ShapeIt is a ROOT GUI application, launched as a ROOT macro:

```bash
cd ShapeIt1.0/Source
root ShapeIt.C
```

This opens the main ShapeIt window. From here, everything is driven through
the menu bar and the control panel on the left — no command-line arguments
are needed.

Continue to [Input data format](data-format.md) to prepare your ROOT matrix,
or jump straight to the [user guide](user-guide/loading-a-matrix.md).
