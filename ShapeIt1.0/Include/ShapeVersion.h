/* ***********************************************************************
* Copyright (C) 2019-2026, Dennis Muecher.                               *
* All rights reserved.                                                   *
*                                                                        *
* This program is free software: you can redistribute it and/or modify   *
* it under the terms of the GNU General Public License as published by   *
* the Free Software Foundation, either version 3 of the License, or      *
* (at your option) any later version.                                    *
* You should have received a copy of the GNU General Public License      *
* along with this program. If not, see  http://www.gnu.org/licenses/.    *
*************************************************************************/

#ifndef SHAPEVERSION_H
#define SHAPEVERSION_H

// Single source of truth for ShapeIt's version. Update this file, commit it,
// then tag the commit to match (see RELEASING.md) -- that keeps the About
// dialog, the GitHub release/tag, and the Zenodo archive all pointing at the
// same version number instead of drifting independently.
//
// Versioning follows semver (semver.org): MAJOR.MINOR.PATCH
//   MAJOR: incompatible settings-file format or analysis-changing behavior
//   MINOR: new features, backwards compatible
//   PATCH: bug fixes, no new features
#define SHAPEIT_VERSION "1.2.0"
#define SHAPEIT_RELEASE_DATE "2026-07-05"
#define SHAPEIT_COPYRIGHT_YEARS "2019-2026"

#endif
