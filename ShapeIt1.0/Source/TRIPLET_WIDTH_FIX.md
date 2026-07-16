# Fix: Triplet Autofit with "Keep Widths Equal" 

## Problem
When using three peaks (triplet) with "keep widths equal" enabled, the second peak (doublet peak) was not fitting correctly. The fit worked fine when widths were allowed to vary independently, but failed when all peaks should share the same width.

## Root Cause
The parameter indexing in `ShapeMatrix.C` FitGauss() was **inconsistent** with the formula used in `ShapeFitFunction.C`.

### ShapeFitFunction Formula
In `ShapeFitFunction.C`, multiplet peaks use this indexing scheme:
```cpp
for (int j = 0; j < multip; j++) {
    // Amplitude: par[3*j + 6]
    // Position:  par[3*j + 7]  
    // Width:     par[3*j + 8] (or par[5] if fix_width=true)
}
```

This means:
- **Doublet (j=0)**: amp=6, pos=7, width=8 (or 5 if fixed)
- **Triplet (j=1)**: amp=9, pos=10, width=11 (or 5 if fixed)

The parameter indices are **always 9, 10, 11** for the triplet, regardless of whether widths are fixed or free.

### The Bug in ShapeMatrix.C
The code was calculating triplet parameter indices **conditionally** based on `fix_multiplet_width`:

```cpp
// OLD BUGGY CODE
int base_idx = fix_multiplet_width ? 8 : 9;  // WRONG!
fit_result[level]->SetParameter(base_idx, amplitude_init_3);      // Set amp
fit_result[level]->SetParameter(base_idx + 1, p_3);               // Set pos
```

This meant:
- When widths were fixed: tried to use parameters 8, 9 for triplet
- When widths were free: correctly used parameters 9, 10, 11 for triplet

**Result:** With fixed widths, the triplet amplitude was being set to parameter 8 (which should be the doublet's width parameter), causing the doublet fit to break.

## Solution
Fixed all triplet parameter references in `ShapeMatrix.C` to use **consistent indices**:
- Triplet amplitude: always parameter **9**
- Triplet position: always parameter **10**  
- Triplet width: always parameter **11** (only used when fix_multiplet_width=false)

Also corrected the `nParams` calculation to always allocate 3 parameters per multiplet peak, even when widths are shared (parameters exist but aren't used in the fit function when fix_width=true).

## Changes Made

### 1. Parameter Count (line ~342)
**Before:**
```cpp
nParams += fix_multiplet_width ? 2 : 3;  // For doublet
nParams += fix_multiplet_width ? 2 : 3;  // For triplet
```

**After:**
```cpp
nParams += 3;  // Always 3 for doublet
nParams += 3;  // Always 3 for triplet
// Parameters are always allocated even if width is shared
```

### 2. Background Fit - Zero Triplet Amplitude (line ~380)
**Before:**
```cpp
int base_idx = fix_multiplet_width ? 8 : 9;
fit_result[level]->FixParameter(base_idx, 0);
```

**After:**
```cpp
// Triplet amplitude is always at index 9
fit_result[level]->FixParameter(9, 0);
```

### 3. Set Initial Triplet Parameters (line ~411)
**Before:**
```cpp
int base_idx = fix_multiplet_width ? 8 : 9;
fit_result[level]->SetParameter(base_idx, amplitude_init_3);
fit_result[level]->SetParameter(base_idx + 1, p_3);
if (!fix_multiplet_width) {
    fit_result[level]->SetParameter(base_idx + 2, dp);
}
```

**After:**
```cpp
// Triplet amplitude at index 9, position at 10, width at 11 (if free)
fit_result[level]->SetParameter(9, amplitude_init_3);
fit_result[level]->SetParameter(10, p_3);
if (!fix_multiplet_width) {
    fit_result[level]->SetParameter(11, dp);
}
```

### 4. Set Triplet Parameter Limits (line ~437)
**Before:**
```cpp
int base_idx = fix_multiplet_width ? 8 : 9;
fit_result[level]->SetParLimits(base_idx, 0.01*amplitude_init_3, 100*amplitude_init_3);
fit_result[level]->SetParLimits(base_idx + 1, sett->levEne[2*level], sett->levEne[2*level+1]);
if (!fix_multiplet_width) {
    fit_result[level]->SetParLimits(base_idx + 2, 0.5*dp, 2*dp);
}
```

**After:**
```cpp
// Triplet amplitude at index 9, position at 10
fit_result[level]->SetParLimits(9, 0.01*amplitude_init_3, 100*amplitude_init_3);
fit_result[level]->SetParLimits(10, sett->levEne[2*level], sett->levEne[2*level+1]);
if (!fix_multiplet_width) {
    fit_result[level]->SetParLimits(11, 0.5*dp, 2*dp);
}
```

### 5. Fix Triplet Position (line ~480)
**Before:**
```cpp
int base_idx = fix_multiplet_width ? 9 : 10;  // BACKWARDS!
fit_result[level]->FixParameter(base_idx, sett->tripletPeakPos[level]);
```

**After:**
```cpp
// Triplet position is always at index 10
fit_result[level]->FixParameter(10, sett->tripletPeakPos[level]);
```

## Expected Behavior After Fix

### With "Keep Widths Equal" Enabled:
✅ All three peaks (main, doublet, triplet) share the **same width parameter (par[5])**  
✅ Doublet fits correctly at its specified position  
✅ Triplet fits correctly at its specified position  
✅ All Gaussians have identical widths

### With "Keep Widths Equal" Disabled:
✅ Main peak uses width parameter 5  
✅ Doublet uses independent width parameter 8  
✅ Triplet uses independent width parameter 11  
✅ Each Gaussian can have a different width

## Testing Recommendations

1. Load a matrix with three visible peaks in a bin projection
2. Enable Autofit mode
3. Add doublet and triplet peaks
4. **With "keep widths equal" OFF**: verify all three peaks fit independently
5. **With "keep widths equal" ON**: verify all three peaks fit with the same width
6. Check that the doublet peak position and amplitude are reasonable (not zero/wrong)
7. Verify width values are identical across all peaks when constraint is enabled

## Files Modified
- `ShapeMatrix.C` - Fixed triplet parameter indexing in FitGauss() function

---

**Date:** 2026-07-15  
**Issue:** Doublet peak not fitting correctly when "keep widths equal" enabled with triplet  
**Status:** Fixed
