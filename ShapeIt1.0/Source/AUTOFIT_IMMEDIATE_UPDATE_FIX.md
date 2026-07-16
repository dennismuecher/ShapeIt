# Fix: Immediate Autofit Updates When Side Panel Values Change

## Problem
When changing values in the side panel (peak positions, doublet/triplet settings, background regions), the autofit was **not updating immediately**. It would only update when you dragged the fit region markers with the mouse. This made the UI feel unresponsive and made it unclear whether your changes were taking effect.

## Root Cause
The message handlers in `WebShapeIt.cxx` had **overly restrictive conditions** for triggering a refit:

```c++
// OLD CODE - only refit when specific state changed
if ((multipletStateChanged || widthFixChanged) && sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
    // ... refit code ...
}
```

This meant:
- **Only** multiplet checkbox or width fix changes triggered refits
- Changes to peak positions, background regions, or fit region boundaries were **ignored**
- You had to manually drag a marker to force an update

## Solution
Modified all relevant message handlers to **always trigger a refit when in Autofit mode** (mode 2) and viewing a bin projection (gDisplayMode 5), regardless of what specific parameter changed:

```c++
// NEW CODE - always refit when in Autofit mode and viewing projection
if (sett->mode == 2 && gDisplayMode == 5 && gCurrentBin > 0) {
    // ... refit code ...
}
```

## Modified Handlers

### 1. `LEVELENERGIES:` Handler (line ~1775)
**Changed from:** Only refit when multiplet checkbox or width fix states changed  
**Changed to:** Always refit when any level energy value changes

**Impact:** Changing fit region boundaries now immediately updates the fit

### 2. `PEAKPOS:` Handler (line ~2526)
**Changed from:** Conditional refit  
**Changed to:** Always refit when peak position changes

**Impact:** Adjusting peak 1 position now immediately refits

### 3. `DOUBLETPEAKPOS:` Handler (line ~2563)
**Changed from:** Conditional refit  
**Changed to:** Always refit when doublet peak position changes

**Impact:** Adjusting peak 2 (doublet) position now immediately refits

### 4. `TRIPLETPEAKPOS:` Handler (line ~2600)
**Changed from:** Conditional refit  
**Changed to:** Always refit when triplet peak position changes

**Impact:** Adjusting peak 3 (triplet) position now immediately refits

### 5. `NEW_PEAK_CONFIG:` Handler (line ~2698)
**Changed from:** Only redraw markers, no refit  
**Changed to:** Refit when peak configuration changes

**Impact:** Toggling doublet/triplet checkboxes or changing "keep widths equal" now immediately refits

### 6. `BGENERGIES:` Handler (line ~1786)
**Changed from:** Only redraw markers, no refit  
**Changed to:** Refit when background regions change

**Impact:** Adjusting background regions now immediately refits

## Behavior After Fix

### In Autofit Mode (mode 2) + Viewing Bin Projection (gDisplayMode 5):
✅ **Changing any peak parameter → immediate autofit**
- Fit region boundaries (Level 1/2 low/high)
- Peak positions (Peak 1, Peak 2 doublet, Peak 3 triplet)
- Background regions (left/right for Level 1 and Level 2)
- Multiplet checkboxes (enable/disable doublet/triplet)
- Width fix toggles ("keep widths equal")

### In Integration Mode (mode 1) or Other Views:
✅ **Changes only update markers visually** (no performance impact)

## Testing Recommendations

1. **Load a matrix** and switch to **Projection view** (not Matrix view)
2. **Enable Autofit mode** (Mode panel)
3. **Change values** in the Levels panel:
   - Adjust fit region boundaries → should see fit update immediately
   - Change peak positions → should see Gaussian peaks shift
   - Toggle doublet/triplet checkboxes → should see extra peaks appear/disappear
   - Adjust background regions → should see background fit change

4. **Verify no performance issues:**
   - Changes should feel responsive (< 1 second for most bins)
   - No delays when typing in number fields
   - Canvas should update smoothly

## Notes

- The fix **only affects Autofit mode** - Integration mode is unchanged
- Axis zoom/pan is preserved across refits (no jumping around)
- The refit uses the same `GetDiagEx()` call that marker dragging already triggered
- No changes to the HTML frontend were needed - this was purely a backend issue

## Files Modified

- `WebShapeIt.cxx` - Updated 6 message handlers

## Compatibility

- ✅ No API changes
- ✅ No settings file format changes  
- ✅ No breaking changes to existing functionality
- ✅ Works with all existing ROOT/ShapeIt code

---

**Date:** 2026-07-14  
**Author:** Xcode Assistant  
**Issue:** Autofit not updating immediately when side panel values change  
**Status:** Fixed
