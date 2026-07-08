# Display Options Implementation Guide

This guide shows you how to add the Display Options controls to WebShapeIt.

## Backend Changes (Already Done ✓)

The WebShapeIt.cxx file has been updated with:

1. **SendSettingsSync** - Now includes displaySingle, displayAvg, and colour settings
2. **DISPLAY_OPTIONS handler** - New message handler for receiving display option changes from the UI

## Frontend Changes (To Add to webshapeit.html)

### 1. Add HTML Controls

Add this section to your webshapeit.html file, preferably after the Options panel but in its own distinct group:

```html
<!-- Display Options Panel -->
<div style="margin: 20px 0; padding: 15px; border: 2px solid #4a90e2; border-radius: 8px; background-color: #f8f9fa;">
  <h3 style="margin-top: 0; color: #2c3e50; border-bottom: 2px solid #4a90e2; padding-bottom: 8px;">
    Display Options
  </h3>
  
  <div style="margin-bottom: 12px;">
    <label style="display: flex; align-items: center; cursor: pointer; padding: 5px 0;">
      <input type="checkbox" id="chk_displaySingle" onchange="sendDisplayOptions()" 
             style="margin-right: 8px; width: 18px; height: 18px; cursor: pointer;">
      <span style="font-weight: 500;">Show Individual Data Points</span>
    </label>
    <div style="margin-left: 26px; font-size: 0.9em; color: #666;">
      Display individual gSF data points from each peak
    </div>
  </div>
  
  <div style="margin-bottom: 12px;">
    <label style="display: flex; align-items: center; cursor: pointer; padding: 5px 0;">
      <input type="checkbox" id="chk_displayAvg" onchange="sendDisplayOptions()" 
             style="margin-right: 8px; width: 18px; height: 18px; cursor: pointer;">
      <span style="font-weight: 500;">Show Average Values</span>
    </label>
    <div style="margin-left: 26px; font-size: 0.9em; color: #666;">
      Display smoothed/averaged gSF data
    </div>
  </div>
  
  <div style="margin-bottom: 8px; padding-top: 10px; border-top: 1px solid #ddd;">
    <label style="display: flex; align-items: center; cursor: pointer; padding: 5px 0;">
      <input type="checkbox" id="chk_colour" onchange="sendDisplayOptions()" 
             style="margin-right: 8px; width: 18px; height: 18px; cursor: pointer;">
      <span style="font-weight: 500;">Use Different Colors for Peaks</span>
    </label>
    <div style="margin-left: 26px; font-size: 0.9em; color: #666;">
      Plot Peak 1 and Peak 2 data in different colors (only affects individual points)
    </div>
  </div>
</div>
```

### 2. Add JavaScript Functions

Add these JavaScript functions to the `<script>` section of webshapeit.html:

```javascript
// Send display options to backend
function sendDisplayOptions() {
    const displaySingle = document.getElementById('chk_displaySingle').checked ? 1 : 0;
    const displayAvg = document.getElementById('chk_displayAvg').checked ? 1 : 0;
    const colour = document.getElementById('chk_colour').checked ? 1 : 0;
    
    const msg = 'DISPLAY_OPTIONS:' + displaySingle + '|' + displayAvg + '|' + colour;
    handle.send(msg);
}

// Update display options UI from settings sync
function updateDisplayOptionsUI(values) {
    // Assumes values array has these indices (update based on your SETTINGS_SYNC order):
    // Index 21: displaySingle
    // Index 22: displayAvg
    // Index 23: colour
    
    if (values.length >= 24) {
        document.getElementById('chk_displaySingle').checked = (values[21] != 0);
        document.getElementById('chk_displayAvg').checked = (values[22] != 0);
        document.getElementById('chk_colour').checked = (values[23] != 0);
    }
}
```

### 3. Update Your Existing SETTINGS_SYNC Handler

Find your existing SETTINGS_SYNC message handler in the JavaScript code and update it to call `updateDisplayOptionsUI`. It should look something like this:

```javascript
// In your message handler (where you handle incoming messages from backend)
if (data.startsWith('SETTINGS_SYNC:')) {
    const values = data.substring(14).split('|').map(parseFloat);
    
    // ... existing code to update other fields ...
    
    // Add this new call:
    updateDisplayOptionsUI(values);
}
```

### 4. Initialize Display Options on Page Load

Add this to your page initialization code (likely in a `window.onload` or similar):

```javascript
// Set default values for display options
document.getElementById('chk_displaySingle').checked = true;  // matches ShapeSetting.h default
document.getElementById('chk_displayAvg').checked = false;    // matches ShapeSetting.h default
document.getElementById('chk_colour').checked = true;         // matches ShapeSetting.h default
```

## How It Works

### Data Flow:

1. **User clicks checkbox** → `sendDisplayOptions()` is called
2. **JavaScript** → Sends `DISPLAY_OPTIONS:1|0|1` message to backend
3. **WebShapeIt.cxx** → Receives message, updates `sett->displaySingle`, `sett->displayAvg`, `sett->colour`
4. **User clicks "ShapeIt!"** → Analysis runs with new display settings
5. **ShapeCollector::getMultGraph()** → Respects `displaySingle` and `displayAvg` flags when building graph
6. **ShapeGSF constructor** → Respects `colour` flag to set marker colors

### Settings Sync Flow:

1. **User loads settings file** → Backend updates sett
2. **Backend** → Calls `SendSettingsSync()`
3. **JavaScript** → Receives `SETTINGS_SYNC:...` with all 24 values
4. **JavaScript** → Calls `updateDisplayOptionsUI()` to update checkboxes
5. **UI** → Now reflects the loaded settings accurately

## Testing

After implementing these changes:

1. **Check defaults**: When you start WebShapeIt, verify the checkboxes show:
   - ✓ Show Individual Data Points (checked)
   - ☐ Show Average Values (unchecked)
   - ✓ Use Different Colors for Peaks (checked)

2. **Test interactions**:
   - Uncheck "Show Individual Data Points" → Run ShapeIt → Should see no individual points
   - Check "Show Average Values" → Run ShapeIt → Should see smoothed average line
   - Uncheck "Use Different Colors" → Run ShapeIt → Both peaks should be same color

3. **Test settings file**:
   - Change display options
   - Save settings
   - Change display options to different values
   - Load settings → Checkboxes should revert to saved values

## Color Codes (Reference)

When `colour = true`:
- Peak 1 (levGraph_1): Color 6 (magenta)
- Peak 2 (levGraph_2): Color 7 (cyan)

When `colour = false`:
- Peak 1 (levGraph_1): Color 6 (magenta)
- Peak 2 (levGraph_2): Color 6 (magenta) - same as Peak 1

## Summary

The implementation is complete in the backend (WebShapeIt.cxx). You just need to:
1. Add the HTML controls to webshapeit.html
2. Add the JavaScript functions
3. Update the SETTINGS_SYNC handler to call updateDisplayOptionsUI
4. Initialize the checkboxes on page load

This will give you full control over displaying individual points, average values, and color differentiation between peaks, with proper sync to/from settings files.
