# ShapeIt 2.0 - Code Organization Guide

This document explains the structure and flow of the webshapeit.html JavaScript code to help you understand and maintain the application.

## 📋 Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Global State](#global-state)
3. [Key Functions Reference](#key-functions-reference)
4. [Message Protocol](#message-protocol)
5. [UI Components](#ui-components)
6. [Data Flow](#data-flow)

---

## Architecture Overview

ShapeIt 2.0 is a **web-based frontend** for nuclear physics analysis, communicating with a **C++ backend** via WebSockets.

```
┌─────────────────────────────────────────────────────────┐
│                    Browser (Frontend)                   │
│                                                         │
│  ┌────────────┐  ┌──────────┐  ┌──────────────────┐  │
│  │   Canvas   │  │ Sidebar  │  │  Settings UI     │  │
│  │  (JSROOT)  │  │  Panels  │  │  & File Browser  │  │
│  └────────────┘  └──────────┘  └──────────────────┘  │
│          │              │                 │            │
│          └──────────────┴─────────────────┘            │
│                         │                               │
│                    WebSocket                            │
│                         │                               │
└─────────────────────────┼───────────────────────────────┘
                          │
                ┌─────────▼──────────┐
                │   C++ Backend      │
                │   (ROOT/ShapeIt)   │
                │                    │
                │  - Matrix analysis │
                │  - GSF calculation │
                │  - Level analysis  │
                └────────────────────┘
```

### Technology Stack

- **JSROOT**: Visualization library for ROOT objects (histograms, plots)
- **WebSockets**: Real-time bidirectional communication with C++ backend
- **ES6 Modules**: Modern JavaScript with import/export
- **ROOT Framework**: Nuclear physics analysis framework (C++ side)

---

## Global State

### Connection State
- `conn_handle` - WebSocket connection to C++ backend
- `startDir` - Initial directory path from server

### Settings State
- `currentSettingsPath` - Full path to loaded settings file
- `settingsDirty` - Boolean tracking unsaved changes

### UI State
- `menuBarActive` - Tracks if any menu dropdown is open
- `fbMode` - Current file browser mode
- `fbCurrentPath` - Current directory in file browser
- `isResizing` - Log panel resize drag state

---

## Key Functions Reference

### Settings Management

#### `markSettingsDirty()`
**Purpose**: Mark settings as modified  
**When called**: Any time a setting is changed by user  
**Effect**: Shows "(modified)" in menubar, enables exit warning

#### `markSettingsClean()`
**Purpose**: Mark settings as saved  
**When called**: After successful save or load  
**Effect**: Shows "(saved)" in menubar, disables exit warning

#### `updateSettingsDisplay()`
**Purpose**: Update menubar display of settings status  
**Updates**: Filename, path, and saved/modified indicator

### Canvas & Visualization

#### `embedCanvas(handle)`
**Purpose**: Initialize JSROOT canvas for displaying plots  
**Parameters**: 
- `handle` - WebSocket connection handle
**Creates**: Live canvas with real-time updates from C++ backend

### UI Management

#### `showTab(name)`
**Purpose**: Switch between sidebar panels  
**Parameters**:
- `name` - Tab name: 'levels', 'excitation', 'projections', 'options'

#### `openFileBrowser(mode, title)`
**Purpose**: Open file selection dialog  
**Parameters**:
- `mode` - What action to take: 'open-matrix', 'open-settings', 'save-settings', 'open-oslo', 'save-settings-and-exit'
- `title` - Dialog title text

### Logging

#### `setStatus(msg)`
**Purpose**: Display message in status bar and log panel  
**Parameters**:
- `msg` - Message string to display  
**Effect**: Adds timestamped entry to log, updates status text

---

## Message Protocol

Communication with C++ backend uses string messages over WebSocket.

### Messages FROM Frontend TO Backend

| Message | Purpose | Example |
|---------|---------|---------|
| `RUN:<params>` | Start ShapeIt analysis | `RUN:685\|900\|1497\|1732\|3500\|6700...` |
| `OPEN:<path>` | Open matrix file | `OPEN:/data/matrix.root` |
| `SAVE_SETTINGS:<path>` | Save settings to file | `SAVE_SETTINGS:/config/mysettings.txt` |
| `LOAD_SETTINGS:<path>` | Load settings from file | `LOAD_SETTINGS:/config/mysettings.txt` |
| `LISTDIR:<path>` | Request directory listing | `LISTDIR:/home/user` |
| `SELECTMATRIX:<index>` | Select which matrix to use | `SELECTMATRIX:1` |
| `MODE:<mode>` | Set analysis mode | `MODE:1` (1=Integration, 2=Autofit) |
| `LEVELENERGIES:<params>` | Update level energies | `LEVELENERGIES:685\|900\|1497\|1732\|0\|0\|0\|0\|0\|0` |
| `EXCITATION:<lo>\|<hi>` | Set excitation energy range | `EXCITATION:3500\|6700` |
| `BGENERGIES:<params>` | Set background regions | `BGENERGIES:260\|360\|700\|800\|850\|950\|1350\|1450` |
| `OPTIONS:<flags>` | Set analysis options | `OPTIONS:1\|0\|1\|1` (interpol\|oslo\|sliding\|bg) |
| `DISPLAY_OPTIONS:<flags>` | Set display options | `DISPLAY_OPTIONS:1\|0\|1` (single\|avg\|colour) |
| `BINSIZE:<lo>\|<hi>\|<var>` | Set bin size | `BINSIZE:400\|400\|0` |
| `NBINSLO:<n>` | Set number of bins (low) | `NBINSLO:1` |
| `NBINSHI:<n>` | Set number of bins (high) | `NBINSHI:1` |
| `INTPARAMS:<params>` | Integration parameters | `INTPARAMS:0\|1\|0\|1` (minCounts\|scaling\|autoScale\|effCorr) |
| `SHOWMATRIX` | Display 2D matrix | `SHOWMATRIX` |
| `SHOWPROJ` | Show diagonal projection | `SHOWPROJ` |
| `SHOWBINPROJ:<bin>` | Show specific bin projection | `SHOWBINPROJ:1` |
| `VERBOSE:<level>` | Set logging verbosity | `VERBOSE:0` (0=Silent, 1=Basic, 2=Debug) |
| `EXIT` | Graceful shutdown | `EXIT` |

### Messages FROM Backend TO Frontend

| Message | Purpose | Handler |
|---------|---------|---------|
| `DIRLIST:<path>\n<entries>` | Directory listing | `renderDirListing()` |
| `DIRERROR:<msg>` | Directory error | `setStatus()` |
| `MATRIXLIST:<index>\n<names>` | Available matrices | Updates matrix dropdown |
| `SETTINGS_SYNC:<values>` | Sync all settings | Updates all UI elements |
| `STARTDIR:<path>` | Initial directory | Stores in `startDir` |
| `NBINS:<params>` | Bin information | Updates bin dropdown |
| `LOG:<msg>` | Single log message | Appends to log panel |
| `LOGBATCH:<msgs>` | Multiple log messages | Batch append to log |
| Plain text | Status message | `setStatus()` |

---

## UI Components

### Menubar Dropdowns

**File Menu**:
- Open matrix - Loads .root file
- Exit - Closes application (checks for unsaved changes)

**Settings Menu**:
- Open settings - Load .txt settings file
- Save settings - Save to current path
- Save settings as - Save to new file
- Set literature file - Load Oslo method reference data

**Display Menu**:
- Verbose level - Controls C++ logging output
  - Silent (0)
  - Level 1 (basic)
  - Level 2 (debug)

### Sidebar Tabs

#### 1. Levels Tab
Configure gamma-ray cascade levels:
- **Level 1**: First gamma-ray energy range
- **Level 2**: Second gamma-ray energy range
- **Doublet Mode**: Add second peak for each level
- **Background**: Two regions per level for background subtraction

#### 2. Excitation & Bins Tab
Control excitation energy binning:
- **Excitation Range**: Energy window for analysis
- **Bin Size**: Width of energy bins (can vary)
- **Nr. of Bins**: How many bins to create
- **Min. Counts**: Threshold for valid bins
- **gSF Scaling**: Scale factor for gamma strength function
- **Efficiency Correction**: Detector efficiency factor

#### 3. Projections Tab
Visualize different views of data:
- **Matrix**: Full 2D coincidence matrix
- **Diagonal Projection**: Sum along diagonal
- **Bin**: Individual excitation energy bin projection

#### 4. Options Tab
Analysis settings:
- **Matrix Selection**: Choose which ROOT matrix to analyze
- **Mode**: Integration vs Autofit
- **Sewing Interpolation**: Smooth between bins
- **Display Expectation**: Show Oslo method prediction
- **Sliding Window**: Use overlapping bins
- **Background Subtraction**: Apply background correction
- **Display Options**:
  - Show individual data points
  - Show average values
  - Different colors for peaks

### Log Panel

- Displays timestamped messages
- Resizable (drag handle between canvas and log)
- Clear button to reset
- Auto-scrolls to newest messages

### File Browser

Modal dialog for file operations:
- Navigate directories
- "Up" button for parent directory
- Manual path entry
- File/directory listing with icons
- Save mode shows filename input

---

## Data Flow

### Application Startup

```
1. Browser loads HTML
2. JSROOT modules imported
3. connectWebWindow() establishes WebSocket
4. onWebsocketOpened() callback fires
   ├─> embedCanvas() creates visualization canvas
   ├─> updateSettingsDisplay() initializes UI
   └─> setStatus('Connected.')
5. Server sends STARTDIR message
6. Server sends MATRIXLIST if file loaded
7. Ready for user interaction
```

### Running Analysis

```
User clicks "ShapeIt!" button
  │
  ├─> Collect all parameters from UI
  │   ├─ Level energies (lvl1-lo, lvl1-hi, lvl2-lo, lvl2-hi)
  │   ├─ Excitation range (exc-lo, exc-hi)
  │   └─ Doublet settings (if enabled)
  │
  ├─> Format as: "RUN:685|900|1497|1732|3500|6700|0|0|0|0|0|0"
  │
  ├─> Send via conn_handle.send()
  │
  └─> Backend processes
        │
        ├─> Performs analysis
        ├─> Sends LOG: messages with progress
        ├─> Sends NBINS: with results
        ├─> Updates canvas with plots
        └─> Sends LOGBATCH: with final summary
```

### Settings File Operations

#### Save Settings
```
User: Settings → Save settings
  │
  ├─ Has current path? 
  │  ├─ YES: conn_handle.send('SAVE_SETTINGS:' + currentSettingsPath)
  │  │       markSettingsClean()
  │  │
  │  └─ NO: openFileBrowser('save-settings', 'Save settings as...')
  │         User selects/enters filename
  │         conn_handle.send('SAVE_SETTINGS:' + fullPath)
  │         currentSettingsPath = fullPath
  │         markSettingsClean()
```

#### Load Settings
```
User: Settings → Open settings
  │
  ├─> openFileBrowser('open-settings', 'Open settings file')
  │
  ├─> User selects file
  │
  ├─> conn_handle.send('LOAD_SETTINGS:' + fullPath)
  │
  ├─> Backend reads file, sends SETTINGS_SYNC: message
  │
  └─> Frontend updates ALL UI elements with synced values
      markSettingsClean()
```

### Exit with Unsaved Changes

```
User clicks File → Exit
  │
  ├─ Are settings dirty?
  │  │
  │  ├─ YES: showExitDialog()
  │  │       User chooses:
  │  │       ├─ Save: Save to current path → doExit()
  │  │       ├─ Save As: Open file browser → doExit()
  │  │       ├─ Discard: doExit()
  │  │       └─ Cancel: hideExitDialog()
  │  │
  │  └─ NO: doExit()
  │         └─> conn_handle.send('EXIT')
  │             window.close()
```

---

## Common Tasks

### Adding a New Setting

1. **Add HTML input element** in appropriate panel
2. **Add to SETTINGS_SYNC handler** to receive value from backend
3. **Add event listener** that calls `markSettingsDirty()` on change
4. **Send to backend** when changed (e.g., `conn_handle.send('NEWSETTING:' + value)`)
5. **Include in RUN message** if needed for analysis

### Adding a New Message Type

1. **Define protocol** (decide message format)
2. **Add handler** in `onWebsocketMsg()` with `if (msg.startsWith('NEWMSG:'))`
3. **Send from appropriate UI element** via `conn_handle.send()`
4. **Update C++ backend** to handle the message

### Debugging Connection Issues

Check browser console for:
- WebSocket connection errors
- Message send/receive logs
- JavaScript errors

Add temporary logging:
```javascript
console.log('Sending:', message);
conn_handle.send(message);
```

---

## Best Practices

### When Modifying Code

1. ✅ **Always call `markSettingsDirty()`** when user changes a setting
2. ✅ **Check `conn_handle`** exists before sending messages
3. ✅ **Use `setStatus()`** for user-visible messages
4. ✅ **Add JSDoc comments** for new functions
5. ✅ **Test with unsaved changes** (exit dialog should work)

### Code Organization

- Keep related functionality in the same section
- Use clear section markers (`// ═══ SECTION NAME ═══`)
- Add inline comments for complex logic
- Use descriptive variable names

---

## Troubleshooting

### Canvas doesn't display
- Check `embedCanvas()` was called successfully
- Verify WebSocket connection is open
- Check browser console for JSROOT errors

### Settings not saving/loading
- Verify file paths are correct
- Check C++ backend has read/write permissions
- Look for DIRERROR messages

### UI not updating after backend message
- Check message format matches protocol
- Verify handler in `onWebsocketMsg()` exists
- Add console.log to see if message received

---

## Further Reading

- [JSROOT Documentation](https://root.cern/js/)
- [ROOT Framework](https://root.cern/)
- [WebSocket API](https://developer.mozilla.org/en-US/docs/Web/API/WebSocket)
- [JSDoc Syntax](https://jsdoc.app/)

---

**Last Updated**: Based on refactoring Phase 1 (CSS extraction complete)
