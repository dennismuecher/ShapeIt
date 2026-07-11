# ShapeIt 2.0 - Quick Reference Card

**Quick navigation guide for common tasks in the codebase**

---

## 📁 File Structure

```
ShapeIt/
├── webshapeit.html          ← Main application file (~610 lines, well-organized)
├── webshapeit-styles.css    ← All CSS styling (~240 lines)
├── CODE_GUIDE.md            ← Complete documentation with diagrams
├── REFACTORING_PROGRESS.md  ← History of improvements
└── js/
    └── state.js             ← (Unused - kept for reference)
```

---

## 🔍 Finding Things in webshapeit.html

The JavaScript code is organized with clear section markers. Search for:

| Looking for... | Search for this |
|----------------|----------------|
| **Global variables** | `GLOBAL STATE` |
| **Settings save/load** | `SETTINGS STATE MANAGEMENT` |
| **Canvas/plotting** | `CANVAS SETUP` |
| **Tab switching** | `SIDEBAR TAB SWITCHING` |
| **Menu dropdowns** | `MENUBAR DROPDOWN` |
| **File browser** | `FILE BROWSER` |
| **Server communication** | `WEBSOCKET COMMUNICATION` |
| **Button click handlers** | Search for the button ID |
| **Settings inputs** | Search for input ID |

---

## 🎯 Common Modifications

### Adding a New Setting Input

**Example: Adding a new checkbox called "use-advanced-mode"**

1. **Add HTML** (in appropriate sidebar panel):
```html
<label>
  <input type="checkbox" id="opt-advanced-mode"> Advanced Mode
</label>
```

2. **Add event listener** (in INPUT EVENT HANDLERS section):
```javascript
document.getElementById('opt-advanced-mode').addEventListener('change', () => {
  const isAdvanced = document.getElementById('opt-advanced-mode').checked ? 1 : 0;
  if (conn_handle) conn_handle.send('ADVANCEDMODE:' + isAdvanced);
  markSettingsDirty();
});
```

3. **Add to settings sync** (in SETTINGS SYNCHRONIZATION section):
```javascript
// In onWebsocketMsg, SETTINGS_SYNC handler, add:
document.getElementById('opt-advanced-mode').checked = v[30] !== 0;  // Next available index
```

4. **Update C++ backend** to handle `ADVANCEDMODE:` message

### Adding a New Menu Item

**Example: Adding "Export Results" to File menu**

1. **Add HTML** (in File menu dropdown):
```html
<button class="menu-item" id="btn-export">Export results...</button>
```

2. **Add click handler** (near other File menu handlers):
```javascript
document.getElementById('btn-export').addEventListener('click', () => {
  if (conn_handle) conn_handle.send('EXPORT');
  closeAllMenus();
});
```

### Adding a New Sidebar Tab

**Example: Adding an "Advanced" tab**

1. **Add tab button** (in sidebar-tabs):
```html
<button class="tab-btn" data-tab="advanced">Advanced</button>
```

2. **Add panel** (in sidebar-scroll):
```html
<div id="panel-advanced" style="display: none;">
  <h3>Advanced Options</h3>
  <!-- Your content here -->
  <button class="btn-run">ShapeIt!</button>
</div>
```

3. **Register in tabs object**:
```javascript
const tabs = {
  levels: { el: document.getElementById('panel-levels') },
  excitation: { el: document.getElementById('panel-excitation') },
  projections: { el: document.getElementById('panel-projections') },
  options: { el: document.getElementById('panel-options') },
  advanced: { el: document.getElementById('panel-advanced') },  // ADD THIS
};
```

---

## 🔌 Message Protocol Cheat Sheet

### Sending to Backend

```javascript
// Format: conn_handle.send('COMMAND:parameter')

// Examples:
conn_handle.send('OPEN:/path/to/file.root');
conn_handle.send('SAVE_SETTINGS:/path/to/settings.txt');
conn_handle.send('VERBOSE:1');
conn_handle.send('MODE:2');
conn_handle.send('LEVELENERGIES:685|900|1497|1732|0|0|0|0|0|0');
```

### Receiving from Backend

```javascript
// In onWebsocketMsg():
if (msg.startsWith('LOG:')) {
  const logText = msg.slice(4);
  // Handle log message
}
```

---

## 🐛 Debugging Tips

### Check WebSocket Connection

```javascript
// Add to onWebsocketOpened():
console.log('WebSocket connected!', handle);

// Add to onWebsocketMsg():
console.log('Received:', msg);
```

### Check if Message Sent

```javascript
// Before sending:
console.log('Sending:', 'MYCOMMAND:' + value);
conn_handle.send('MYCOMMAND:' + value);
```

### Check if Settings Dirty Flag Works

```javascript
// Add to markSettingsDirty():
console.log('Settings marked dirty, settingsDirty =', settingsDirty);
```

### Test Message Handlers

```javascript
// Manually trigger:
onWebsocketMsg(conn_handle, 'TEST:data');
```

---

## 📝 Code Style Guidelines

### Section Markers

```javascript
  // ═══════════════════════════════════════════════════════════════════════════════
  // SECTION NAME - Brief description
  // ═══════════════════════════════════════════════════════════════════════════════
```

### Function Documentation

```javascript
  /**
   * Brief description of what function does
   * @param {type} paramName - Parameter description
   * @returns {type} Return value description
   */
  function myFunction(paramName) {
    // Implementation
  }
```

### Variable Documentation

```javascript
  /**
   * Description of variable purpose
   * @type {string|null}
   */
  let myVariable = null;
```

### Inline Comments

```javascript
// Short comment for single line

// Longer explanation for complex logic
// can span multiple lines
const result = complexCalculation();
```

---

## ⚠️ Important Gotchas

### Always Check conn_handle

```javascript
// ❌ DON'T DO THIS:
conn_handle.send('COMMAND');

// ✅ DO THIS:
if (conn_handle) conn_handle.send('COMMAND');
```

### Remember markSettingsDirty()

```javascript
// ❌ MISSING:
document.getElementById('my-input').addEventListener('change', (e) => {
  conn_handle.send('NEWVALUE:' + e.target.value);
});

// ✅ CORRECT:
document.getElementById('my-input').addEventListener('change', (e) => {
  conn_handle.send('NEWVALUE:' + e.target.value);
  markSettingsDirty();  // Don't forget this!
});
```

### Settings Sync Array Indices

When adding new settings, use the **next available index** in the SETTINGS_SYNC handler:

```javascript
// Current settings use indices 0-29
// Your new setting should be index 30+
if (v.length >= 31) {
  document.getElementById('my-new-setting').value = v[30];
}
```

---

## 🚀 Performance Tips

### Batch DOM Updates

```javascript
// ❌ SLOW (reflows on each line):
for (const item of items) {
  list.appendChild(createItem(item));
}

// ✅ FAST (single reflow):
const fragment = document.createDocumentFragment();
for (const item of items) {
  fragment.appendChild(createItem(item));
}
list.appendChild(fragment);
```

### Debounce Frequent Events

```javascript
let timeout;
inputElement.addEventListener('input', (e) => {
  clearTimeout(timeout);
  timeout = setTimeout(() => {
    // Send to server after user stops typing
    conn_handle.send('UPDATE:' + e.target.value);
  }, 300);
});
```

---

## 📚 Key Documentation Files

- **CODE_GUIDE.md** - Full documentation with architecture, data flow, and troubleshooting
- **REFACTORING_PROGRESS.md** - History of code improvements
- **This file** - Quick reference for common tasks

---

## 🆘 Getting Help

### In the Code

1. **Search for section markers** (`═══`) to find related code
2. **Read JSDoc comments** above functions
3. **Check inline comments** for complex logic

### In Documentation

1. **CODE_GUIDE.md** - Comprehensive guide
2. **Message Protocol section** - All commands explained
3. **Data Flow diagrams** - Understand how things connect

### External Resources

- [JSROOT docs](https://root.cern/js/) - Canvas visualization
- [ROOT framework](https://root.cern/) - Backend analysis
- [MDN WebSocket](https://developer.mozilla.org/en-US/docs/Web/API/WebSocket) - Communication

---

**Pro tip**: Use your code editor's "Go to Symbol" or "Outline" feature to navigate between sections quickly!
