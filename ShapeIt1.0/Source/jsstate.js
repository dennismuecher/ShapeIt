/**
 * ShapeIt 2.0 - Settings State Management
 * 
 * Manages application state including:
 * - Current settings file path
 * - Dirty/clean state tracking
 * - Settings display updates
 */

export class SettingsState {
  constructor() {
    this.currentPath = null;
    this.isDirty = false;
  }

  /**
   * Mark settings as modified (unsaved changes)
   */
  markDirty() {
    this.isDirty = true;
    this.updateDisplay();
  }

  /**
   * Mark settings as clean (just saved or loaded)
   */
  markClean() {
    this.isDirty = false;
    this.updateDisplay();
  }

  /**
   * Set the current settings file path
   */
  setPath(path) {
    this.currentPath = path;
    this.updateDisplay();
  }

  /**
   * Get the current settings file path
   */
  getPath() {
    return this.currentPath;
  }

  /**
   * Check if settings have unsaved changes
   */
  hasDirtyChanges() {
    return this.isDirty;
  }

  /**
   * Update the menubar display with current path and status
   */
  updateDisplay() {
    const pathEl = document.getElementById('settings-path');
    const statusEl = document.getElementById('settings-status');
    
    if (this.currentPath) {
      // Show just the filename, not the full path
      const filename = this.currentPath.split('/').pop();
      pathEl.textContent = filename;
      
      if (this.isDirty) {
        statusEl.textContent = '(modified)';
        statusEl.className = 'modified';
      } else {
        statusEl.textContent = '(saved)';
        statusEl.className = '';
      }
    } else {
      pathEl.textContent = 'No settings file';
      statusEl.textContent = '';
      statusEl.className = '';
    }
  }
}
