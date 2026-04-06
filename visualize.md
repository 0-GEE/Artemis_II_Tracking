## Stellarium Directions

Use Stellarium to visualize Orion's position in the sky at your local position on Earth using a custom script.

**USE THIS TO PLAN BEST TIMES AND COORDINATES TO IMAGE** for specific location (even better if setting up custom landscape for your location)

### Shortcut Controls in Stellarium

- **F5** = time
- **F6** = location  
- **F12** = load scripts

### Setup Instructions

1. Run this script to generate needed files
2. Copy 3 files from `Stellarium/` subfolder to Stellarium's Resources/scripts folder
   - Two `.inc` files and one `.ssc` file
   - *Note: For data updates, only the `Orion_Artemis_II_Data.inc` needs updating*
   
   **File destinations:**
   - **MacOS:** `/Applications/Astronomy/Stellarium.app/Contents/Resources/scripts`
     (Right-click on app and select "Show Package Contents" to access)
   - **Windows:** `C:\Program Files\Stellarium\scripts\`

3. In Stellarium, ensure your location matches the location in this script
   - Press **F6** to open location dialog
   
4. Set Date: Change Date/Time to a time during mission
   - Press **F5** to set time/date
   
5. Pause time: Press **K** to toggle

6. *(Optional)* To manually run script or see debug messages, press **F12** to open script window

7. Run Script:
   - Open Configuration Window (press **F2**)
   - Select "Scripts" tab
   - Find "Orion_Artemis_II_Stellarium" in the list
   - Push the "Play" button at the bottom to start script (window will close)
   - Manually step to different times (F5 window)
   - Spacecraft should jump to new position after time updates
   - If spacecraft doesn't jump on large time changes (>1 hr), the script may not be running—restart Stellarium or ensure the time is within the mission timeline

8. When done:
   - Open Configuration Window (F2)
   - Open Scripts Tab
   - Find "Orion_Artemis_II..." and manually stop script

---

## Celestia Directions

To visualize trajectory in 3D, download the [Celestia application](https://celestiaproject.space/)

### Setup Instructions

1. Run this script to generate needed files

2. Copy files from `Celestia/` subdirectory to application's Resources Data directory
   - Two files: `.ssc` + `.xyzv`
   
   **File destinations:**
   - **MacOS:** `/Applications/Astronomy/Celestia.app/Contents/Resources/CelestiaResources/data`
     (Or wherever app is installed; right-click on app and select "Show Package Contents")
   - **Windows:** `C:\Program Files\Celestia\data\`

3. Edit `celestia.cfg` (one time only):
   - Add `"data/orion_artemis_ii_celestia.ssc"` to `SolarSystemCatalogs`
   
   **File locations:**
   - **MacOS:** `/Applications/Astronomy/Celestia.app/Contents/Resources/CelestiaResources/celestia.cfg`
   - **Windows:** `C:\Program Files\Celestia\celestia.cfg`

4. Restart Celestia Application (if open)

5. Open Celestia Application:

   - **Set Time:** Time (menu) → Set Time → Enter a time within the trajectory range
   
   - **Time Controls:**
     - Spacebar = pause/unpause
     - L = faster time speed
     - K = slower time speed
     - J = reverse time (toggle)
   
   - **Display Options:**
     - **(Mac)** Display (menu) → Under Labels and Orbits, ensure moons and spacecraft are **ON**
     - **(Windows)** Render (menu) → View Options → (1) Orbits selected on Right and (2) ensure moons and spacecraft are **checked ON**
   
   - **Follow Orion:**
     - **(Mac)** Location (menu) → Browser → Select Planets → Earth → Spacecraft → Orion → Push Select, Center, Follow buttons
     - **(Windows)** Navigation (menu) → Solar System Browser → Earth → Orion → Push Center and Go To buttons
   
   - **Mouse Controls:**
     - Orbit View (3D): Right-click drag
     - Zoom: Mouse wheel
     - Pan: Left-click drag
   
   - **Enhanced Earth Imagery:** See [celestiamotherlode.net](http://celestiamotherlode.net/) for improved earth imagery and other textures

6. **NOTE:** Each time the script is re-run with new data, the files need to be re-copied into the Celestia Data directories