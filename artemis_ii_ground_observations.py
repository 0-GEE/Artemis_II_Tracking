# Original script by Shawn Gano (https://www.gano.name/shawn/artemis2/)
#   
# -------- STELLARIUM DIRECTIONS: ----------------
# Use Stellarium to visualize Orion's position in the sky at your local position on Earth using a custom script.
#  USE THIS TO PLAN BEST TIMES AND COORDINATES TO IMAGE for specific LOCATION (even better if setting up custom landscape for your location)
# >> Shortcut Controls in Stellarium:  F5 is time, F6 is location, F12 allows you to load scripts
#  - run this script to generate needed files
#  - Copy: 3 files from Stellarium/ subfolder to Stellarium's Resources/scripts folder (two .inc files, one .ssc file) -- note for data updates, only the "Orion_Artemis_II_Data.inc" needs updating
#      MacOS:  /Applications/Astronomy/ Stellarium.app/Contents/Resources/scripts   (right click on app and select Show Package Contents - to get to this path)
#      Win:   C:\Program Files\Stellarium\scripts\
#  - In Stellarium, be sure your location matches the location in this script (hit F6 to open location dialog)
#  - Set Date: Change Date/Time to a time during mission (press F5 to set time/date)
#  - Pause time [k - toggle]
#  - (optional) To manually run script or see debug messages PRESS F12 to open script window
#  - Run Script:  Open Configuration Window (press F2), select "Scripts" tab, find "Orion_Artemis_II_Setllarium" in the list.  
#        -> push the "Play" button at the bottom to start script (window will close)
#        -> Manually step to different times (F5 window), Spacecraft should jump to new position after time updates
#           If spacecraft doesn't jump (on large time changes like > 1hr -- then script may not be running - double check this if this happens and restart. Or ensure the time with withing mission timeline.
# - When done - open Configuration Window (F2) - open Scripts Tab, find "Orion_Artemis_II..." and manually stop script.
#
# -------- CELESTIA DIRECTIONS: ----------------
#  To visualize trajectory in 3D download Celestia application:   https://celestiaproject.space/ 
#  - Run this script to generate needed files
#  - Copy files from Celestia/ sub directory (two files: .ssc + .xyzv) to application's Resources Data directory
#       MacOS:  /Applications/Astronomy/Celestia.app/Contents/Resources/CelestiaResources/data   (MacOS) [Or wherever app is] (right click on app and select Show Package Contents - to get to this path)
#       Win:  C:\Program Files\Celestia\data\
#  - Edit (one time) celestia.cfg  :  --> ADD  "data/orion_artemis_ii_celestia.ssc" ]  to SolarSystemCatalogs 
#       MacOS:  /Applications/Astronomy/Celestia.app/Contents/Resources/CelestiaResources/celestia.cfg   (MacOs) [Or wherever app is] (right click on app and select Show Package Contents - to get to this path)
#       Win:  C:\Program Files\Celestia\celestia.cfg
#  - Restart Application (if open)
#  - Open Celestia Application:
#       - Set Time:   Time (menu) - > Set Time  and enter a time within the trajectory range
#           - Time Controls:  Spacebar = pause/unpause,  L = faster time speed,  K = slower,  J = reverse time (toggle)
#       - Display Options:  (Mac) Display (menu) -> under Labels and Orbits be SURE  moons and spacecraft are turned ON
#                        Windows: Render (menu) -> View Options  ->   (1) Orbits selected on Right and (2) Labels and Orbits be SURE  moons and spacecraft are checked ON
#       - Follow Orion:  (Mac) Location (menu) ->  Browser, select Planets -> Earth -> Spacecraft -> Orion and push Select, Center, Follow buttons
#                        windows: Navigation (menu) -> Solar System Browser ->  Earth -> Orion - then push center and Go To
#       - Mouse Controls:  Orbit View (3d): right-click drag,  Zoom: mouse wheel, Pan: left-click drag
#       - Improved earth imagery (and other textures) see:  http://celestiamotherlode.net/ 
#  - NOTE: Each time script is re-ran with new data, the files need to be re-copied into the Celestia Data directories
# --------------------------

import os
import numpy as np  # may not be used anymore? 
import math
import matplotlib.pyplot as plt
import tzlocal

from datetime import datetime
from skyfield.api import load, wgs84, utc, Distance, Velocity
from skyfield.positionlib import Geocentric
from dotenv import load_dotenv, find_dotenv

# -----------------------------------
# ----- Customizable Data: ---------
# -----------------------------------

load_dotenv(find_dotenv())

# Orion ephemeris OEM file to load:
orionEphmerisFileName = os.getenv("EPHEMERIS_FILE", "Artemis_II_OEM_2026_04_04_to_EI.asc")

# Earth Obersrvation Location (Ground Site)
observationLocationName = os.getenv("LOCATION_NAME", "Unnamed Location")
observation_lat         = float(os.getenv("LATITUDE"))
observation_lon         = float(os.getenv("LONGITUDE"))
observationElevation_m  = float(os.getenv("ELEVATION")) # metres

# Observation equipment specs (telescope and camera)
observationTelescopeFocalLengh_mm      = float(os.getenv("FOCAL_LENGTH", 1422))  # focal length of telescope in mm
observationTelescopeCameraPixelSize_um = float(os.getenv("PIXEL_SIZE", 3.76))  # pixel size of camera in micro meters (um) [1e-6 meters], assuming square pixels


# Output Table Display Options for showing date/time
showUTCtime = bool(os.getenv("SHOW_UTC_TIME", False))
showJD = bool(os.getenv("SHOW_JD", False))
showMET = bool(os.getenv("SHOW_MET", True))

# Output filename for Celestia .xyzv data:
celestiaOutputTrajFileName = os.getenv("CELESTIA_FILENAME", "Celestia/Orion_Artemis_II_Traj_Celestia.xyzv") 
# Output filename for Celestia .ssc Spacecraft file:
celestiaOutputSpaceCraftFileName = os.getenv("CELESTIA_SPACECRAFT_NAME", "Celestia/orion_artemis_ii_celestia.ssc")

# Output filename for Stellarium data:
stellariumDataOutputFileName = os.getenv("SETALLARIUM_FILENAME", "Stellarium/Orion_Artemis_II_Data.inc")


# ----------------------------------------------------
# ----- Constraints and Assumptions Data:  ----------
# ----------------------------------------------------

# Orion Size (meters) : 19m solar array span, or  5.2m in diameter main spacecraft
orionSize_m = 19.0   

minPixelsToSeeObject = .10  # minumum number of pixels an object needs to cover to be considered "seeable"
                            #   Started at 2 pixels (seemed logical) -- but then noted I had seen starlink and geo satelights... used them to better widen this (as this is very restrictive!)
minGroundAltForViewing_deg = 15 # lowest elevation that can be readily viewed from ground observation site (in all directions)

# Check if sun is above the horizon
# Using -0.8333 degrees for standard sunrise/sunset
#  standard sunsite  when sun is -0.8333 degrees below horizon
#  civil twilight 0 to -6   
#  Nautical twilight -6 to -12
#  astronautical twilight -12 to -18
#  Astronomical night sun is more than 18 degrees below horizon
sunAltAboveWhichIsLight_deg = -12.0 # -12 is start of astronimical twilight


# --------------------------------------
# ----- Constants and Constructors -----
# ---------------------------------------

au2km = 149597870.7 # number of km in 1 AU

#constructors
ts = load.timescale()


def ft_to_m(ft: float) -> float:
    return ft * 0.3048

# --------------------------------------------
# ------ start of main script       ----------
# --------------------------------------------


# calulate Geographic position of ground site (static)
earthGroundLoc = wgs84.latlon(observation_lat, observation_lon, elevation_m = observationElevation_m) 

# calculate the angular size per pixel  (the last bit converts radians to degrees and then to arc seconds and a factor to correct um/mm )
arcSecPerPixel = (observationTelescopeCameraPixelSize_um / observationTelescopeFocalLengh_mm) * 3600 * (180/math.pi) / 1000.0
print(f"Given Telescope with a focal length of {observationTelescopeFocalLengh_mm}mm and Camera with pixels of size {observationTelescopeCameraPixelSize_um}um:" )
print(f"  -- Angular Size per pixel is:  {arcSecPerPixel} arcsec/pixel")



## ======== LOAD OEM (.asc) Trajectory File ==============
#
#   Note this is a very rudimentary file reading method - it just skips a fixed number of lines then reads in the data.
#         if the data format changes this must be updated

launchInfoCommentStr = "" # variable to save comment data about the launch time for this data set
orionEphemArray = [] 

with open(orionEphmerisFileName) as f:

    f.readline() # first line is the CCSDS OEM Format Version 
    launchInfoCommentStr = f.readline()  # second line contains a comment on the launch time for this data set
    
    #skip the next 18 lines  (i= 0 to 17)
    for i in range(18):  
        f.readline()

    # Read in the tabular data:   Date/Time[UTC] X Y Z dX dY dZ  
    #      Time format: 2026-03-07T04:59:02.099
    #      Earth J2000 coordinates and velocity (distances in km, velocity in km/s?)
    # Store in 2D array, each row:   [JulianDate (JD) in UTC , X, Y, Z, dX, dY, dZ]
    for line in f:
        
        #--debug test to be sure we have the correct first line
        #print(line)
        #sys.exit(0)

        lineData = line.split()  # splits line at every space, ignores leading and trailing spaces, and multi spaces are treated as one

        # process the line if it has exactly 7 data items
        if len(lineData) == 7:
            #print(lineData)
            
            # convert the time data (in ISO 8601 format) to Julian Date
            dateTimeFromISO = datetime.fromisoformat(lineData[0])  # parse date and time from ISO 8601 string
            dateTimeFromISO = dateTimeFromISO.replace(tzinfo=utc) # specifically set UTC - as datetime.fromisoformat() doesn't always handle this well

            # Convert the datetime object to a Skyfield Time object
            t = ts.from_datetime(dateTimeFromISO)

            #  t now can be accessed to get different dates including: 
            # t.tt - Terrestrial Time (TT)
            # t.tai - International Atomic Time (TAI), continuous time scale without leap seconds
            # t.ut1 - UT1 is a time scale based on the Earth's actual rotation, UT1 is slightly irregular; UTC is uniform (if accuracy < 1s is needed - can be treated like UTC)
            # careful UT1 is not exactly UTC - but using it here as it is the closesed option
            #print(f"ISO String: {lineData[0]}")
            #print(f"Python Datetime Object: {dateTimeFromISO}")
            #print(f"Julian Date (TT): {t.tt}")
            #print(f"Julian Date (TAI): {t.tai}")  
            #print(f"Julian Date (UT1): {t.ut1}")   

            # add to the ephemeris array, using UTC JD, cartesian J2000 coords, and velocity
            orionEphemArray.append( [t.ut1, float(lineData[1]), float(lineData[2]), float(lineData[3]), float(lineData[4]), float(lineData[5]), float(lineData[6]) ] )

            #sys.exit(0) # used for debug - stop after a single iteration
        # end if lineData.shape ==7

    # end for line in f

# end with open file f

## ======== END LOADING TRAJECTORY FILE ==================
#sys.exit(0) # debug - stop after loading data file


print("*** Ephemeris Loaded: **********************")
#print(f"Number of timesteps: {orionEphemArray.shape[0]}" ) # old numpy version reading in data file (.xyzv)
print(f"Number of timesteps: {len(orionEphemArray)}" )  # new plain python 2D array size
print(f"Total data time span: {orionEphemArray[-1][0]-orionEphemArray[0][0]} days." )


earthCode = 399 # GEOCENTRIC (Center of Earth) - is the designation in JPL data as Earth see: https://rhodesmill.org/skyfield/positions.html
eph = load('de430t.bsp') # load ephemeris for solar system bodies, sun Earth moon, etc.
earthEph, sunEph = eph['earth'], eph['sun'] # save for use to determine if a ground location is in sunlight or not


# define ground with respect to the sun-Earth (used for calculating if ground site is in day/night)
groundObserver = earthEph + earthGroundLoc

groundToOrionDistance_km_array = []
angularSizeArray_arcSec = []
MET_days_array = []
altitudeFromGroundSiteArray_deg = []
isOrionSunlitArray_bool = []
isGroundSiteSunlitArray_bool = []

orionRightAscensionArray_hms = []
orionDeclinationArray_dms= []

# data for Stellarium
orionRightAscensionArray_hrs = []
orionDeclinationArray_deg= []

i = 0

for row in orionEphemArray:
    
    # data for this iteration
    julDate = row[0]

    # create position (Distance [Skyfield object]) J2000 in km   [Note: no need to convert to AU if using this method]
    orion_J2000 = Distance( km=[row[1] ,  row[2] ,  row[3]] )

    # create Velocity [Skyfield object]) km/s   [Note: no need to convert to AU / day  if using this method]
    orion_Velocity = Velocity.km_per_s( [row[4] ,  row[5] ,  row[6]] )

    # create time object from Julian Date
    time = ts.ut1_jd(julDate)  # careful UT1 is not exactly UTC
    #print('Julian Date: ', julDate )

    # The Geocentric J200.0 position  NOTE:  this is technically in the Geocentric Celestial Reference System (GCRS) -- but inputing J2000 (slight error 23 milliarcseconds - may want to fix in the future)
    #  API:   https://rhodesmill.org/skyfield/api-position.html#skyfield.positionlib.Geocentric
    orionGeocentricLoc = Geocentric(
                    position_au = orion_J2000.au,
                    velocity_au_per_d = orion_Velocity.au_per_d,
                    center = earthCode,  # 399 is the Earth
                    t = time)
    

    #Ground Observation point on Earth
    # observer = earthGroundLoc.at(time) # earth ground location at the given time.
    #distanceFromEarthCenter = math.sqrt(observer.position.km[0]*observer.position.km[0]+observer.position.km[1]*observer.position.km[1]+observer.position.km[2]*observer.position.km[2])
    #print("Distance From Earth Center [km] = ", distanceFromEarthCenter)

    # NOTE: not accounting for light travel time when imaging....  

    ## topocentric = Difference in position from satellite and Earth point topocentric (at given time)
    topocentric = orionGeocentricLoc - earthGroundLoc.at(time)

    # calculate the Right Asc. and Declination of the topocentric vector:  (J2000)
    ra, dec, distance = topocentric.radec()  # GRCF?  ICRF ("J2000") --  WILL NEED THIS TO POINT TELESCOPE!
    #print('-- ORION --')
    #print('J2000.0')
    #print(ra._degrees)
    #print(dec._degrees)
    #print(distance.km)
    #print(ra.hms)
    #print(dec.dms)
    orionRightAscensionArray_hms.append( ra.hstr() )
    orionDeclinationArray_dms.append( dec.dstr() )
    orionRightAscensionArray_hrs.append(ra.hours )
    orionDeclinationArray_deg.append( dec.degrees )

    # caluclate the altiude (above horizeon), azimuth and distance from Ground Site to Orion:
    alt, az, dist2 = topocentric.altaz()

    # save data -------
    MET_days_array.append( julDate - orionEphemArray[0][0] )  # mission elapse time in days (just subtract initial JD)
    groundToOrionDistance_km_array.append( distance.km )

    angularSizeRadians = 2.0 * math.atan( orionSize_m / (2.0*distance.km*1000.0) ); 
    angularSizeArcSec = angularSizeRadians * 180.0/math.pi * 60.0 * 60.0 # radians to degrees -> 60 minutes per degree -> 60 seconds per minute

    angularSizeArray_arcSec.append (  angularSizeArcSec )   

    altitudeFromGroundSiteArray_deg.append( alt._degrees)

    
    # Future? calcuate "pixel size" based on distance -- after calcuating angular size? 


    # =======  Check if Orion is sunlit  ===============
    #    NOTE: is_sunlit() only works for satellites of earth (not locations on Earth)
    # https://rhodesmill.org/skyfield/earth-satellites.html#satellite-is-sunlit
    isOrionSunlitArray_bool.append( orionGeocentricLoc.is_sunlit(eph) )
    
    
  
    # ====  calulate if it is night time at Ground Site =======
    
    # Compute sun's position relative to ground site location
    altGroundToSun, azGroundToSun, distanceGroundToSun = groundObserver.at(time).observe(sunEph).apparent().altaz()
    
    # Check if sun is above the horizon
    if altGroundToSun.degrees > sunAltAboveWhichIsLight_deg :   # -0.8333 (normal sunset) , -12.0 cut off is start of astronimical twilight
        #print("Sunlit (Daytime)")
        isGroundSiteSunlitArray_bool.append(True)
    else:
        #print("Shadowed (Nighttime)")
        isGroundSiteSunlitArray_bool.append(False)
  

    i= i + 1 # increment counter

# END for each timestamp in the Orion Ephemeris Array  =======================


# ======================================================================
#      Search Results for Times Orion can be imaged from Ground site
# ======================================================================

# function to convert UTC julian date to a string representation of the time in the local timezone (works across time changes) - also can append UTC time
def convertJulianDateToLocalTimeStr(julianDate, includeUTCtime = False, includeJD = False, MET = 0, includeMET = False):
    # create time object from Julian Date [UTC] / Convert JD to Skyfield Time object (UTC)
    time = ts.ut1_jd(julianDate) # careful UT1 is not exactly UTC
                
    # local timezone   ## NOTE: I tested this across time change boundaries and it still works, always returns current daylight savings or not for whatever date is being considered
    local_timezone = tzlocal.get_localzone()
                
    # Convert the Skyfield Time object to a timezone-aware Python datetime object
    local_datetime = time.astimezone(local_timezone)

    formatted_date = local_datetime.strftime('%Y-%m-%d %H:%M:%S %Z%z')

    if(includeUTCtime):
        formatted_date = formatted_date + "  (" + time.utc_iso()+ ")"

    if(includeJD):
        formatted_date = formatted_date + "  {JD: " +f"{julianDate:0.5f}" + "}"

    if(includeMET):
        formatted_date = formatted_date + "  [MET: " +f"{MET:0.4f}" + "]"

    return formatted_date
#----  End convertJulianDateToLocalTimeStr() function ----------

print("-------------------------------------------------------------")
print(f"  VIEWING SIZE CONSTRAINT  ")
print(f"  Opportunities when Orion is >= {minPixelsToSeeObject} pixels from {observationLocationName}")
print(f"     Focal Length = {observationTelescopeFocalLengh_mm} and Pixel Size = {observationTelescopeCameraPixelSize_um} um")
print("-------------------------------------------------------------")


validTimeSpan = False # bool
lastVal = 0
lastIndex = 0
lastDervativePositive = False
startCount = 0 # keep track of valid time periods


for i in range(0, len(angularSizeArray_arcSec)-1):
    

    if(angularSizeArray_arcSec[i] >= arcSecPerPixel*minPixelsToSeeObject ): # is Orion big enough (covers enough pixels)

        if(validTimeSpan):  #already in valid time space
            # Orion is large enough - search for local maxima
            if( (angularSizeArray_arcSec[i] < lastVal) and lastDervativePositive): 
                # deivative changed from positive to negative - the last point was a local maxima
                julDate = orionEphemArray[i-1][0] # old NumPy orionEphemArray[i-1,0]
                formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i-1], includeMET=showMET)
                print(f"     MAX : "  + formatted_date + f",  RA/Dec: {orionRightAscensionArray_hms[i-1]} / {orionDeclinationArray_dms[i-1]}"  + ",  Size: " + "{:.4f}".format(angularSizeArray_arcSec[i-1]) + " arc-sec"  + f" or {angularSizeArray_arcSec[i-1]/arcSecPerPixel:.2f} pixels")
            #end if - last point was a maximum

            # save current values to last values
            lastVal = angularSizeArray_arcSec[i]
            lastIndex = i

            if(angularSizeArray_arcSec[i] > angularSizeArray_arcSec[i-1]):
                lastDervativePositive = True
            else:
                lastDervativePositive = False

        else: # not valid timeSpan before
                startCount = startCount + 1 #increment number or starts:
                
                # get the corresponding Julian Date (JD)
                julDate = orionEphemArray[i][0] # old NumPy julDate = orionEphemArray[i,0]
                #convert JD to string (using custom function)
                formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)
                
                print(">> START : " + formatted_date + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}"  + ",  Size: " + "{:.4f}".format(angularSizeArray_arcSec[i]) + " arc-sec"  + f" or {angularSizeArray_arcSec[i]/arcSecPerPixel:.2f} pixels")

                lastVal = angularSizeArray_arcSec[i]
                lastIndex = i

                if(i==0):
                    lastDervativePositive = False
                else:  # else see if derivative is currently positive or negative - and set flag accordingly
                    if(angularSizeArray_arcSec[i] > angularSizeArray_arcSec[i-1]):
                        lastDervativePositive = True
                    else:
                        lastDervativePositive = False

                validTimeSpan = True # mark 

    else: # Orion is not big enough 

        if(validTimeSpan): # last time was valid - this is an END
            # get the corresponding Julian Date (JD)
            julDate = orionEphemArray[i][0] # old NumPy julDate = orionEphemArray[i,0]
            #convert JD to string (using custom function)
            formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)

            print("   END   : " + formatted_date + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}" + ",  Size: " + "{:.4f}".format(angularSizeArray_arcSec[i]) + " arc-sec" + f" or {angularSizeArray_arcSec[i]/arcSecPerPixel:.2f} pixels")
            print("-------------------------------------------------------------")
            validTimeSpan = False
    # end else (if angular size of Orion is large enough)

# END for each timestamp in the Orion Ephemeris Array  =======================

# if ended seach in a valid time  -- then close it with last element
if(validTimeSpan):
    validTimeSpan = False
    i = -1 # last element
    # get the corresponding Julian Date (JD)
    julDate = orionEphemArray[-1][0] # old NumPy orionEphemArray[-1,0]
    #convert JD to string (using custom function)
    formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)

    print("   END   : " + formatted_date  + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}" + ",  Size: " + "{:.4f}".format(angularSizeArray_arcSec[i]) + " arc-sec" + f" or {angularSizeArray_arcSec[i]/arcSecPerPixel:.2f} pixels" + " [No More Data]")
    print("-------------------------------------------------------------")
# end if valid time at end of data

# final step print number of valid time periods:
print(f"Total Number of Feasible Time Periods: {startCount} ")
print("") # empty line
    


print("-------------------------------------------------------------")
print(f"  ELEVATION/ALTITUDE CONSTRAINT")
print(f"  Opportunities when Orion is >= {minGroundAltForViewing_deg} degrees above the horizon in {observationLocationName}")
print("-------------------------------------------------------------")

validTimeSpan = False # bool
lastVal = 0
lastIndex = 0
lastDervativePositive = False
startCount = 0 # keep track of valid time periods

for i in range(0, len(altitudeFromGroundSiteArray_deg)-1):
    

    if(altitudeFromGroundSiteArray_deg[i] >= minGroundAltForViewing_deg ): # is Orion above the ground site's horizon

        if(validTimeSpan):  #already in valid time space
            # Orion is above horizon - search for local maxima
            if( (altitudeFromGroundSiteArray_deg[i] < lastVal) and lastDervativePositive): 
                # deivative changed from positive to negative - the last point was a local maxima
                julDate = orionEphemArray[i-1][0] #old NumPy rionEphemArray[i-1,0]
                formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i-1], includeMET=showMET)
                print(f"     MAX : "  + formatted_date + f",  RA/Dec: {orionRightAscensionArray_hms[i-1]} / {orionDeclinationArray_dms[i-1]}"  + ",  Alt: " + "{:.1f}".format(altitudeFromGroundSiteArray_deg[i-1]) + " deg")
            #end if - last point was a maximum

            # save current values to last values
            lastVal = altitudeFromGroundSiteArray_deg[i]
            lastIndex = i

            if(altitudeFromGroundSiteArray_deg[i] > altitudeFromGroundSiteArray_deg[i-1]):
                lastDervativePositive = True
            else:
                lastDervativePositive = False

        else: # not valid timeSpan before
                startCount = startCount + 1 #increment number or starts:
                
                # get the corresponding Julian Date (JD)
                julDate = orionEphemArray[i][0] # old NumPy orionEphemArray[i,0]
                #convert JD to string (using custom function)
                formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)
                
                print(">> START : " + formatted_date + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}"  + ",  Alt: " + "{:.1f}".format(altitudeFromGroundSiteArray_deg[i]) + " deg")

                lastVal = altitudeFromGroundSiteArray_deg[i]
                lastIndex = i

                if(i==0):
                    lastDervativePositive = False
                else:  # else see if derivative is currently positive or negative - and set flag accordingly
                    if(altitudeFromGroundSiteArray_deg[i] > altitudeFromGroundSiteArray_deg[i-1]):
                        lastDervativePositive = True
                    else:
                        lastDervativePositive = False

                validTimeSpan = True # mark 

    else: # Orion is not above the horizon

        if(validTimeSpan): # last time was valid - this is an END
            # get the corresponding Julian Date (JD)
            julDate = orionEphemArray[i][0] # old NumPy orionEphemArray[i,0]
            #convert JD to string (using custom function)
            formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)

            print("   END   : " + formatted_date + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}"  + ",  Alt: " + "{:.1f}".format(altitudeFromGroundSiteArray_deg[i]) + " deg")
            print("-------------------------------------------------------------")
            validTimeSpan = False
    # end else (if angular size of Orion is large enough)

# END for each timestamp in the Orion Ephemeris Array  =======================

# if ended seach in a valid time  -- then close it with last element
if(validTimeSpan):
    validTimeSpan = False
    i = -1 # last element
    # get the corresponding Julian Date (JD)
    julDate = orionEphemArray[-1][0] # old NumPy orionEphemArray[-1,0]
    #convert JD to string (using custom function)
    formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)

    print("   END   : " + formatted_date  + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}"  + ",  Alt: " + "{:.1f}".format(altitudeFromGroundSiteArray_deg[i]) + " deg")
    print("-------------------------------------------------------------")
# end if valid time at end of data

# final step print number of valid time periods:
print(f"Total Number of Feasible Time Periods: {startCount} ")
print("") # empty line







print("-------------------------------------------------------------")
print(f"  LIGHTING CONSTRAINTS ")
print(f"  Opportunities when Orion is Sunlit and {observationLocationName} is dark (sun <= {sunAltAboveWhichIsLight_deg} deg Alt)")
print("-------------------------------------------------------------")

validTimeSpan = False # bool
lastVal = 0
lastIndex = 0
startCount = 0 # keep track of valid time periods

#isOrionSunlitArray_bool
#isGroundSiteSunlitArray_bool

for i in range(0, len(isOrionSunlitArray_bool)-1):
    
    if( isOrionSunlitArray_bool[i]==True and isGroundSiteSunlitArray_bool[i]==False): # is Orion is lit AND ground site is NOT

        if(validTimeSpan):  #already in valid time space
            pass

        else: # not valid timeSpan before
                startCount = startCount + 1 #increment number or starts:
                
                # get the corresponding Julian Date (JD)
                julDate = orionEphemArray[i][0] #old NumPy orionEphemArray[i,0]
                #convert JD to string (using custom function)
                formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)
                
                print(">> START : " + formatted_date + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}")

                validTimeSpan = True # mark 

    else: # Lighting constraints NOT met (one or both)

        if(validTimeSpan): # last time was valid - this is an END
            # get the corresponding Julian Date (JD)
            julDate = orionEphemArray[i][0] # old NumPy orionEphemArray[i,0]
            #convert JD to string (using custom function)
            formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)

            print("   END   : " + formatted_date + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}"  )
            print("-------------------------------------------------------------")
            validTimeSpan = False
    # end else (if angular size of Orion is large enough)

# END for each timestamp in the Orion Ephemeris Array  =======================

# if ended seach in a valid time  -- then close it with last element
if(validTimeSpan):
    validTimeSpan = False
    i = -1 # last element
    # get the corresponding Julian Date (JD)
    julDate = orionEphemArray[-1][0] # old NumPy orionEphemArray[-1,0]
    #convert JD to string (using custom function)
    formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)

    print("   END   : " + formatted_date  + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}" )
    print("-------------------------------------------------------------")
# end if valid time at end of data

# final step print number of valid time periods:
print(f"Total Number of Feasible Time Periods: {startCount} ")
print("") # empty line






print("-------------------------------------------------------------")
print(f"  *********  ALL CONSTRAINTS *********")
print(f"  + Orion is >= {minPixelsToSeeObject} pixels from {observationLocationName}")
print(f"          Focal Length = {observationTelescopeFocalLengh_mm} and Pixel Size = {observationTelescopeCameraPixelSize_um} um")
print(f"  + Orion is >= {minGroundAltForViewing_deg} degrees above the horizon in {observationLocationName}")
print(f"  + Orion is Sunlit ")
print(f"  + {observationLocationName} is dark (sun <= {sunAltAboveWhichIsLight_deg} deg Alt)")
print("-------------------------------------------------------------")

validTimeSpan = False # bool
lastVal = 0
lastIndex = 0
startCount = 0 # keep track of valid time periods

#isOrionSunlitArray_bool
#isGroundSiteSunlitArray_bool

for i in range(0, len(isOrionSunlitArray_bool)-1):
    
    # if all cosntraints met: 
    if( (angularSizeArray_arcSec[i] >= arcSecPerPixel*minPixelsToSeeObject) and (altitudeFromGroundSiteArray_deg[i] >= minGroundAltForViewing_deg) and (isOrionSunlitArray_bool[i]==True) and (isGroundSiteSunlitArray_bool[i]==False) ): # all constraints

        if(validTimeSpan):  #already in valid time space
            pass

        else: # not valid timeSpan before
                startCount = startCount + 1 #increment number or starts:
                
                # get the corresponding Julian Date (JD)
                julDate = orionEphemArray[i][0] # old NumPy orionEphemArray[i,0]
                #convert JD to string (using custom function)
                formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)
                
                print(">> START : " + formatted_date + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}")

                validTimeSpan = True # mark 

    else: # Lighting constraints NOT met (one or both)

        if(validTimeSpan): # last time was valid - this is an END
            # get the corresponding Julian Date (JD)
            julDate = orionEphemArray[i][0] # old NumPy orionEphemArray[i,0]
            #convert JD to string (using custom function)
            formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)

            print("   END   : " + formatted_date + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}"  )
            print("-------------------------------------------------------------")
            validTimeSpan = False
    # end else (if angular size of Orion is large enough)

# END for each timestamp in the Orion Ephemeris Array  =======================

# if ended seach in a valid time  -- then close it with last element
if(validTimeSpan):
    validTimeSpan = False
    i = -1 # last element
    # get the corresponding Julian Date (JD)
    julDate = orionEphemArray[-1][0] # old NumPy  orionEphemArray[-1,0]
    #convert JD to string (using custom function)
    formatted_date = convertJulianDateToLocalTimeStr(julDate, includeUTCtime=showUTCtime, includeJD=showJD, MET = MET_days_array[i], includeMET=showMET)

    print("   END   : " + formatted_date  + f",  RA/Dec: {orionRightAscensionArray_hms[i]} / {orionDeclinationArray_dms[i]}" )
    print("-------------------------------------------------------------")
# end if valid time at end of data

# final step print number of valid time periods:
print(f"Total Number of Feasible Time Periods: {startCount} ")
print("    For Launch: " + launchInfoCommentStr)
print("") # empty line






# ====================================
#      Create Plots
# ====================================

# Create a line plot for ground distance to Orion
plt.figure(1)
plt.plot(MET_days_array, groundToOrionDistance_km_array, label='Ground to Orion', color='blue', linestyle='-')
plt.title('Distance from Orion to ' + observationLocationName)
plt.xlabel('MET (days from data start)')
plt.ylabel('Dist [km]')
plt.legend()
plt.grid(True)
#plt.show()


# Create subplots with shared x axes
fig, axs = plt.subplots(3, 1, sharex=True, sharey=False)

fig.suptitle(f"Observation Site: {observationLocationName}") # , fontsize=16   (add main title)

# Angular Size Plot
axisNum = 0
#plt.figure(2)
axs[axisNum].plot(MET_days_array, angularSizeArray_arcSec, label='Orion Angular Size', color='blue', linestyle='-')
minSizeLabel = "Min Size to Image (" + str(minPixelsToSeeObject) + " pixels wide)"
axs[axisNum].plot([MET_days_array[0], MET_days_array[-1]], [arcSecPerPixel*minPixelsToSeeObject,arcSecPerPixel*minPixelsToSeeObject], label=minSizeLabel, color='red', linestyle='-')
#axs[axisNum].xlabel('MET (days from data start)')
axs[axisNum].set_ylabel('Angular Size [arc-sec]')
axs[axisNum].legend()
axs[axisNum].grid(True)

# Alitude above Horizon plot 
axisNum = 1
#plt.figure(3)
axs[axisNum].plot(MET_days_array, altitudeFromGroundSiteArray_deg, label='Altitude above horizon', color='blue', linestyle='-')
minLocalElevationLabel = "Min Altidue :" + str(minGroundAltForViewing_deg) + " deg"
axs[axisNum].plot([MET_days_array[0], MET_days_array[-1]], [minGroundAltForViewing_deg,minGroundAltForViewing_deg], label=minLocalElevationLabel, color='red', linestyle='-')
#axs[axisNum].title('Altitude Above Horizon of Orion from Ground Site')
#axs[0, axisNum].xlabel('MET (days from data start)')
axs[axisNum].set_ylabel('Altitude [degrees]')
axs[axisNum].legend()
axs[axisNum].grid(True)


# Sun Light-- combined plot
axisNum = 2
axs[axisNum].plot(MET_days_array, isOrionSunlitArray_bool, label='Orion Sunlit', color='blue', linestyle='-')
new_list = [i + 2 for i in isGroundSiteSunlitArray_bool]
axs[axisNum].plot(MET_days_array, new_list, label='Ground Sunlit (>Astro Twilight)', color='green', linestyle='-')
#axs[axisNum].title('Angular Size of Orion from Ground')
#axs[axisNum].xlabel('MET (days from data start)')
axs[axisNum].set_ylabel('1,3 = sun lit')
axs[axisNum].legend()
axs[axisNum].grid(True)


fig.supxlabel('MET (days from data start)')
#fig.supylabel('Y Value')
plt.tight_layout()
plt.subplots_adjust(bottom=0.08) # Default is 0.2
plt.subplots_adjust(left=0.10) #



############################################################################################
###############################   Write Data Files For Stellarium ##########################

# write Stellarium data file (javascript included)
with open(stellariumDataOutputFileName, "w") as file:
    file.write("// Shawn Gano - Stellarium Data (for Artemis II)\n")
    file.write("// File Generated: " + datetime.now(tzlocal.get_localzone()).strftime("%Y-%b-%d %H:%M:%S %Z%z") + "\n" ) 
    file.write("// " + launchInfoCommentStr) 
    file.write("//\n")
    file.write("// Data only good for specific location:  " + observationLocationName + "\n")
    file.write("//\n")
    file.write("//   place in same folder as Orion_ArtemisII.ssc script\n")
    file.write("//\n")
    # ------ RA data ----------
    file.write("// Right Ascention Data [JD, RA (hours)] (Note: must be sorted by x-value)\n")
    file.write("const raLookupTable = [\n")
    for i in range(0, len(orionRightAscensionArray_hrs)):
        #file.write( "[" + str(orionEphemArray[i,0]) + ", " + str(orionRightAscensionArray_hrs[i]) + "],\n" ) #old NumPy
        file.write( "[" + str(orionEphemArray[i][0]) + ", " + str(orionRightAscensionArray_hrs[i]) + "],\n" )
        i = i + 1
    # for i  - right ascension 
    file.write("];\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    # ------ end RA -----------

    # ------ Dec data ----------
    file.write("// Declination Data [JD, Dec (deg)] (Note: must be sorted by x-value)\n")
    file.write("const decLookupTable = [\n")
    for i in range(0, len(orionDeclinationArray_deg)):
        #file.write( "[" + str(orionEphemArray[i,0]) + ", " + str(orionDeclinationArray_deg[i]) + "],\n" ) #old NumPy
        file.write( "[" + str(orionEphemArray[i][0]) + ", " + str(orionDeclinationArray_deg[i]) + "],\n" )
        i = i + 1
    # for i  - declination 
    file.write("];\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    # ------ end Dec -----------


    file.write("\n")
    file.write("\n")

# end with writing to file (Stellarium)
print(" ")
print(">> Stellarium data file: Orion_Artemis_II_Data.inc - successfully written.")

############################################################################################
###############################   Write Data Files For Celestia ##########################

# ====  Write .xyzv data file for Celestia Visualzation:
# format is JD X Y Z dX dY dZ
#   JD is in UTC
#   X,Y,Z are in km
#   dX,dY,dZ are in km/s
with open(celestiaOutputTrajFileName, "w") as file:
    for i in range(0, len(orionEphemArray)):
        file.write( str(orionEphemArray[i][0]) + " " + str(orionEphemArray[i][1]) + " " + str(orionEphemArray[i][2]) + " " + str(orionEphemArray[i][3]) + " " + str(orionEphemArray[i][4]) + " " + str(orionEphemArray[i][5]) + " " + str(orionEphemArray[i][6]) + "\n" )
        i = i + 1
    # for i  - declination 
# end with writing to file (Celestia)
print(f">> Celstia data file: {celestiaOutputTrajFileName} - successfully written.")
print(" ")


# ====  Write .ssc spacecraft file for Celestia Visualzation: [must add correct start and stop time to file]
with open(celestiaOutputSpaceCraftFileName, "w") as file:
    file.write("\n")
    file.write("\"Orion\" \"Sol/Earth\" \n")
    file.write("{\n")
    file.write("	# Orion - Artemis II  Celestia Spacecraft File\n")
    file.write("	#\n")
    file.write("	# " + launchInfoCommentStr ) # note this string has it's own line return
    file.write("	#\n")
    file.write("	#   --> copy this file along with .xyzv file to: /  (then restart Celestia to load)\n")
    file.write("	#       Celestia.app/Contents/Resources/CelestiaResources/data   (MacOS) [Or wherever app is]\n")
    file.write("	#   --> MUST UPDATE  celestia.cfg   --> ADD  \"data/orion_artemis_ii_celestia.ssc\" ]  to SolarSystemCatalogs  \n")
    file.write("	#         /Users/USERNAME/Desktop/Celestia.app/Contents/Resources/CelestiaResources/celestia.cfg   (MacOs) [Or wherever app is]\n")
    file.write("	#\n")
    file.write("	#   --> be sure to have correct capitalization on files!\n")
    file.write("	#\n")
    file.write("    \n")
    file.write("	Class \"spacecraft\"\n")
    file.write("	#Radius 0.009                  # 18 meters across - Solar array wings\n")

    file.write("	Beginning         " + str(orionEphemArray[0][0]) + " \n") # data start Julian Date
    file.write("	Ending            " + str(orionEphemArray[-1][0]) + " \n")  # end Julian Date
    
    file.write("\n")
    file.write("	OrbitFrame { EquatorJ2000 { Center \"Sol/Earth\" } }	\n")
    file.write("	SampledTrajectory { Source \"Orion_Artemis_II_Traj_Celestia.xyzv\" }\n")
    file.write("\n")
    file.write("	Albedo       1.0\n")
    file.write("}\n")
    file.write("\n")

# end with writing to file (Celestia spacecraft)
print(f">> Celstia spacecraft file: {celestiaOutputSpaceCraftFileName} - successfully written.")
print(" ")



###################################################################
#####  Display Plots --- must be last as it "blocks" execution of the rest of the script
###################################################################
plt.show()
