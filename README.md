# Artemis II Tracking

This project is a fork of [Shawn Gano's script](https://www.gano.name/shawn/artemis2/).
It generates ground-observation planning data and visualization files for Artemis II / Orion using an OEM ephemeris file which can be downloaded from [NASA's website](https://www.nasa.gov/missions/artemis/artemis-2/track-nasas-artemis-ii-mission-in-real-time/) (scroll to the bottom of the article under "Download Artemis II ephemeris data here").

The included script is essentially a carbon copy of the original, except the configuration has been
moved to a separate file so that the data generation can be configured without having to
modify the script itself.

Python 3 is required to use. Download it [Here](https://www.python.org/downloads/) if you do not have it already.


## Setup

1. Clone this repository:

```bash
git clone https://github.com/0-GEE/Artemis_II_Tracking.git
cd Artemis_II_Tracking
```

2. Create a virtual environment:

- Windows:

```powershell
python -m venv env
env\Scripts\activate
```

- macOS/Linux:

```bash
python3 -m venv env
source env/bin/activate
```

3. Install dependencies:

```bash
pip install -r requirements.txt
```


## Configuration

Create a file called `.env` (use that exact name) and copy the contents of `.env.example`
into this file. This file stores the configuration that the script will use at runtime.
**Update the configuration variables as needed.**

`.env` configuration template:

```ini
# (REQUIRED) Path to the ephemeris file. Download the latest file from https://www.nasa.gov/missions/artemis/artemis-2/track-nasas-artemis-ii-mission-in-real-time/
EPHEMERIS_FILE='Path/to/ephemeris/file'



# Name of your observation site
LOCATION_NAME='My Location'

# (REQUIRED) Latitude of your observation site (decimal degrees)
LATITUDE=0.0

# (REQUIRED) Longitude of your observation site (decimal degrees)
LONGITUDE=0.0

# (REQUIRED) Elevation of your observation site (decimal metres)
ELEVATION=0.0



# Focal length of telescope in mm
FOCAL_LENGTH=1422

# Pixel size of camera in micro meters (um) [1e-6 meters], assuming square pixels
PIXEL_SIZE=3.76



# Optional output table display options. Change these if you like.
SHOW_UTC_TIME=False
SHOW_JD=False
SHOW_MET=True




# Output filename for Celestia .xyzv data
CELESTIA_FILENAME='Celestia/Orion_Artemis_II_Traj_Celestia.xyzv'

# Output filename for Celestia .ssc Spacecraft file
CELESTIA_SPACECRAFT_NAME='Celestia/orion_artemis_ii_celestia.ssc'

# Output filename for Stellarium data
SETALLARIUM_FILENAME='Stellarium/Orion_Artemis_II_Data.inc'
```


## Run

Run the script:

```bash
python artemis_ii_ground_observations.py
```

If you want to use the generated data with Celestia and/or Stellarium, see [this guide](/visualize.md)
