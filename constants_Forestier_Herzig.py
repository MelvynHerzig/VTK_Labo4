# Labo 4 from VTK.
#
# Goal of the file: This file lists some constants used by the actors creation or interactor
#                   creation file.
#
# Authors: Forestier Quentin & Herzig Melvyn
#
# Date: 14.06.2022

# --------------------------------------------------------------------------------
# -                                    Bil                                       -
# --------------------------------------------------------------------------------

# File name
BIL_NAME = "EarthEnv-DEM90_N60E010.bil"

# Latitudes exposed
BIL_LAT_START = 60
BIL_LAT_END = 65

# Longitudes exposed
BIL_LON_START = 10
BIL_LON_END = 15

# Sizes (amount of measures)
BIL_WIDTH = 6000
BIL_HEIGHT = 6000

# --------------------------------------------------------------------------------
# -                                Area to display                               -
# --------------------------------------------------------------------------------

# Corners of the area to display (in RT90)
AREA_BL = [1349602, 7005969]
AREA_BR = [1371835, 7006362]
AREA_TR = [1371573, 7022967]
AREA_TL = [1349340, 7022573]

# File name of the texture
AREA_TEXTURE_NAME = "glider_map.jpg"

# --------------------------------------------------------------------------------
# -                                    Glider                                    -
# --------------------------------------------------------------------------------

# File name that contains the glider measures.
GLIDER_MEASURES_FILENAME = "vtkgps.txt"

# --------------------------------------------------------------------------------
# -                                   Earth                                      -
# --------------------------------------------------------------------------------

# Radius of the earth in km to cut the terrain or to compute point position.
EARTH_RADIUS = 6_371_009

# --------------------------------------------------------------------------------
# -                                  Colors                                      -
# --------------------------------------------------------------------------------

# For background color
WHITE = (1, 1, 1)

# For font color
BLACK = (0, 0, 0)

# --------------------------------------------------------------------------------
# -                                  Render                                      -
# --------------------------------------------------------------------------------

# Width and height of the window
WINDOW_SIZE = [800, 800]
