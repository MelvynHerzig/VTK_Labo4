# Labo 4 from VTK.
#
# File goal: This file aims to provide functions to generate the terrain and the glider path.
#             In addition, it provides useful utility functions like the possibility
#             to convert a set of points from RT90 to WGS84 coordinates or a geographical coordinate
#             into a vtkPoint.
#
# Authors: Forestier Quentin & Herzig Melvyn
#
# Date: 14.06.2022

import math
import numpy as np
from pyproj import Transformer
import constants_Forestier_Herzig as constants
import vtk


def to_vtk_point(altitude, latitude, longitude):
    """
    Convert a geographical coordinate into a vtk point.
    :param altitude: Point's altitude
    :param latitude: Point's latitude
    :param longitude: Point's longitude
    :return: Return the corresponding vtk point.
    """

    coordinate_transform = vtk.vtkTransform()
    coordinate_transform.RotateY(longitude)
    coordinate_transform.RotateX(-latitude)
    coordinate_transform.Translate(0, 0, constants.EARTH_RADIUS + altitude)

    return coordinate_transform.TransformPoint(0, 0, 0)


def convert_rt90_wgs84(x, y):
    """
    Convert an RT90 (swedish) coordinate into WGS84 coordinate.
    :param x: First part of the coordinate
    :param y: Second part of the coordinate
    :return: Return the resulting conversion.
    """
    return convert_rt90_wgs84.transformer.transform(y, x)


# Static attribute of convert_rt90_wgs84 function to make coordinates conversion.
convert_rt90_wgs84.transformer = Transformer.from_crs('epsg:3021', 'epsg:4326')


def quadrilateral_interpolation(x, y, a, b):
    """
    Get a quadrilateral interpolation: converts physical (x,y) to logical (l,m).
    Source: https://www.particleincell.com/2012/quad-interpolation/
    :param x: Position x
    :param y: Position y
    :param a: Alphas vector
    :param b: Betas vector
    :return: Logical coordinates (l, m)
    """

    # quadratic equation coeffs, aa*mm^2+bb*m+cc=0
    aa = a[3] * b[2] - a[2] * b[3]
    bb = a[3] * b[0] - a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + x * b[3] - y * a[3]
    cc = a[1] * b[0] - a[0] * b[1] + x * b[1] - y * a[1]

    # compute m = (-b+sqrt(b^2-4ac))/(2a)
    det = math.sqrt(bb ** 2 - 4 * aa * cc)
    m = (-bb - det) / (2 * aa)

    # compute l
    l = (x - a[0] - a[2] * m) / (a[1] + a[3] * m)

    return l, m


def extract_terrain_data(lat_start, lat_end, lon_start, lon_end, nb_rows, nb_cols, bil_file):
    """
    Extract the altitudes from the bil file in argument. Reshape it into 2d in regard to nb_rows and nb_cols. Finally,
    produces the latitudes and longitudes as vector that correspond to the latitudes and longitudes associated to
    2d array of altitudes.
    :param lat_start: Starting latitude.
    :param lat_end: Ending latitude.
    :param lon_start: Starting longitudes.
    :param lon_end: Ending longitudes.
    :param nb_rows: Row count when resizing the elevations array into 2d.
    :param nb_cols: Column count when resizing the elevations array into 2d.
    :param bil_file: bil file name.
    :return: Returns the 2d array of altitudes, the vector of latitudes and the vector of longitudes.
    """

    # Gets the elevations from the bil file as a big array from 6000 * 6000 altitudes.
    altitudes = np.fromfile(bil_file, dtype=np.int16).reshape((nb_rows, nb_cols))

    # Making the latitudes and longitudes vectors in regard to the starting and
    # ending (latitudes and longitudes) with the number of division (nb_rows, nb_cols).
    latitudes_vector = np.linspace(lat_end, lat_start, nb_rows)
    longitudes_vector = np.linspace(lon_start, lon_end, nb_cols)

    return altitudes, latitudes_vector, longitudes_vector


def clipping_plane(wgs84_1, wgs84_2):
    """
    Compute a plane that passes through (0,0,0), wgs84_1 and wgs84_2.
    :param wgs84_1: WGS84 coordinate number 1 to pass through.
    :param wgs84_2: WGS84 coordinate number 2 to pass through.
    :return: Return the resulting vtkPlane.
    """

    # Get x,y,z positions of the two additional points.
    p1 = np.array(to_vtk_point(10, wgs84_1[0], wgs84_1[1]))
    p2 = np.array(to_vtk_point(10, wgs84_2[0], wgs84_2[1]))

    # Compute normal for plane orientation
    n = np.cross(p1, p2)

    # Plane creation
    plane = vtk.vtkPlane()
    plane.SetNormal(n)

    return plane


def extract_glider_data(filename):
    """
    Extract data from the glider path file.
    :param filename: Name of the file with glider data.
    :return: An array with the list x y z of each measure.
    """
    with open(filename) as file:
        file.readline()  # First line is not useful for us.

        # Array that stores each position read.
        coordinates = []

        # Read each measure and extract coordinates
        for line in file.readlines():
            values = line.split()
            coordinates.append((int(values[1]), int(values[2]), float(values[3])))

    return coordinates


def make_terrain_actor():
    """
    This function created the terrain actor that is around the lake of Ottsjön.
    :return: Returns the corresponding vtkActor.
    """

    # WSG84 map corners: 0:  bottom left, 1: bottom right, 2: top right, 3: top left
    wsg84_corners = np.array([
        convert_rt90_wgs84(constants.AREA_BL[0], constants.AREA_BL[1]),
        convert_rt90_wgs84(constants.AREA_BR[0], constants.AREA_BR[1]),
        convert_rt90_wgs84(constants.AREA_TR[0], constants.AREA_TR[1]),
        convert_rt90_wgs84(constants.AREA_TL[0], constants.AREA_TL[1]),
    ])

    # Matrix for quadrilateral interpolation
    # https://www.particleincell.com/2012/quad-interpolation/
    interpolation_matrix = np.array([[1, 0, 0, 0], [-1, 1, 0, 0], [-1, 0, 0, 1], [1, -1, 1, -1]])

    # Matrix multiplication to get alphas and betas for quadrilateral interpolation
    interpolation_alphas = interpolation_matrix.dot(wsg84_corners[:, 0])
    interpolation_betas = interpolation_matrix.dot(wsg84_corners[:, 1])

    # Limits to get the bounding box of the map area to display.
    smallest_latitude = wsg84_corners[:, 0].min()  # Bottom most
    biggest_latitude = wsg84_corners[:, 0].max()  # Top most
    smallest_longitude = wsg84_corners[:, 1].min()  # Left most
    biggest_longitude = wsg84_corners[:, 1].max()  # Right most

    altitudes, latitudes, longitudes = extract_terrain_data(
        constants.BIL_LAT_START,
        constants.BIL_LAT_END,
        constants.BIL_LON_START,
        constants.BIL_LON_END,
        constants.BIL_WIDTH,
        constants.BIL_HEIGHT,
        constants.BIL_NAME)

    # Coordinates of the points in the bounding box of the area to display
    xyz_points = vtk.vtkPoints()
    # Altitude of the points in the bounding box of the area to display
    altitude_points = vtk.vtkIntArray()
    # Texture coordinates of the points in the bounding box of the area to display
    texture_coordinates = vtk.vtkFloatArray()
    texture_coordinates.SetNumberOfComponents(2)

    rows = 0  # Number of inserted rows in the bounding box
    cols = 0  # Number of inserted columns in the bounding box

    # For each latitude, is it in the bounding box of the area to display?
    for i, row in enumerate(altitudes):
        if smallest_latitude <= latitudes[i] <= biggest_latitude:
            rows += 1
            # For each longitude, is it in the bounding box of the area to display?
            for j, altitude in enumerate(row):
                if smallest_longitude <= longitudes[j] <= biggest_longitude:

                    # We increase number of columns ony for one row since they have all the
                    # same amount of columns.
                    if rows == 1:
                        cols += 1

                    # At this point, the lat, long pair is inside the bounding box of the area to display,
                    # so we add it to the structured grid.
                    xyz_points.InsertNextPoint(to_vtk_point(altitude, latitudes[i], longitudes[j]))
                    altitude_points.InsertNextValue(altitude)

                    l, m = quadrilateral_interpolation(latitudes[i],
                                                       longitudes[j],
                                                       interpolation_alphas,
                                                       interpolation_betas)

                    texture_coordinates.InsertNextTuple((l, m))

    # Preparing structured grid to display the area
    terrain_grid = vtk.vtkStructuredGrid()
    terrain_grid.SetDimensions(cols, rows, 1)
    terrain_grid.SetPoints(xyz_points)
    terrain_grid.GetPointData().SetScalars(altitude_points)
    terrain_grid.GetPointData().SetTCoords(texture_coordinates)

    # Making union boolean implicit function to cut exact borders of the area to display
    # In order to cut the bounding box to match the area to display, we create planes that
    # go through each edge.
    terrain_implicit_boolean = vtk.vtkImplicitBoolean()
    terrain_implicit_boolean.SetOperationTypeToUnion()
    terrain_implicit_boolean.AddFunction(clipping_plane(wsg84_corners[0], wsg84_corners[1]))
    terrain_implicit_boolean.AddFunction(clipping_plane(wsg84_corners[1], wsg84_corners[2]))
    terrain_implicit_boolean.AddFunction(clipping_plane(wsg84_corners[2], wsg84_corners[3]))
    terrain_implicit_boolean.AddFunction(clipping_plane(wsg84_corners[3], wsg84_corners[0]))

    # Clipped terrain
    terrain_clipped = vtk.vtkClipDataSet()
    terrain_clipped.SetInputData(terrain_grid)
    terrain_clipped.SetClipFunction(terrain_implicit_boolean)
    terrain_clipped.Update()

    terrain_mapper = vtk.vtkDataSetMapper()
    terrain_mapper.SetInputConnection(terrain_clipped.GetOutputPort())
    terrain_mapper.ScalarVisibilityOff()

    # Loading texture
    jpeg_reader = vtk.vtkJPEGReader()
    jpeg_reader.SetFileName(constants.AREA_TEXTURE_NAME)
    terrain_texture = vtk.vtkTexture()
    terrain_texture.SetInputConnection(jpeg_reader.GetOutputPort())

    # Terrain actor
    terrain_actor = vtk.vtkActor()
    terrain_actor.SetMapper(terrain_mapper)
    terrain_actor.SetTexture(terrain_texture)

    return terrain_actor


def make_glider_path_actor():
    """
    This function create the plane gps path actor
    :return: Return the corresponding vtkActor
    """

    # Retrieving coordinates from file.
    coords = extract_glider_data(constants.GLIDER_MEASURES_FILENAME)

    # Points that will be used to form the "tube".
    path_points = vtk.vtkPoints()

    # Array to store the altitude delta between points.
    # That allows to color the "tube" accordingly to the difference with the
    # previous point.
    delta_altitudes = vtk.vtkFloatArray()

    last_elev = coords[0][2]

    for i, (x, y, altitude) in enumerate(coords):
        lat, long = convert_rt90_wgs84(x, y)  # Convert from rt90 to wsg84.

        path_points.InsertNextPoint(to_vtk_point(altitude, lat, long))  # Insert point position.
        delta_altitudes.InsertNextValue(last_elev - altitude)  # Insert point altitude delta

        last_elev = altitude

    # Making lines out of the measures
    path_lines = vtk.vtkLineSource()
    path_lines.SetPoints(path_points)
    path_lines.Update()

    # Setting altitudes.
    path_lines.GetOutput().GetPointData().SetScalars(delta_altitudes)

    # Making the lines thicker with tube filter.
    tube = vtk.vtkTubeFilter()
    tube.SetRadius(25)
    tube.SetInputConnection(path_lines.GetOutputPort())

    # Mapping and actor creation.
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())
    mapper.SetScalarRange((-5, 5))

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    return actor


def make_altitude_text_actor():
    """
    Creates the text actor that will display the intersected altitude.
    We set a white background to make it easy to read even on top of the map.
    :return: A text actor.
    """
    altitude_actor = vtk.vtkTextActor()
    altitude_actor.GetTextProperty().SetColor(0, 0, 0)
    altitude_actor.GetTextProperty().SetBackgroundColor(1, 1, 1)
    altitude_actor.GetTextProperty().SetBackgroundOpacity(1)
    altitude_actor.SetInput("")
    altitude_actor.GetTextProperty().SetFontSize(20)
    altitude_actor.SetPosition((40, 40))

    return altitude_actor
