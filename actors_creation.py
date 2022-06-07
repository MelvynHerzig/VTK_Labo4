#!/usr/bin/env python
#
# Labo 4 from VTK.
#
# File goal: This file aims to provide functions to generate the terrain and the glider path.
#             In addition, it provides useful utility functions like the possibility
#             to convert a set of points from RT90 to WGS84 coordinates or a geographical coordinate
#             into a vtkPoint.
#
# Authors: Forestier Quentin & Herzig Melvyn
#
# Date: 31.05.2022

import math
import numpy as np
from pyproj import Transformer
import vtk


def to_vtkPoint(elevation, latitude, longitude):
    """
    Converts a geographical coordinate into a vtk point.
    :param elevation: Point's elevation
    :param latitude: Point's latitude
    :param longitude: Point's longitudes
    :return: Return the corresponding vtk point.
    """

    cartesian = vtk.vtkTransform()
    cartesian.RotateY(longitude)
    cartesian.RotateX(-latitude)
    cartesian.Translate(0, 0, to_vtkPoint.earth_radius + elevation)

    return cartesian.TransformPoint(0, 0, 0)


to_vtkPoint.earth_radius = 6_371_009


def convert_rt90_wgs84(x, y):
    """
    Converts an RT90 (swedish) coordinate into WGS84 coordinate.
    :param x: First coordinate (x)
    :param y: Second coordinate (y)
    :return: Return the resulting conversion.
    """
    return convert_rt90_wgs84.transformer.transform(y, x)


# Attribute of convert_rt90_wgs84 function to make coordinates conversion.
convert_rt90_wgs84.transformer = Transformer.from_crs('epsg:3021', 'epsg:4326')


def quadrilateral_interpolation(x, y, a, b):
    """
    Gets a quadrilateral interpolation: converts physical (x,y) to logical (l,m).
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


def extract_datas(lat_start, lat_end, lon_start, lon_end, nb_rows, nb_cols, bil_file):
    """
    Extracts the elevations from the bil file in argument. Reshape it into 2d in regard to nb_row and nb_col. Finally,
    produces the latitudes and longitudes as vector that correspond to the latitudes and longitudes associated to
    2d array of elevations.
    :param lat_start: Starting latitude.
    :param lat_end: Ending latitude.
    :param lon_start: Starting longitudes.
    :param lon_end: Ending longitudes.
    :param nb_rows: Row count when resizing the elevations array into 2d.
    :param nb_cols: Column count when resizing the elevations array into 2d.
    :param bil_file: bil file name.
    :return: Returns the 2d array of elevations, the vector of latitudes and the vector of longitudes.
    """

    # Gets the elevations from the bil file as a big array from 6000 * 6000 elevations.
    elevations = np.fromfile(bil_file, dtype=np.int16).reshape((nb_rows, nb_cols))

    # Making the latitudes and longitudes vectors in regard to the starting and
    # ending position with the number of division (nb_rows, nb_cols).
    latitudes_vector = np.linspace(lat_end, lat_start, nb_rows)
    longitudes_vector = np.linspace(lon_start, lon_end, nb_cols)

    return elevations, latitudes_vector, longitudes_vector


def make_terrain_actor():
    '''
    This function create the terrain actor. The are around the lake of Ottsjön.
    :return: Returns the corresponding vtkActor.
    '''

    # WSG84 map corners: 0:  bottom left, 1: bottom right, 2: top right, 3: top left
    WSG84_CORNERS = np.array([
        convert_rt90_wgs84(1349602, 7005969),
        convert_rt90_wgs84(1371835, 7006362),
        convert_rt90_wgs84(1371573, 7022967),
        convert_rt90_wgs84(1349340, 7022573),
    ])

    # Matrix for quadrilateral interpolation
    # https://www.particleincell.com/2012/quad-interpolation/
    INTERPOLATION_MATRIX = np.array([[1, 0, 0, 0], [-1, 1, 0, 0], [-1, 0, 0, 1], [1, -1, 1, -1]])

    # Matrix multiplication to get alphas and betas of quadrilateral interpolation
    INTERPOLATION_ALPHAS = INTERPOLATION_MATRIX.dot(WSG84_CORNERS[:, 0])
    INTERPOLATION_BETAS = INTERPOLATION_MATRIX.dot(WSG84_CORNERS[:, 1])

    # Limits to get the bounding box of the map area.
    smallest_latitude = WSG84_CORNERS[:, 0].min()
    biggest_latitude = WSG84_CORNERS[:, 0].max()
    smallest_longitude = WSG84_CORNERS[:, 1].min()
    biggest_longitude = WSG84_CORNERS[:, 1].max()

    elevations, latitudes, longitudes = extract_datas(60, 65, 10, 15, 6000, 6000, "EarthEnv-DEM90_N60E010.bil")

    xyz_points = vtk.vtkPoints()
    elevation_points = vtk.vtkIntArray()
    texture_coordinates_float_array = vtk.vtkFloatArray()
    texture_coordinates_float_array.SetNumberOfComponents(2)

    rows = 0  # Number of inserted rows
    cols = 0  # Number of inserted columns

    # For each latitude, is it in the bounding box of the map to display?
    for i, row in enumerate(elevations):
        if smallest_latitude <= latitudes[i] <= biggest_latitude:
            rows += 1
            # For each longitude, is it in the bounding box of the map to display?
            for j, altitude in enumerate(row):
                if smallest_longitude <= longitudes[j] <= biggest_longitude:

                    # We increase number of columns ony for one row since they have all the
                    # same amount of columns.
                    if rows == 1:
                        cols += 1

                    # At this point, the lat, long pair is inside the bouding box of the zone to display,
                    # so we add it to the structured grid.
                    xyz_points.InsertNextPoint(to_vtkPoint(altitude, latitudes[i], longitudes[j]))
                    elevation_points.InsertNextValue(altitude)
                    texture_coordinates_float_array.InsertNextTuple(
                        quadrilateral_interpolation(latitudes[i],
                                                    longitudes[j],
                                                    INTERPOLATION_ALPHAS,
                                                    INTERPOLATION_BETAS)
                    )

    terrain_grid = vtk.vtkStructuredGrid()
    terrain_grid.SetDimensions(cols, rows, 1)
    terrain_grid.SetPoints(xyz_points)
    terrain_grid.GetPointData().SetScalars(elevation_points)
    terrain_grid.GetPointData().SetTCoords(texture_coordinates_float_array)

    terrain_mapper = vtk.vtkDataSetMapper()
    terrain_mapper.SetInputData(terrain_grid)
    terrain_mapper.ScalarVisibilityOff()

    jpeg_reader = vtk.vtkJPEGReader()
    jpeg_reader.SetFileName("glider_map.jpg")
    terrain_texture = vtk.vtkTexture()
    terrain_texture.SetInputConnection(jpeg_reader.GetOutputPort())

    terrain_actor = vtk.vtkActor()
    terrain_actor.SetMapper(terrain_mapper)
    terrain_actor.SetTexture(terrain_texture)

    return terrain_actor


# --------- Render ---------
renderer = vtk.vtkRenderer()
renderer.AddActor(make_terrain_actor())

renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(renderer)
renWin.SetSize(800, 800)

# --------- Interactor ---------
intWin = vtk.vtkRenderWindowInteractor()
intWin.SetRenderWindow(renWin)

style = vtk.vtkInteractorStyleTrackballCamera()
intWin.SetInteractorStyle(style)

# --------- Print image ---------
renWin.Render()
intWin.Initialize()
intWin.Start()
