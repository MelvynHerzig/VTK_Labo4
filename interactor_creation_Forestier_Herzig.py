# Labo 4 from VTK.
#
# Goal of the file: This file creates the interactor style that allows to cut the terrain
#                   where the mouse intersects the terrain. It displays the related
#                   altitude curves and updates a text actor.
#
# Authors: Forestier Quentin & Herzig Melvyn
#
# Date: 14.06.2022

import vtk
import constants_Forestier_Herzig as constants


class TerrainInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
    """
    Create a custom interactor that detects the mouse position to cut the terrain (terrain_actor) and displays
    an altitude curve. It also displays the numeric altitude into the altitude_actor.
    """

    def __init__(self, terrain_actor, altitude_actor, renderer):
        # Function to trigger when the mouse moves.
        self.AddObserver("MouseMoveEvent", self.mouse_move_event)

        # Actors to work with.
        self.terrain_actor = terrain_actor  # Terrain tu cut
        self.altitude_text_actor = altitude_actor  # Text to update

        # Sphere that will be sized to earth radius + altitude intercepted in order to cut the terrain.
        self.sphere = vtk.vtkSphere()

        # Cutter that will use a sphere to cut the terrain.
        self.cutter = vtk.vtkCutter()
        self.cutter.SetCutFunction(self.sphere)
        self.cutter.SetInputData(self.terrain_actor.GetMapper().GetInput())

        # The cut data will be transformed into strips.
        self.stripper = vtk.vtkStripper()
        self.stripper.SetInputConnection(self.cutter.GetOutputPort())

        # And finally the altitude strips are joint into tubes.
        self.tube_filter = vtk.vtkTubeFilter()
        self.tube_filter.SetRadius(40)
        self.tube_filter.SetInputConnection(self.stripper.GetOutputPort())

        # Altitude curve actor that wil display the altitude curve made out of the cut
        self.altitude_curve_data_set_mapper = vtk.vtkDataSetMapper()
        self.altitude_curve_data_set_mapper.SetInputConnection(self.tube_filter.GetOutputPort())

        self.altitude_curve_actor = vtk.vtkActor()
        self.altitude_curve_actor.SetMapper(self.altitude_curve_data_set_mapper)

        # Picker
        self.point_picker = vtk.vtkPointPicker()
        self.point_picker.PickFromListOn()
        self.point_picker.AddPickList(self.terrain_actor)

        # PLacing the altitude curve actor in the renderer.
        self.SetDefaultRenderer(renderer)
        renderer.AddActor(self.altitude_curve_actor)

    def mouse_move_event(self, obj, event):

        # Get the selected actor
        mouse_pos = self.GetInteractor().GetEventPosition()
        self.point_picker.Pick(mouse_pos[0], mouse_pos[1], 0, self.GetDefaultRenderer())
        picked_actor = self.point_picker.GetActor()

        # If we have selected an actor (since we have only one actor, it's the terrain)
        if picked_actor:

            # Retrieving altitude.
            altitude = self.point_picker.GetDataSet().GetPointData().GetScalars().GetValue(self.point_picker.GetPointId())

            # Update the sphere radius to cut at the corresponding altitude.
            self.sphere.SetRadius(constants.EARTH_RADIUS + altitude)
            self.cutter.Update()

            # Update the content from the text actor.
            self.altitude_text_actor.SetInput("altitude: " + str(altitude) + "m")

            # In case the curves where hidden, set visibility to on.
            self.altitude_curve_actor.VisibilityOn()

        # No actor intercepted, we must hide the elements.
        else:
            self.altitude_curve_actor.VisibilityOff()
            self.altitude_text_actor.SetInput("")

        self.GetInteractor().Render()
        self.OnMouseMove()
