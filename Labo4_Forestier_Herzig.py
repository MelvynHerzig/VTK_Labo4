#!/usr/bin/env python
#
# Labo 4 from VTK.
#
# Goal of this file: The goal of this laboratory is to visualize the path of a glider taking off
#                    from the frozen lake of Ottsjön in Sweden. The flight of this glider will also
#                    allow us to visualize the ascending and descending airflow under the influence
#                    of the surrounding mountain relief.
#
# Authors: Forestier Quentin & Herzig Melvyn
#
# Date: 14.06.2022

import vtk
import actors_creation_Forestier_Herzig as actorCreation
import interactor_creation_Forestier_Herzig as customInteractor
import constants_Forestier_Herzig as constants


# --------- Actors ---------
terrain_actor = actorCreation.make_terrain_actor()
glider_path_actor = actorCreation.make_glider_path_actor()
altitude_text_actor = actorCreation.make_altitude_text_actor()

# --------- Render ---------
renderer = vtk.vtkRenderer()
renderer.AddActor(terrain_actor)
renderer.AddActor(glider_path_actor)
renderer.AddActor(altitude_text_actor)
renderer.SetBackground(constants.WHITE)

renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(renderer)
renWin.SetSize(constants.WINDOW_SIZE)

# --------- Interactor ---------
intWin = vtk.vtkRenderWindowInteractor()
intWin.SetRenderWindow(renWin)

style = customInteractor.TerrainInteractorStyle(terrain_actor, altitude_text_actor, renderer)
intWin.SetInteractorStyle(style)

# --------- Display ---------
intWin.Initialize()
renWin.Render()
intWin.Start()
