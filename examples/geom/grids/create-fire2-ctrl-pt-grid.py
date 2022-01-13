# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

import sys
vtkFile = sys.argv[1]
vtkPtsFile = sys.argv[2]
pngFile = sys.argv[3]

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2847, 1171]

# get layout
layout1 = GetLayout()

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# get the material library
materialLibrary1 = GetMaterialLibrary()

# create a new 'Legacy VTK Reader'
fire2ctrlptpatchgridvtk = LegacyVTKReader(FileNames=[vtkFile])

# show data in view
fire2ctrlptpatchgridvtkDisplay = Show(fire2ctrlptpatchgridvtk, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
fire2ctrlptpatchgridvtkDisplay.Representation = 'Surface'
fire2ctrlptpatchgridvtkDisplay.ColorArrayName = [None, '']
fire2ctrlptpatchgridvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
fire2ctrlptpatchgridvtkDisplay.SelectOrientationVectors = 'None'
fire2ctrlptpatchgridvtkDisplay.ScaleFactor = 0.045406922698020935
fire2ctrlptpatchgridvtkDisplay.SelectScaleArray = 'None'
fire2ctrlptpatchgridvtkDisplay.GlyphType = 'Arrow'
fire2ctrlptpatchgridvtkDisplay.GlyphTableIndexArray = 'None'
fire2ctrlptpatchgridvtkDisplay.GaussianRadius = 0.002270346134901047
fire2ctrlptpatchgridvtkDisplay.SetScaleArray = [None, '']
fire2ctrlptpatchgridvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
fire2ctrlptpatchgridvtkDisplay.OpacityArray = [None, '']
fire2ctrlptpatchgridvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
fire2ctrlptpatchgridvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
fire2ctrlptpatchgridvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
fire2ctrlptpatchgridvtkDisplay.ScalarOpacityUnitDistance = 0.05481861671452429

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.03642173111438751, 0.22703461349010468, 10000.0]
renderView1.CameraFocalPoint = [0.03642173111438751, 0.22703461349010468, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change representation type
fire2ctrlptpatchgridvtkDisplay.SetRepresentationType('Outline')

# change solid color
fire2ctrlptpatchgridvtkDisplay.AmbientColor = [0.0, 0.3333333333333333, 1.0]
fire2ctrlptpatchgridvtkDisplay.DiffuseColor = [0.0, 0.3333333333333333, 1.0]

# Properties modified on fire2ctrlptpatchgridvtkDisplay
fire2ctrlptpatchgridvtkDisplay.LineWidth = 4.0

# create a new 'Legacy VTK Reader'
fire2ctrlptpatchgridvtk_1 = LegacyVTKReader(FileNames=[vtkFile])

# show data in view
fire2ctrlptpatchgridvtk_1Display = Show(fire2ctrlptpatchgridvtk_1, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
fire2ctrlptpatchgridvtk_1Display.Representation = 'Surface'
fire2ctrlptpatchgridvtk_1Display.ColorArrayName = [None, '']
fire2ctrlptpatchgridvtk_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
fire2ctrlptpatchgridvtk_1Display.SelectOrientationVectors = 'None'
fire2ctrlptpatchgridvtk_1Display.ScaleFactor = 0.045406922698020935
fire2ctrlptpatchgridvtk_1Display.SelectScaleArray = 'None'
fire2ctrlptpatchgridvtk_1Display.GlyphType = 'Arrow'
fire2ctrlptpatchgridvtk_1Display.GlyphTableIndexArray = 'None'
fire2ctrlptpatchgridvtk_1Display.GaussianRadius = 0.002270346134901047
fire2ctrlptpatchgridvtk_1Display.SetScaleArray = [None, '']
fire2ctrlptpatchgridvtk_1Display.ScaleTransferFunction = 'PiecewiseFunction'
fire2ctrlptpatchgridvtk_1Display.OpacityArray = [None, '']
fire2ctrlptpatchgridvtk_1Display.OpacityTransferFunction = 'PiecewiseFunction'
fire2ctrlptpatchgridvtk_1Display.DataAxesGrid = 'GridAxesRepresentation'
fire2ctrlptpatchgridvtk_1Display.PolarAxes = 'PolarAxesRepresentation'
fire2ctrlptpatchgridvtk_1Display.ScalarOpacityUnitDistance = 0.05481861671452429

# update the view to ensure updated data information
renderView1.Update()

# change representation type
fire2ctrlptpatchgridvtk_1Display.SetRepresentationType('Wireframe')

# change solid color
fire2ctrlptpatchgridvtk_1Display.AmbientColor = [0.0, 0.0, 0.0]
fire2ctrlptpatchgridvtk_1Display.DiffuseColor = [0.0, 0.0, 0.0]

# create a new 'Legacy VTK Reader'
fire2ctrlptsvtk = LegacyVTKReader(FileNames=[vtkPtsFile])

# show data in view
fire2ctrlptsvtkDisplay = Show(fire2ctrlptsvtk, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
fire2ctrlptsvtkDisplay.Representation = 'Surface'
fire2ctrlptsvtkDisplay.ColorArrayName = [None, '']
fire2ctrlptsvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
fire2ctrlptsvtkDisplay.SelectOrientationVectors = 'None'
fire2ctrlptsvtkDisplay.ScaleFactor = 0.04540737427842032
fire2ctrlptsvtkDisplay.SelectScaleArray = 'None'
fire2ctrlptsvtkDisplay.GlyphType = 'Arrow'
fire2ctrlptsvtkDisplay.GlyphTableIndexArray = 'None'
fire2ctrlptsvtkDisplay.GaussianRadius = 0.002270368713921016
fire2ctrlptsvtkDisplay.SetScaleArray = [None, '']
fire2ctrlptsvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
fire2ctrlptsvtkDisplay.OpacityArray = [None, '']
fire2ctrlptsvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
fire2ctrlptsvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
fire2ctrlptsvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
fire2ctrlptsvtkDisplay.ScalarOpacityUnitDistance = 0.14642180418341505

# update the view to ensure updated data information
renderView1.Update()

# change representation type
fire2ctrlptsvtkDisplay.SetRepresentationType('Wireframe')

# change solid color
fire2ctrlptsvtkDisplay.AmbientColor = [1.0, 0.0, 0.0]
fire2ctrlptsvtkDisplay.DiffuseColor = [1.0, 0.0, 0.0]

# Properties modified on fire2ctrlptsvtkDisplay
fire2ctrlptsvtkDisplay.LineWidth = 2.0

# Properties modified on fire2ctrlptsvtkDisplay
fire2ctrlptsvtkDisplay.LineWidth = 3.0

# Properties modified on fire2ctrlptsvtkDisplay
fire2ctrlptsvtkDisplay.LineWidth = 2.0

# create a new 'Legacy VTK Reader'
fire2ctrlptsvtk_1 = LegacyVTKReader(FileNames=[vtkPtsFile])

# show data in view
fire2ctrlptsvtk_1Display = Show(fire2ctrlptsvtk_1, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
fire2ctrlptsvtk_1Display.Representation = 'Surface'
fire2ctrlptsvtk_1Display.ColorArrayName = [None, '']
fire2ctrlptsvtk_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
fire2ctrlptsvtk_1Display.SelectOrientationVectors = 'None'
fire2ctrlptsvtk_1Display.ScaleFactor = 0.04540737427842032
fire2ctrlptsvtk_1Display.SelectScaleArray = 'None'
fire2ctrlptsvtk_1Display.GlyphType = 'Arrow'
fire2ctrlptsvtk_1Display.GlyphTableIndexArray = 'None'
fire2ctrlptsvtk_1Display.GaussianRadius = 0.002270368713921016
fire2ctrlptsvtk_1Display.SetScaleArray = [None, '']
fire2ctrlptsvtk_1Display.ScaleTransferFunction = 'PiecewiseFunction'
fire2ctrlptsvtk_1Display.OpacityArray = [None, '']
fire2ctrlptsvtk_1Display.OpacityTransferFunction = 'PiecewiseFunction'
fire2ctrlptsvtk_1Display.DataAxesGrid = 'GridAxesRepresentation'
fire2ctrlptsvtk_1Display.PolarAxes = 'PolarAxesRepresentation'
fire2ctrlptsvtk_1Display.ScalarOpacityUnitDistance = 0.14642180418341505

# update the view to ensure updated data information
renderView1.Update()

# change representation type
fire2ctrlptsvtk_1Display.SetRepresentationType('Points')

# Properties modified on fire2ctrlptsvtk_1Display
fire2ctrlptsvtk_1Display.PointSize = 8.0

# change solid color
fire2ctrlptsvtk_1Display.AmbientColor = [1.0, 0.0, 0.0]
fire2ctrlptsvtk_1Display.DiffuseColor = [1.0, 0.0, 0.0]

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.03642173111438751, 0.22703461349010468, 10000.0]
renderView1.CameraFocalPoint = [0.03642173111438751, 0.22703461349010468, 0.0]
renderView1.CameraParallelScale = 0.24173495193661831

# save screenshot
SaveScreenshot(pngFile, renderView1, ImageResolution=[2847, 1171])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.03642173111438751, 0.22703461349010468, 10000.0]
renderView1.CameraFocalPoint = [0.03642173111438751, 0.22703461349010468, 0.0]
renderView1.CameraParallelScale = 0.24173495193661831

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
