# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

import sys
vtkFile = sys.argv[1]
pngFile = sys.argv[2]

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
ductcoonspatchgridvtk = LegacyVTKReader(FileNames=[vtkFile])

# show data in view
ductcoonspatchgridvtkDisplay = Show(ductcoonspatchgridvtk, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
ductcoonspatchgridvtkDisplay.Representation = 'Surface'
ductcoonspatchgridvtkDisplay.ColorArrayName = [None, '']
ductcoonspatchgridvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ductcoonspatchgridvtkDisplay.SelectOrientationVectors = 'None'
ductcoonspatchgridvtkDisplay.ScaleFactor = 1.2000000000000002
ductcoonspatchgridvtkDisplay.SelectScaleArray = 'None'
ductcoonspatchgridvtkDisplay.GlyphType = 'Arrow'
ductcoonspatchgridvtkDisplay.GlyphTableIndexArray = 'None'
ductcoonspatchgridvtkDisplay.GaussianRadius = 0.06
ductcoonspatchgridvtkDisplay.SetScaleArray = [None, '']
ductcoonspatchgridvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
ductcoonspatchgridvtkDisplay.OpacityArray = [None, '']
ductcoonspatchgridvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ductcoonspatchgridvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
ductcoonspatchgridvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
ductcoonspatchgridvtkDisplay.ScalarOpacityUnitDistance = 1.1885640454213802

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [6.0, 3.509006977081299, 10000.0]
renderView1.CameraFocalPoint = [6.0, 3.509006977081299, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change representation type
ductcoonspatchgridvtkDisplay.SetRepresentationType('Outline')

# change solid color
ductcoonspatchgridvtkDisplay.AmbientColor = [0.0, 0.3333333333333333, 1.0]
ductcoonspatchgridvtkDisplay.DiffuseColor = [0.0, 0.3333333333333333, 1.0]

# Properties modified on ductcoonspatchgridvtkDisplay
ductcoonspatchgridvtkDisplay.LineWidth = 4.0

# create a new 'Legacy VTK Reader'
ductcoonspatchgridvtk_1 = LegacyVTKReader(FileNames=[vtkFile])

# show data in view
ductcoonspatchgridvtk_1Display = Show(ductcoonspatchgridvtk_1, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
ductcoonspatchgridvtk_1Display.Representation = 'Surface'
ductcoonspatchgridvtk_1Display.ColorArrayName = [None, '']
ductcoonspatchgridvtk_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
ductcoonspatchgridvtk_1Display.SelectOrientationVectors = 'None'
ductcoonspatchgridvtk_1Display.ScaleFactor = 1.2000000000000002
ductcoonspatchgridvtk_1Display.SelectScaleArray = 'None'
ductcoonspatchgridvtk_1Display.GlyphType = 'Arrow'
ductcoonspatchgridvtk_1Display.GlyphTableIndexArray = 'None'
ductcoonspatchgridvtk_1Display.GaussianRadius = 0.06
ductcoonspatchgridvtk_1Display.SetScaleArray = [None, '']
ductcoonspatchgridvtk_1Display.ScaleTransferFunction = 'PiecewiseFunction'
ductcoonspatchgridvtk_1Display.OpacityArray = [None, '']
ductcoonspatchgridvtk_1Display.OpacityTransferFunction = 'PiecewiseFunction'
ductcoonspatchgridvtk_1Display.DataAxesGrid = 'GridAxesRepresentation'
ductcoonspatchgridvtk_1Display.PolarAxes = 'PolarAxesRepresentation'
ductcoonspatchgridvtk_1Display.ScalarOpacityUnitDistance = 1.1885640454213802

# update the view to ensure updated data information
renderView1.Update()

# change representation type
ductcoonspatchgridvtk_1Display.SetRepresentationType('Wireframe')

# change solid color
ductcoonspatchgridvtk_1Display.AmbientColor = [0.0, 0.0, 0.0]
ductcoonspatchgridvtk_1Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on ductcoonspatchgridvtk_1Display
ductcoonspatchgridvtk_1Display.LineWidth = 2.0

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [6.0, 3.509006977081299, 10000.0]
renderView1.CameraFocalPoint = [6.0, 3.509006977081299, 0.0]
renderView1.CameraParallelScale = 6.950764703628316

# save screenshot
SaveScreenshot(pngFile, renderView1, ImageResolution=[2847, 1171])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [6.0, 3.509006977081299, 10000.0]
renderView1.CameraFocalPoint = [6.0, 3.509006977081299, 0.0]
renderView1.CameraParallelScale = 6.950764703628316

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
