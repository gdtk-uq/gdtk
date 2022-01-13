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
nozzlecoonspatchgridvtk = LegacyVTKReader(FileNames=[vtkFile])

# show data in view
nozzlecoonspatchgridvtkDisplay = Show(nozzlecoonspatchgridvtk, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
nozzlecoonspatchgridvtkDisplay.Representation = 'Surface'
nozzlecoonspatchgridvtkDisplay.ColorArrayName = [None, '']
nozzlecoonspatchgridvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
nozzlecoonspatchgridvtkDisplay.SelectOrientationVectors = 'None'
nozzlecoonspatchgridvtkDisplay.ScaleFactor = 0.02286000028252602
nozzlecoonspatchgridvtkDisplay.SelectScaleArray = 'None'
nozzlecoonspatchgridvtkDisplay.GlyphType = 'Arrow'
nozzlecoonspatchgridvtkDisplay.GlyphTableIndexArray = 'None'
nozzlecoonspatchgridvtkDisplay.GaussianRadius = 0.0011430000141263007
nozzlecoonspatchgridvtkDisplay.SetScaleArray = [None, '']
nozzlecoonspatchgridvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
nozzlecoonspatchgridvtkDisplay.OpacityArray = [None, '']
nozzlecoonspatchgridvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
nozzlecoonspatchgridvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
nozzlecoonspatchgridvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
nozzlecoonspatchgridvtkDisplay.ScalarOpacityUnitDistance = 0.019351093878880326

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.038100000470876694, 0.034650806337594986, 10000.0]
renderView1.CameraFocalPoint = [0.038100000470876694, 0.034650806337594986, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change representation type
nozzlecoonspatchgridvtkDisplay.SetRepresentationType('Outline')

# change solid color
nozzlecoonspatchgridvtkDisplay.AmbientColor = [0.0, 0.3333333333333333, 1.0]
nozzlecoonspatchgridvtkDisplay.DiffuseColor = [0.0, 0.3333333333333333, 1.0]

# Properties modified on nozzlecoonspatchgridvtkDisplay
nozzlecoonspatchgridvtkDisplay.LineWidth = 4.0

# create a new 'Legacy VTK Reader'
nozzlecoonspatchgridvtk_1 = LegacyVTKReader(FileNames=[vtkFile])

# show data in view
nozzlecoonspatchgridvtk_1Display = Show(nozzlecoonspatchgridvtk_1, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
nozzlecoonspatchgridvtk_1Display.Representation = 'Surface'
nozzlecoonspatchgridvtk_1Display.ColorArrayName = [None, '']
nozzlecoonspatchgridvtk_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
nozzlecoonspatchgridvtk_1Display.SelectOrientationVectors = 'None'
nozzlecoonspatchgridvtk_1Display.ScaleFactor = 0.02286000028252602
nozzlecoonspatchgridvtk_1Display.SelectScaleArray = 'None'
nozzlecoonspatchgridvtk_1Display.GlyphType = 'Arrow'
nozzlecoonspatchgridvtk_1Display.GlyphTableIndexArray = 'None'
nozzlecoonspatchgridvtk_1Display.GaussianRadius = 0.0011430000141263007
nozzlecoonspatchgridvtk_1Display.SetScaleArray = [None, '']
nozzlecoonspatchgridvtk_1Display.ScaleTransferFunction = 'PiecewiseFunction'
nozzlecoonspatchgridvtk_1Display.OpacityArray = [None, '']
nozzlecoonspatchgridvtk_1Display.OpacityTransferFunction = 'PiecewiseFunction'
nozzlecoonspatchgridvtk_1Display.DataAxesGrid = 'GridAxesRepresentation'
nozzlecoonspatchgridvtk_1Display.PolarAxes = 'PolarAxesRepresentation'
nozzlecoonspatchgridvtk_1Display.ScalarOpacityUnitDistance = 0.019351093878880326

# update the view to ensure updated data information
renderView1.Update()

# change representation type
nozzlecoonspatchgridvtk_1Display.SetRepresentationType('Wireframe')

# change solid color
nozzlecoonspatchgridvtk_1Display.AmbientColor = [0.0, 0.0, 0.0]
nozzlecoonspatchgridvtk_1Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on nozzlecoonspatchgridvtk_1Display
nozzlecoonspatchgridvtk_1Display.LineWidth = 2.0

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.038100000470876694, 0.034650806337594986, 10000.0]
renderView1.CameraFocalPoint = [0.038100000470876694, 0.034650806337594986, 0.0]
renderView1.CameraParallelScale = 0.11943688166882435

# save screenshot
SaveScreenshot(pngFile, renderView1, ImageResolution=[2847, 1171])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.038100000470876694, 0.034650806337594986, 10000.0]
renderView1.CameraFocalPoint = [0.038100000470876694, 0.034650806337594986, 0.0]
renderView1.CameraParallelScale = 0.11943688166882435

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
