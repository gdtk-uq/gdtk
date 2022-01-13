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
nozzlectrlptpatchgridvtk = LegacyVTKReader(FileNames=[vtkFile])

# show data in view
nozzlectrlptpatchgridvtkDisplay = Show(nozzlectrlptpatchgridvtk, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
nozzlectrlptpatchgridvtkDisplay.Representation = 'Surface'
nozzlectrlptpatchgridvtkDisplay.ColorArrayName = [None, '']
nozzlectrlptpatchgridvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
nozzlectrlptpatchgridvtkDisplay.SelectOrientationVectors = 'None'
nozzlectrlptpatchgridvtkDisplay.ScaleFactor = 0.02286000028252602
nozzlectrlptpatchgridvtkDisplay.SelectScaleArray = 'None'
nozzlectrlptpatchgridvtkDisplay.GlyphType = 'Arrow'
nozzlectrlptpatchgridvtkDisplay.GlyphTableIndexArray = 'None'
nozzlectrlptpatchgridvtkDisplay.GaussianRadius = 0.0011430000141263007
nozzlectrlptpatchgridvtkDisplay.SetScaleArray = [None, '']
nozzlectrlptpatchgridvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
nozzlectrlptpatchgridvtkDisplay.OpacityArray = [None, '']
nozzlectrlptpatchgridvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
nozzlectrlptpatchgridvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
nozzlectrlptpatchgridvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
nozzlectrlptpatchgridvtkDisplay.ScalarOpacityUnitDistance = 0.019351093878880326

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.038100000470876694, 0.034650806337594986, 10000.0]
renderView1.CameraFocalPoint = [0.038100000470876694, 0.034650806337594986, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change representation type
nozzlectrlptpatchgridvtkDisplay.SetRepresentationType('Outline')

# change solid color
nozzlectrlptpatchgridvtkDisplay.AmbientColor = [0.0, 0.3333333333333333, 1.0]
nozzlectrlptpatchgridvtkDisplay.DiffuseColor = [0.0, 0.3333333333333333, 1.0]

# Properties modified on nozzlectrlptpatchgridvtkDisplay
nozzlectrlptpatchgridvtkDisplay.LineWidth = 4.0

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

# update the view to ensure updated data information
renderView1.Update()

# change representation type
nozzlecoonspatchgridvtkDisplay.SetRepresentationType('Wireframe')

# change solid color
nozzlecoonspatchgridvtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
nozzlecoonspatchgridvtkDisplay.DiffuseColor = [0.0, 0.0, 0.0]

# create a new 'Legacy VTK Reader'
nozzlectrlptsvtk = LegacyVTKReader(FileNames=[vtkPtsFile])

# show data in view
nozzlectrlptsvtkDisplay = Show(nozzlectrlptsvtk, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
nozzlectrlptsvtkDisplay.Representation = 'Surface'
nozzlectrlptsvtkDisplay.ColorArrayName = [None, '']
nozzlectrlptsvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
nozzlectrlptsvtkDisplay.SelectOrientationVectors = 'None'
nozzlectrlptsvtkDisplay.ScaleFactor = 0.02286000028252602
nozzlectrlptsvtkDisplay.SelectScaleArray = 'None'
nozzlectrlptsvtkDisplay.GlyphType = 'Arrow'
nozzlectrlptsvtkDisplay.GlyphTableIndexArray = 'None'
nozzlectrlptsvtkDisplay.GaussianRadius = 0.0011430000141263007
nozzlectrlptsvtkDisplay.SetScaleArray = [None, '']
nozzlectrlptsvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
nozzlectrlptsvtkDisplay.OpacityArray = [None, '']
nozzlectrlptsvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
nozzlectrlptsvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
nozzlectrlptsvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
nozzlectrlptsvtkDisplay.ScalarOpacityUnitDistance = 0.07962458777921622

# update the view to ensure updated data information
renderView1.Update()

# change representation type
nozzlectrlptsvtkDisplay.SetRepresentationType('Points')

# change solid color
nozzlectrlptsvtkDisplay.AmbientColor = [1.0, 0.0, 0.0]
nozzlectrlptsvtkDisplay.DiffuseColor = [1.0, 0.0, 0.0]

# Properties modified on nozzlectrlptsvtkDisplay
nozzlectrlptsvtkDisplay.PointSize = 8.0

# create a new 'Legacy VTK Reader'
nozzlectrlptsvtk_1 = LegacyVTKReader(FileNames=[vtkPtsFile])

# show data in view
nozzlectrlptsvtk_1Display = Show(nozzlectrlptsvtk_1, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
nozzlectrlptsvtk_1Display.Representation = 'Surface'
nozzlectrlptsvtk_1Display.ColorArrayName = [None, '']
nozzlectrlptsvtk_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
nozzlectrlptsvtk_1Display.SelectOrientationVectors = 'None'
nozzlectrlptsvtk_1Display.ScaleFactor = 0.02286000028252602
nozzlectrlptsvtk_1Display.SelectScaleArray = 'None'
nozzlectrlptsvtk_1Display.GlyphType = 'Arrow'
nozzlectrlptsvtk_1Display.GlyphTableIndexArray = 'None'
nozzlectrlptsvtk_1Display.GaussianRadius = 0.0011430000141263007
nozzlectrlptsvtk_1Display.SetScaleArray = [None, '']
nozzlectrlptsvtk_1Display.ScaleTransferFunction = 'PiecewiseFunction'
nozzlectrlptsvtk_1Display.OpacityArray = [None, '']
nozzlectrlptsvtk_1Display.OpacityTransferFunction = 'PiecewiseFunction'
nozzlectrlptsvtk_1Display.DataAxesGrid = 'GridAxesRepresentation'
nozzlectrlptsvtk_1Display.PolarAxes = 'PolarAxesRepresentation'
nozzlectrlptsvtk_1Display.ScalarOpacityUnitDistance = 0.07962458777921622

# update the view to ensure updated data information
renderView1.Update()

# change representation type
nozzlectrlptsvtk_1Display.SetRepresentationType('Wireframe')

# change solid color
nozzlectrlptsvtk_1Display.AmbientColor = [1.0, 0.0, 0.0]
nozzlectrlptsvtk_1Display.DiffuseColor = [1.0, 0.0, 0.0]

# Properties modified on nozzlectrlptsvtk_1Display
nozzlectrlptsvtk_1Display.LineWidth = 2.0

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.038100000470876694, 0.034650806337594986, 10000.0]
renderView1.CameraFocalPoint = [0.038100000470876694, 0.034650806337594986, 0.0]
renderView1.CameraParallelScale = 0.11943688166882435

# save screenshot
SaveScreenshot(pngFile, renderView1, ImageResolution=[2847, 1171])

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
