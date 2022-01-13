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
ductctrlptpatchgridvtk = LegacyVTKReader(FileNames=[vtkFile])

# show data in view
ductctrlptpatchgridvtkDisplay = Show(ductctrlptpatchgridvtk, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
ductctrlptpatchgridvtkDisplay.Representation = 'Surface'
ductctrlptpatchgridvtkDisplay.ColorArrayName = [None, '']
ductctrlptpatchgridvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ductctrlptpatchgridvtkDisplay.SelectOrientationVectors = 'None'
ductctrlptpatchgridvtkDisplay.ScaleFactor = 1.2000000000000002
ductctrlptpatchgridvtkDisplay.SelectScaleArray = 'None'
ductctrlptpatchgridvtkDisplay.GlyphType = 'Arrow'
ductctrlptpatchgridvtkDisplay.GlyphTableIndexArray = 'None'
ductctrlptpatchgridvtkDisplay.GaussianRadius = 0.06
ductctrlptpatchgridvtkDisplay.SetScaleArray = [None, '']
ductctrlptpatchgridvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
ductctrlptpatchgridvtkDisplay.OpacityArray = [None, '']
ductctrlptpatchgridvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ductctrlptpatchgridvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
ductctrlptpatchgridvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
ductctrlptpatchgridvtkDisplay.ScalarOpacityUnitDistance = 1.1885640454213802

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [6.0, 3.509006977081299, 10000.0]
renderView1.CameraFocalPoint = [6.0, 3.509006977081299, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change representation type
ductctrlptpatchgridvtkDisplay.SetRepresentationType('Outline')

# change solid color
ductctrlptpatchgridvtkDisplay.AmbientColor = [0.0, 0.3333333333333333, 1.0]
ductctrlptpatchgridvtkDisplay.DiffuseColor = [0.0, 0.3333333333333333, 1.0]

# Properties modified on ductctrlptpatchgridvtkDisplay
ductctrlptpatchgridvtkDisplay.LineWidth = 4.0

# create a new 'Legacy VTK Reader'
ductctrlptpatchgridvtk_1 = LegacyVTKReader(FileNames=[vtkFile])

# show data in view
ductctrlptpatchgridvtk_1Display = Show(ductctrlptpatchgridvtk_1, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
ductctrlptpatchgridvtk_1Display.Representation = 'Surface'
ductctrlptpatchgridvtk_1Display.ColorArrayName = [None, '']
ductctrlptpatchgridvtk_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
ductctrlptpatchgridvtk_1Display.SelectOrientationVectors = 'None'
ductctrlptpatchgridvtk_1Display.ScaleFactor = 1.2000000000000002
ductctrlptpatchgridvtk_1Display.SelectScaleArray = 'None'
ductctrlptpatchgridvtk_1Display.GlyphType = 'Arrow'
ductctrlptpatchgridvtk_1Display.GlyphTableIndexArray = 'None'
ductctrlptpatchgridvtk_1Display.GaussianRadius = 0.06
ductctrlptpatchgridvtk_1Display.SetScaleArray = [None, '']
ductctrlptpatchgridvtk_1Display.ScaleTransferFunction = 'PiecewiseFunction'
ductctrlptpatchgridvtk_1Display.OpacityArray = [None, '']
ductctrlptpatchgridvtk_1Display.OpacityTransferFunction = 'PiecewiseFunction'
ductctrlptpatchgridvtk_1Display.DataAxesGrid = 'GridAxesRepresentation'
ductctrlptpatchgridvtk_1Display.PolarAxes = 'PolarAxesRepresentation'
ductctrlptpatchgridvtk_1Display.ScalarOpacityUnitDistance = 1.1885640454213802

# update the view to ensure updated data information
renderView1.Update()

# change representation type
ductctrlptpatchgridvtk_1Display.SetRepresentationType('Points')

# change representation type
ductctrlptpatchgridvtk_1Display.SetRepresentationType('Wireframe')

# change solid color
ductctrlptpatchgridvtk_1Display.AmbientColor = [0.0, 0.0, 0.0]
ductctrlptpatchgridvtk_1Display.DiffuseColor = [0.0, 0.0, 0.0]

# create a new 'Legacy VTK Reader'
ductctrlptsvtk = LegacyVTKReader(FileNames=[vtkPtsFile])

# show data in view
ductctrlptsvtkDisplay = Show(ductctrlptsvtk, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
ductctrlptsvtkDisplay.Representation = 'Surface'
ductctrlptsvtkDisplay.ColorArrayName = [None, '']
ductctrlptsvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ductctrlptsvtkDisplay.SelectOrientationVectors = 'None'
ductctrlptsvtkDisplay.ScaleFactor = 1.2000000000000002
ductctrlptsvtkDisplay.SelectScaleArray = 'None'
ductctrlptsvtkDisplay.GlyphType = 'Arrow'
ductctrlptsvtkDisplay.GlyphTableIndexArray = 'None'
ductctrlptsvtkDisplay.GaussianRadius = 0.06
ductctrlptsvtkDisplay.SetScaleArray = [None, '']
ductctrlptsvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
ductctrlptsvtkDisplay.OpacityArray = [None, '']
ductctrlptsvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ductctrlptsvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
ductctrlptsvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
ductctrlptsvtkDisplay.ScalarOpacityUnitDistance = 5.513220050815476

# update the view to ensure updated data information
renderView1.Update()

# change representation type
ductctrlptsvtkDisplay.SetRepresentationType('Points')

# change solid color
ductctrlptsvtkDisplay.AmbientColor = [1.0, 0.0, 0.0]
ductctrlptsvtkDisplay.DiffuseColor = [1.0, 0.0, 0.0]

# Properties modified on ductctrlptsvtkDisplay
ductctrlptsvtkDisplay.LineWidth = 8.0

# Properties modified on ductctrlptsvtkDisplay
ductctrlptsvtkDisplay.PointSize = 8.0

# create a new 'Legacy VTK Reader'
ductctrlptsvtk_1 = LegacyVTKReader(FileNames=[vtkPtsFile])

# show data in view
ductctrlptsvtk_1Display = Show(ductctrlptsvtk_1, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
ductctrlptsvtk_1Display.Representation = 'Surface'
ductctrlptsvtk_1Display.ColorArrayName = [None, '']
ductctrlptsvtk_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
ductctrlptsvtk_1Display.SelectOrientationVectors = 'None'
ductctrlptsvtk_1Display.ScaleFactor = 1.2000000000000002
ductctrlptsvtk_1Display.SelectScaleArray = 'None'
ductctrlptsvtk_1Display.GlyphType = 'Arrow'
ductctrlptsvtk_1Display.GlyphTableIndexArray = 'None'
ductctrlptsvtk_1Display.GaussianRadius = 0.06
ductctrlptsvtk_1Display.SetScaleArray = [None, '']
ductctrlptsvtk_1Display.ScaleTransferFunction = 'PiecewiseFunction'
ductctrlptsvtk_1Display.OpacityArray = [None, '']
ductctrlptsvtk_1Display.OpacityTransferFunction = 'PiecewiseFunction'
ductctrlptsvtk_1Display.DataAxesGrid = 'GridAxesRepresentation'
ductctrlptsvtk_1Display.PolarAxes = 'PolarAxesRepresentation'
ductctrlptsvtk_1Display.ScalarOpacityUnitDistance = 5.513220050815476

# update the view to ensure updated data information
renderView1.Update()

# change representation type
ductctrlptsvtk_1Display.SetRepresentationType('Wireframe')

# change solid color
ductctrlptsvtk_1Display.AmbientColor = [1.0, 0.0, 0.0]
ductctrlptsvtk_1Display.DiffuseColor = [1.0, 0.0, 0.0]

# Properties modified on ductctrlptsvtk_1Display
ductctrlptsvtk_1Display.LineWidth = 2.0

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
