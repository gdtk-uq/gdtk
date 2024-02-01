# script-version: 2.0
from paraview.simple import *
from paraview import catalyst
import time

# registrationName must match the channel name used in the
# 'CatalystAdaptor'.
producer = TrivialProducer(registrationName="grid")
#groupDatasets1 = GroupDatasets(registrationName='GroupDatasets1', Input=producer)
#mergeBlocks1 = MergeBlocks(registrationName='MergeBlocks1', Input=groupDatasets1)
clip1 = Clip(registrationName='Clip1', Input=producer)
clip1.Invert = 0
clip1.ClipType.Origin = [0.0, 0.5, 0.0]
clip2 = Clip(registrationName='Clip2', Input=clip1)
clip2.ClipType.Origin = [4.0, 0.5, 0.0]
#cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=clip2)
#cleantoGrid1 = CleantoGrid(registrationName='CleantoGrid1', Input=cellDatatoPointData1)

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1600,800]
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [2.000000476837158, 0.12827871219978723, 113.93579144477845]
renderView1.CameraFocalPoint = [2.000000476837158, 0.12827871219978723, 0.0]
renderView1.CameraParallelScale = 1.1814243632270562
renderView1.CameraFocalDisk = 1.0
#renderView1.CameraParallelScale = 54.99504523136608

# get color transfer function/color map for 'velocity'
velocityLUT = GetColorTransferFunction('temperature')
velocityLUT.RGBPoints = [0.0, 0.0, 0.0, 0.34902, 21.875, 0.039216, 0.062745, 0.380392, 43.75, 0.062745, 0.117647, 0.411765, 65.625, 0.090196, 0.184314, 0.45098, 87.5, 0.12549, 0.262745, 0.501961, 109.375, 0.160784, 0.337255, 0.541176, 131.25, 0.2, 0.396078, 0.568627, 153.125, 0.239216, 0.454902, 0.6, 175.0, 0.286275, 0.521569, 0.65098, 196.875, 0.337255, 0.592157, 0.701961, 218.75, 0.388235, 0.654902, 0.74902, 240.625, 0.466667, 0.737255, 0.819608, 262.5, 0.572549, 0.819608, 0.878431, 284.375, 0.654902, 0.866667, 0.909804, 306.25, 0.752941, 0.917647, 0.941176, 328.125, 0.823529, 0.956863, 0.968627, 350.0, 0.988235, 0.960784, 0.901961, 350.0, 0.941176, 0.984314, 0.988235, 364.0, 0.988235, 0.945098, 0.85098, 378.0, 0.980392, 0.898039, 0.784314, 393.75, 0.968627, 0.835294, 0.698039, 415.625, 0.94902, 0.733333, 0.588235, 437.5, 0.929412, 0.65098, 0.509804, 459.375, 0.909804, 0.564706, 0.435294, 481.25, 0.878431, 0.458824, 0.352941, 503.125, 0.839216, 0.388235, 0.286275, 525.0, 0.760784, 0.294118, 0.211765, 546.875, 0.701961, 0.211765, 0.168627, 568.75, 0.65098, 0.156863, 0.129412, 590.625, 0.6, 0.094118, 0.094118, 612.5, 0.54902, 0.066667, 0.098039, 634.375, 0.501961, 0.05098, 0.12549, 656.25, 0.45098, 0.054902, 0.172549, 678.125, 0.4, 0.054902, 0.192157, 700.0, 0.34902, 0.070588, 0.211765]
velocityLUT.ScalarRangeInitialized = 1.0
velocityLUT.RescaleTransferFunction(0.0, 700.0)
velocityLUT.ColorSpace = 'Lab'
velocityLUT.NanColor = [0.25, 0.0, 0.0]
velocityLUT.ScalarRangeInitialized = 1.0

# show data from grid
gridDisplay = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

gridDisplay.Representation = 'Surface'
gridDisplay.ColorArrayName = ['CELLS', 'temperature']
gridDisplay.LookupTable = velocityLUT

# get color legend/bar for velocityLUT in view renderView1
velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)
velocityLUTColorBar.Title = 'Temperature (K)'
velocityLUTColorBar.ComponentTitle = ''

# change scalar bar placement
velocityLUTColorBar.Orientation = 'Horizontal'
velocityLUTColorBar.WindowLocation = 'Any Location'
velocityLUTColorBar.Position = [0.312338177014531, 0.90]
velocityLUTColorBar.ScalarBarLength = 0.3300000000000001
velocityLUTColorBar.AddRangeLabels = 0
velocityLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
velocityLUTColorBar.LabelColor = [0.0, 0.0, 0.0]


# Properties modified on renderView1
renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.OrientationAxesVisibility = 0

# set color bar visibility
velocityLUTColorBar.Visibility = 1

# show color legend
gridDisplay.SetScalarBarVisibility(renderView1, True)
gridDisplay.Ambient = 0.2

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

SetActiveView(renderView1)
# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'TimeStep'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'screenshot_{timestep:06d}.png'
pNG1.Writer.ImageResolution = [1600,800]
pNG1.Writer.Format = 'PNG'

# ------------------------------------------------------------------------------
# Catalyst options
options = catalyst.Options()
if "--enable-live" in catalyst.get_args():
  options.EnableCatalystLive = 1


# Greeting to ensure that ctest knows this script is being imported
#print("# executing catalyst_pipeline...")
def catalyst_execute(info):
    global producer
    producer.UpdatePipeline()
    #print("-----------------------------------")
    #print("executing (cycle={}, time={})".format(info.cycle, info.time))
    #print("bounds:", producer.GetDataInformation().GetBounds())
    #print("velocity-magnitude-range:", producer.CellData["velx"].GetRange(-1))
    #print("pressure-range:", producer.CellData["pressure"].GetRange(0))
    # In a real simulation sleep is not needed. We use it here to slow down the
    # "simulation" and make sure ParaView client can catch up with the produced
    # results instead of having all of them flashing at once.
    #if options.EnableCatalystLive:
    #    time.sleep(1)
