-- 2D Zero Pressure Gradient High Mach Number Flat Plate Validation Case
-- https://turbmodels.larc.nasa.gov/ZPGflatplateSS_val_sa.html
-- @author: Nick N. Gibbons (n.gibbons@uq.edu.au)
--
-- 2025-01-18, RJG
--    Ported to Eilmer 5
--
config.dimensions = 2
gridFile = "flatplate_clust2.p2dfmt.gz"
nptsInX = 545
nptsInY = 385
nptsOnSolidPlate = 449

-- 1. Read in grid from Plot3D format, as gzipped file (arg4 = true)
singleGrid = importPlot3DGrid(gridFile, config.dimensions, 1.0, true)

-- 2. Grid is a single piece. Now split at x = 0
--    We can use the information about number of points on the
--    solid plate to determine the break point in the grid.
--    Then we form two subgrids from the single grid.
nptsOnSymm = nptsInX - nptsOnSolidPlate + 1
subgrid0 = singleGrid[1]:subgrid(0, nptsOnSymm, 0, nptsInY)
subgrid1 = singleGrid[1]:subgrid(nptsOnSymm-1, nptsOnSolidPlate, 0, nptsInY)

-- 3. Use a GridArray to chop these up for parallel processing on 10 processors
registerFluidGridArray{grid=subgrid0, nib=1, njb=2, fsTag='initial',
   bcTags={north='supersonic', west='supersonic', south='symm'}}
registerFluidGridArray{grid=subgrid1, nib=4, njb=2, fsTag='initial',
   bcTags={north='supersonic', east='outflow', south='fixedT'}}
identifyGridConnections()
