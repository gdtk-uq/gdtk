all: fire2-examples

fire2-examples: fire2-coons-patch-grid.png \
	fire2-AO-patch-grid.png \
	fire2-channel-patch-grid.png \
	fire2-ruled-surface-grid.png \
	fire2-ctrl-pt-patch-grid.png \
	fire2-ctrl-pt-patch-grid-p-pts.png

fire2-coons-patch-grid.png : fire2-coons-patch-grid.vtk
	pvbatch create-fire2-grid.py fire2-coons-patch-grid.vtk fire2-coons-patch-grid.png
	mogrify -trim fire2-coons-patch-grid.png

fire2-coons-patch-grid.vtk : fire2-examples.lua shock-shape.dat
	e4shared --custom-script --script-file=fire2-examples.lua

fire2-AO-patch-grid.png : fire2-AO-patch-grid.vtk
	pvbatch create-fire2-grid.py fire2-AO-patch-grid.vtk fire2-AO-patch-grid.png
	mogrify -trim fire2-AO-patch-grid.png

fire2-AO-patch-grid.vtk : fire2-examples.lua shock-shape.dat
	e4shared --custom-script --script-file=fire2-examples.lua

fire2-channel-patch-grid.png : fire2-channel-patch-grid.vtk
	pvbatch create-fire2-grid.py fire2-channel-patch-grid.vtk fire2-channel-patch-grid.png
	mogrify -trim fire2-channel-patch-grid.png

fire2-channel-patch-grid.vtk : fire2-examples.lua shock-shape.dat
	e4shared --custom-script --script-file=fire2-examples.lua

fire2-ruled-surface-grid.png : fire2-ruled-surface-grid.vtk
	pvbatch create-fire2-grid.py fire2-ruled-surface-grid.vtk fire2-ruled-surface-grid.png
	mogrify -trim fire2-ruled-surface-grid.png

fire2-ruled-surface-grid.vtk : fire2-examples.lua shock-shape.dat
	e4shared --custom-script --script-file=fire2-examples.lua

fire2-ctrl-pt-patch-grid.png : fire2-ctrl-pt-patch-grid.vtk
	pvbatch create-fire2-grid.py fire2-ctrl-pt-patch-grid.vtk fire2-ctrl-pt-patch-grid.png
	mogrify -trim fire2-ctrl-pt-patch-grid.png

fire2-ctrl-pt-patch-grid.vtk : fire2-examples.lua shock-shape.dat
	e4shared --custom-script --script-file=fire2-examples.lua

fire2-ctrl-pt-patch-grid-p-pts.png : fire2-ctrl-pts.vtk fire2-ctrl-pt-patch-grid.vtk
	pvbatch create-fire2-ctrl-pt-grid.py fire2-ctrl-pt-patch-grid.vtk fire2-ctrl-pts.vtk fire2-ctrl-pt-patch-grid-p-pts.png
	mogrify -trim fire2-ctrl-pt-patch-grid-p-pts.png

fire2-ctrl-pts.vtk : fire2-examples.lua shock-shape.dat
	e4shared --custom-script --script-file=fire2-examples.lua

