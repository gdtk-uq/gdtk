all: fire2-examples \
	nozzle-examples \
	duct-examples

fire2-examples: fire2-coons-patch-grid.png \
	fire2-AO-patch-grid.png \
	fire2-channel-patch-grid.png \
	fire2-ruled-surface-grid.png \
	fire2-ctrl-pt-patch-grid.png \
	fire2-ctrl-pt-patch-grid-p-pts.png

nozzle-examples: nozzle-coons-patch-grid.png \
	nozzle-AO-patch-grid.png \
	nozzle-AO-patch-grid-2.png \
	nozzle-channel-patch-grid.png \
	nozzle-ruled-surface-grid.png \
	nozzle-ctrl-pt-patch-grid.png \
	nozzle-ctrl-pt-patch-grid-p-pts.png

duct-examples: duct-coons-patch-grid.png \
	duct-AO-patch-grid.png \
	duct-ctrl-pt-patch-grid.png \
	duct-ctrl-pt-patch-grid-p-pts.png

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


nozzle-coons-patch-grid.png : nozzle-coons-patch-grid.vtk
	pvbatch create-nozzle-grid-image.py nozzle-coons-patch-grid.vtk nozzle-coons-patch-grid.png
	mogrify -trim nozzle-coons-patch-grid.png

nozzle-coons-patch-grid.vtk : nozzle-examples.lua
	e4shared --custom-script --script-file=nozzle-examples.lua

nozzle-AO-patch-grid.png : nozzle-AO-patch-grid.vtk
	pvbatch create-nozzle-grid-image.py nozzle-AO-patch-grid.vtk nozzle-AO-patch-grid.png
	mogrify -trim nozzle-AO-patch-grid.png

nozzle-AO-patch-grid.vtk : nozzle-examples.lua
	e4shared --custom-script --script-file=nozzle-examples.lua

nozzle-AO-patch-grid-2.png : nozzle-AO-patch-grid-2.vtk
	pvbatch create-nozzle-grid-image.py nozzle-AO-patch-grid-2.vtk nozzle-AO-patch-grid-2.png
	mogrify -trim nozzle-AO-patch-grid-2.png

nozzle-AO-patch-grid-2.vtk : nozzle-examples.lua
	e4shared --custom-script --script-file=nozzle-examples.lua

nozzle-channel-patch-grid.png : nozzle-channel-patch-grid.vtk
	pvbatch create-nozzle-grid-image.py nozzle-channel-patch-grid.vtk nozzle-channel-patch-grid.png
	mogrify -trim nozzle-channel-patch-grid.png

nozzle-channel-patch-grid.vtk : nozzle-examples.lua
	e4shared --custom-script --script-file=nozzle-examples.lua

nozzle-ruled-surface-grid.png : nozzle-ruled-surface-grid.vtk
	pvbatch create-nozzle-grid-image.py nozzle-ruled-surface-grid.vtk nozzle-ruled-surface-grid.png
	mogrify -trim nozzle-ruled-surface-grid.png

nozzle-ruled-surface-grid.vtk : nozzle-examples.lua
	e4shared --custom-script --script-file=nozzle-examples.lua

nozzle-ctrl-pt-patch-grid.png : nozzle-ctrl-pt-patch-grid.vtk
	pvbatch create-nozzle-grid-image.py nozzle-ctrl-pt-patch-grid.vtk nozzle-ctrl-pt-patch-grid.png
	mogrify -trim nozzle-ctrl-pt-patch-grid.png

nozzle-ctrl-pt-patch-grid.vtk : nozzle-examples.lua
	e4shared --custom-script --script-file=nozzle-examples.lua

nozzle-ctrl-pt-patch-grid-p-pts.png : nozzle-ctrl-pt-patch-grid.vtk nozzle-ctrl-pts.vtk
	pvbatch create-nozzle-ctrl-pt-grid-image.py nozzle-ctrl-pt-patch-grid.vtk nozzle-ctrl-pts.vtk nozzle-ctrl-pt-patch-grid-p-pts.png
	mogrify -trim nozzle-ctrl-pt-patch-grid-p-pts.png

duct-coons-patch-grid.png : duct-coons-patch-grid.vtk
	pvbatch create-duct-grid-image.py duct-coons-patch-grid.vtk duct-coons-patch-grid.png
	mogrify -trim duct-coons-patch-grid.png

duct-coons-patch-grid.vtk : duct-examples.lua
	e4shared --custom-script --script-file=duct-examples.lua

duct-AO-patch-grid.png : duct-AO-patch-grid.vtk
	pvbatch create-duct-grid-image.py duct-AO-patch-grid.vtk duct-AO-patch-grid.png
	mogrify -trim duct-AO-patch-grid.png

duct-AO-patch-grid.vtk : duct-examples.lua
	e4shared --custom-script --script-file=duct-examples.lua

duct-ctrl-pt-patch-grid.png : duct-ctrl-pt-patch-grid.vtk
	pvbatch create-duct-grid-image.py duct-ctrl-pt-patch-grid.vtk duct-ctrl-pt-patch-grid.png
	mogrify -trim duct-ctrl-pt-patch-grid.png

duct-ctrl-pt-patch-grid.vtk : duct-examples.lua
	e4shared --custom-script --script-file=duct-examples.lua

duct-ctrl-pt-patch-grid-p-pts.png : duct-ctrl-pt-patch-grid.vtk duct-ctrl-pts.vtk
	pvbatch create-duct-ctrl-pt-grid-image.py duct-ctrl-pt-patch-grid.vtk duct-ctrl-pts.vtk duct-ctrl-pt-patch-grid-p-pts.png
	mogrify -trim duct-ctrl-pt-patch-grid-p-pts.png
