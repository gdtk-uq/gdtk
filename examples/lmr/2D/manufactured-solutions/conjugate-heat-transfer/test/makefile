help:
	echo "Main targets are:"
	echo "   prep"
	echo "   run"
	echo "   post"

grid: job.lua
	lmr prep-grid

init: job.lua
	lmr prep-sim

prep: grid init

run:
	lmr run

post:
	lmr snapshot2vtk --all --subtract-ref-solution=ref-soln.lua --subtract-solid-ref-solution=ref-soln.lua


