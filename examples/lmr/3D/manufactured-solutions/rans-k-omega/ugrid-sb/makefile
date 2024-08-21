# Look for pytest under different names
PYTEST := $(shell which pytest)
PYTEST := $(or $(PYTEST), pytest-3)

COMMON_DIR := ../common/
FILES_IN_COMMON := analytic_solution.py \
	constants.txt \
	fill-fn.lua \
	ref-soln.lua \
	udf-bc-ugrid.lua \
	udf-source-terms.lua \
	very-viscous-air.lua

prep-area:
	cd $(COMMON_DIR); cp $(FILES_IN_COMMON) $(CURDIR) 

prep-grid:
	lmr prep-grid

prep-sim:
	lmr prep-sim
run:
	lmr run

test:
	PYTHONDONTWRITEBYTEFILE=1 $(PYTEST) -v

clean:
	- rm -rf lmrsim
	- rm $(FILES_IN_COMMON)


