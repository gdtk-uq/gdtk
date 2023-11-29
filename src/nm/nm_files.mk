NM_DIR ?= .
NM_FILES := $(NM_DIR)/package.d \
	$(NM_DIR)/nm_exception.d \
	$(NM_DIR)/number.d \
	$(NM_DIR)/bbla.d \
	$(NM_DIR)/bdfLU.d \
	$(NM_DIR)/bracketing.d \
	$(NM_DIR)/brent.d \
	$(NM_DIR)/secant.d \
	$(NM_DIR)/gaussquad.d \
	$(NM_DIR)/linesearch.d \
	$(NM_DIR)/nelmin.d \
	$(NM_DIR)/newton.d \
	$(NM_DIR)/newtoncotes.d \
	$(NM_DIR)/ridder.d \
	$(NM_DIR)/rungekutta.d \
	$(NM_DIR)/rsla.d \
	$(NM_DIR)/schedule.d \
	$(NM_DIR)/smla.d \
	$(NM_DIR)/stmatrix.d \
	$(NM_DIR)/tree_patch.d \
	$(NM_DIR)/univariate_lut.d \
	$(NM_DIR)/limiters.d \
	$(NM_DIR)/spline.d \
	$(NM_DIR)/splinelsq.d

NM_LUA_FILES := $(NM_DIR)/luabbla.d
