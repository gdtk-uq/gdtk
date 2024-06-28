# makefile for building the HTML Reference Manuals

LMR_DIR := eilmer
GEOM_DIR := geom
EXAMPLES_DIR := ../examples/eilmer
EILMER_DOC_FILES := $(LMR_DIR)/config-options-user.adoc \
	$(LMR_DIR)/flowsolution_user.adoc \
	$(LMR_DIR)/flow-state-user.adoc \
	$(LMR_DIR)/fluid-block-user.adoc \
	$(LMR_DIR)/boundary-conditions.adoc \
	$(LMR_DIR)/running-a-simulation.adoc \
	$(LMR_DIR)/run-time-utilities.adoc \
	$(LMR_DIR)/restarting-from-a-snapshot.adoc \
	$(LMR_DIR)/extract-boundary-loads.adoc \
	$(GEOM_DIR)/surfaces/surfaces.adoc     \
	$(GEOM_DIR)/cluster_functions/cluster_functions.adoc


GEOM_DOC_FILES := $(GEOM_DIR)/elements.adoc \
	$(GEOM_DIR)/sgrid.adoc

LIBRARY_DOC_FILES := flow/gas-flow.adoc \
	flow/ideal-gas.adoc \
	gas/gas-model.adoc \
	gas/kinetics.adoc \
	geom/elements.adoc \
	geom/sgrid.adoc

EILMER_CATALOG_ENTRIES := \
	$(EXAMPLES_DIR)/2D/foil-circular-arc-transonic/readme.adoc \
	$(EXAMPLES_DIR)/2D/sphere-lehr/m355/readme.adoc \
	$(EXAMPLES_DIR)/2D/nozzle-conical-back/readme.adoc \
    $(EXAMPLES_DIR)/2D/forward-facing-step/readme.adoc \
    $(EXAMPLES_DIR)/2D/richtmyer-meshkov/readme.adoc \
    $(EXAMPLES_DIR)/2D/underexpanded-jet/readme.adoc \
    $(EXAMPLES_DIR)/2D/shear-layer-periodic/readme.adoc \
	$(EXAMPLES_DIR)/2D/flat-plate-hakkinen-swbli/readme.adoc \
	$(EXAMPLES_DIR)/2D/cylinder-type-iv-shock/readme.adoc \
	$(EXAMPLES_DIR)/2D/flat-plate-turbulent-larc/nk-5.45Tw-sa/readme.adoc \
	$(EXAMPLES_DIR)/2D/oblique-detonation-wave/readme.adoc \
	$(EXAMPLES_DIR)/2D/moving-grid/piston-w-const-vel/simple/readme.adoc \
	$(EXAMPLES_DIR)/2D/sphere-nonaka/readme.adoc \
	$(EXAMPLES_DIR)/2D/wall-catalysis/readme.adoc \
	$(EXAMPLES_DIR)/2D/solid/hollow-cylinder/readme.adoc


all: eilmer-reference-manual.html \
	eilmer-catalog-of-examples.html \
	eilmer-faq.html \
	l1d-reference-manual.html \
	library-reference-manual.html \
	geometry-reference-manual.html \
	nm-reference-manual.html \
	estcn-manual.html \
	nenzf1d-manual.html \
	build-lut-manual.html
	@echo "Done."

eilmer-reference-manual.html: eilmer-reference-manual.adoc $(EILMER_DOC_FILES)
	asciidoctor eilmer-reference-manual.adoc

eilmer-catalog-of-examples.html: eilmer-catalog-of-examples.adoc $(EILMER_CATALOG_ENTRIES)
	asciidoctor eilmer-catalog-of-examples.adoc

eilmer-faq.html: eilmer-faq.adoc
	asciidoctor eilmer-faq.adoc

l1d-reference-manual.html: l1d-reference-manual.adoc
	asciidoctor l1d-reference-manual.adoc

estcn-manual.html: estcn-manual.adoc
	asciidoctor estcn-manual.adoc

nenzf1d-manual.html: nenzf1d-manual.adoc
	asciidoctor nenzf1d-manual.adoc

build-lut-manual.html: build-lut-manual.adoc
	asciidoctor build-lut-manual.adoc

library-reference-manual.html: library-reference-manual.adoc $(LIBRARY_DOC_FILES)
	asciidoctor library-reference-manual.adoc

geometry-reference-manual.html: geometry-reference-manual.adoc $(GEOM_DOC_FILES)
	asciidoctor geometry-reference-manual.adoc

nm-reference-manual.html: nm-reference-manual.adoc
	asciidoctor nm-reference-manual.adoc

lmr-reference-manual.html : lmr-reference-manual.adoc
	cd lmr/commands; make
	asciidoctor lmr-reference-manual.adoc

clean:
	- rm eilmer-reference-manual.html \
		eilmer-catalog-of-examples.html \
		eilmer-faq.html \
		l1d-reference-manual.html \
		library-reference-manual.html \
		geometry-reference-manual.html \
		nm-reference-manual.html \
		estcn-manual.html \
		nenzf1d-manual.html \
		build-lut-manual.html
