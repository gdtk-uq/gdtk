# makefile for building the HTML Reference Manuals

LMR_DIR := eilmer
EXAMPLES_DIR := ../examples/eilmer
EILMER_DOC_FILES := $(LMR_DIR)/config-options-user.adoc \
	$(LMR_DIR)/flowsolution_user.adoc \
	$(LMR_DIR)/flow-state-user.adoc \
	$(LMR_DIR)/fluid-block-user.adoc \
	$(LMR_DIR)/boundary-conditions.adoc \
	$(LMR_DIR)/running-a-simulation.adoc \
	$(LMR_DIR)/run-time-utilities.adoc \
	$(LMR_DIR)/restarting-from-a-snapshot.adoc

GEOM_DIR := geom
GEOM_DOC_FILES := $(GEOM_DIR)/elements.adoc \
	$(GEOM_DIR)/sgrid.adoc

LIBRARY_DOC_FILES := flow/gas-flow.adoc \
	flow/ideal-gas.adoc \
	gas/gas-model.adoc \
	gas/kinetics.adoc \
	geom/elements.adoc \
	geom/sgrid.adoc

EILMER_CATALOG_ENTRIES := $(EXAMPLES_DIR)/2D/forward-facing-step/readme.adoc \
	$(EXAMPLES_DIR)/2D/sphere-nonaka/readme.adoc \
	$(EXAMPLES_DIR)/2D/foil-circular-arc-transonic/readme.adoc \
	$(EXAMPLES_DIR)/2D/flat-plate-hakkinen-swbli/readme.adoc \
	$(EXAMPLES_DIR)/2D/nozzle-conical-back/readme.adoc


all: eilmer-reference-manual.html \
	eilmer-catalog-of-examples.html \
	eilmer-faq.html \
	l1d-reference-manual.html \
	library-reference-manual.html \
	geometry-reference-manual.html \
	nm-reference-manual.html \
	estcn-manual.html \
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

build-lut-manual.html: build-lut-manual.adoc
	asciidoctor build-lut-manual.adoc

library-reference-manual.html: library-reference-manual.adoc $(LIBRARY_DOC_FILES)
	asciidoctor library-reference-manual.adoc

geometry-reference-manual.html: geometry-reference-manual.adoc $(GEOM_DOC_FILES)
	asciidoctor geometry-reference-manual.adoc

nm-reference-manual.html: nm-reference-manual.adoc
	asciidoctor nm-reference-manual.adoc

clean:
	- rm eilmer-reference-manual.html \
		eilmer-catalog-of-examples.html \
		eilmer-faq.html \
		l1d-reference-manual.html \
		library-reference-manual.html \
		geometry-reference-manual.html \
		nm-reference-manual.html \
		estcn-manual.html \
		build-lut-manual.html
