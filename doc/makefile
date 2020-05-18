# makefile for building the HTML Reference Manuals

LMR_DIR := eilmer
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
	$(GEOM_DIR)/grid/grid-user.adoc \
	$(GEOM_DIR)/grid/sgrid-user.adoc

LIBRARY_DOC_FILES := flow/gas-flow.adoc \
	flow/ideal-gas.adoc \
	gas/gas-model.adoc \
	gas/kinetics.adoc \
	geom/elements.adoc

all: eilmer-reference-manual.html \
	l1d-reference-manual.html \
	library-reference-manual.html \
	geometry-reference-manual.html
	@echo "Done."

eilmer-reference-manual.html: eilmer-reference-manual.adoc $(EILMER_DOC_FILES)
	asciidoctor eilmer-reference-manual.adoc

l1d-reference-manual.html: l1d-reference-manual.adoc
	asciidoctor l1d-reference-manual.adoc

library-reference-manual.html: library-reference-manual.adoc $(LIBRARY_DOC_FILES)
	asciidoctor library-reference-manual.adoc

geometry-reference-manual.html: geometry-reference-manual.adoc $(GEOM_DOC_FILES)
	asciidoctor geometry-reference-manual.adoc

clean:
	- rm eilmer-reference-manual.html \
		l1d-reference-manual.html \
		library-reference-manual.html \
		geometry-reference-manual.html
