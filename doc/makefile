# makefile for building the HTML users-reference-manual

LMR_DIR := eilmer
EILMER_DOC_FILES := $(LMR_DIR)/config-options-user.adoc \
	$(LMR_DIR)/flowsolution_user.adoc \
	$(LMR_DIR)/flow-state-user.adoc \
	$(LMR_DIR)/fluid-block-user.adoc \
	$(LMR_DIR)/monitoring-a-calculation.adoc \
	$(LMR_DIR)/restarting-from-a-snapshot.adoc

GEOM_DIR := geom
GEOM_DOC_FILES := (GEOM_DIR)/grid/grid-user.adoc \
	$(GEOM_DIR)/grid/sgrid-user.adoc

all: users-reference-manual.html l1d-reference-manual.html
	@echo "Done."

users-reference-manual.html: users-reference-manual.adoc $(EILMER_DOC_FILES)
	asciidoctor users-reference-manual.adoc

l1d-reference-manual.html: l1d-reference-manual.adoc
	asciidoctor l1d-reference-manual.adoc

clean:
	- rm users-reference-manual.html l1d-reference-manual.html
