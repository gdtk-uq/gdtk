# makefile for building the HTML users-reference-manual

LMR_DIR := eilmer
EILMER_DOC_FILES := $(LMR_DIR)/config-options-user.adoc \
	$(LMR_DIR)/flowsolution_user.adoc \
	$(LMR_DIR)/flow-state-user.adoc \
	$(LMR_DIR)/fluid-block-user.adoc

users-reference-manual.html: users-reference-manual.adoc $(EILMER_DOC_FILES)
	asciidoctor users-reference-manual.adoc

clean:
	- rm users-reference-manual.html
