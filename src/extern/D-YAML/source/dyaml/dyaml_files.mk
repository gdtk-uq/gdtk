# dyaml_files.mk
# PJ 2020-09-26

DYAML_DIR ?= .
DYAML_FILES := $(DYAML_DIR)/composer.d \
	$(DYAML_DIR)/constructor.d \
	$(DYAML_DIR)/dumper.d \
	$(DYAML_DIR)/emitter.d \
	$(DYAML_DIR)/encoding.d \
	$(DYAML_DIR)/escapes.d \
	$(DYAML_DIR)/event.d \
	$(DYAML_DIR)/exception.d \
	$(DYAML_DIR)/linebreak.d \
	$(DYAML_DIR)/loader.d \
	$(DYAML_DIR)/node.d \
	$(DYAML_DIR)/package.d \
	$(DYAML_DIR)/parser.d \
	$(DYAML_DIR)/queue.d \
	$(DYAML_DIR)/reader.d \
	$(DYAML_DIR)/representer.d \
	$(DYAML_DIR)/resolver.d \
	$(DYAML_DIR)/scanner.d \
	$(DYAML_DIR)/serializer.d \
	$(DYAML_DIR)/style.d \
	$(DYAML_DIR)/tagdirective.d \
	$(DYAML_DIR)/token.d

