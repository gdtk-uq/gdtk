species-files := $(shell awk '$$1 !~ "\#" {print $$1}' < species-list.txt)

species-database.lua: species-list.txt $(species-files)
	$(info Generating species database...)
	@cat defaults.lua $(species-files) > $@
