species-database.lua : species-list.txt *.lua
	awk '$$1 !~ "#"  {print $$1}' < species-list.txt | xargs cat defaults.lua > $@	
