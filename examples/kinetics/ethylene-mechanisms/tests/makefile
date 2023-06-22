BASEDIR := ..

help:
	@echo "Main targets available are:"
	@echo -e "\t prep"
	@echo -e "\t test"
	@echo -e "\t plot"


vw-2002.gas: $(BASEDIR)/varatharajan-williams-2002-gas-model.lua
	prep-gas $(BASEDIR)/varatharajan-williams-2002-gas-model.lua vw-2002.gas

vw-2002.chem: vw-2002.gas $(BASEDIR)/varatharajan-williams-2002-short-mechanism.lua
	prep-chem vw-2002.gas $(BASEDIR)/varatharajan-williams-2002-short-mechanism.lua vw-2002.chem

prep: vw-2002.gas vw-2002.chem

test: isobaric-reactor.lua vw-2002.gas vw-2002.chem
	gas-calc isobaric-reactor.lua

isobaric-reactor.data : test

plot: isobaric-reactor.data
	gnuplot plot-history.gplot

clean:
	-rm vw-2002.gas vw-2002.chem isobaric-reactor.data ethylene-ignition-history.pdf



