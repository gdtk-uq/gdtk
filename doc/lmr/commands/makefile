COMMANDS := prep-gas \
	prep-reactions \
	prep-energy-exchange \
	prep-grid \
	prep-sim \
	run \
	snapshot2vtk \
	probe-flow \
	slice-flow \
	extract-line \
	compute-norms
	

all: $(COMMANDS)

clean:
	- rm *.txt

$(COMMANDS):
	lmr help $@ > $@.txt

