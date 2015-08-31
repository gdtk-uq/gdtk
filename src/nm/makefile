# makefile for the Numerical Methods package
# Builds all of the demonstration programs.

DEMO_PROGRAMS := ridder_demo linesearch_demo nelmin_demo rungekutta_demo \
	newtoncotes_demo gaussquad_demo bbla_demo

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

ridder_demo: ridder_demo.d ridder.d
	dmd -unittest -g ridder_demo.d

# Don't know why we have to include linesearch.d on the command line.
linesearch_demo: linesearch_demo.d linesearch.d
	dmd -unittest -g linesearch_demo.d linesearch.d

nelmin_demo: nelmin_demo.d nelmin.d
	dmd -unittest -g nelmin_demo.d nelmin.d

rungekutta_demo: rungekutta_demo.d rungekutta.d
	dmd -unittest -g rungekutta_demo.d rungekutta.d

newtoncotes_demo: newtoncotes_demo.d newtoncotes.d
	dmd -unittest -g newtoncotes_demo.d newtoncotes.d

gaussquad_demo: gaussquad_demo.d gaussquad.d
	dmd -unittest -g gaussquad_demo.d gaussquad.d

bbla_demo: bbla_demo.d bbla.d
	dmd -unittest -g bbla_demo.d bbla.d
