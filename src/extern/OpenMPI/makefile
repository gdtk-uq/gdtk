# makefile for cleaning out OpenMPI package
# PJ 2017-06-05

dummy:
	echo "do nothing"

clean:
	- mv dub.json.back dub.json
	- rm -r .dub/
	- rm mpi_configure mpi_splice
	- rm source/mpi/package.d

