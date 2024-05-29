# makefile for the Chapel tests
# PJ 2024-01-28 Now that we will have more than one test, let's tidy up.

test: vector3_test rsla_test spline_test \
	gas_test flow_test face_test cell_test
	./vector3_test --numLocales=1
	./rsla_test --numLocales=1
	./spline_test --numLocales=1
	./gas_test --numLocales=1
	./flow_test --numLocales=1
	./face_test --numLocales=1
	./cell_test --numLocales=1

clean:
	- rm ./vector3_test ./rsla_test ./spline_test
	- rm ./gas_test ./flow_test ./face_test ./cell_test
	- rm ./*_test_real

vector3_test: vector3_test.chpl ../vector3.chpl
	chpl vector3_test.chpl ../vector3.chpl

rsla_test: rsla_test.chpl ../rsla.chpl
	chpl rsla_test.chpl ../rsla.chpl

spline_test: spline_test.chpl ../spline.chpl ../rsla.chpl
	chpl spline_test.chpl ../spline.chpl ../rsla.chpl

gas_test: gas_test.chpl ../gas.chpl
	chpl gas_test.chpl ../gas.chpl

flow_test: flow_test.chpl ../vector3.chpl ../gas.chpl ../flow.chpl
	chpl flow_test.chpl ../vector3.chpl ../gas.chpl ../flow.chpl

face_test: face_test.chpl ../face.chpl ../flow.chpl ../gas.chpl ../vector3.chpl
	chpl face_test.chpl ../face.chpl ../flow.chpl ../gas.chpl ../vector3.chpl

cell_test: cell_test.chpl ../cell.chpl ../face.chpl ../flow.chpl ../gas.chpl ../vector3.chpl
	chpl cell_test.chpl ../cell.chpl ../face.chpl ../flow.chpl ../gas.chpl ../vector3.chpl
