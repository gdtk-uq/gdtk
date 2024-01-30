# makefile for the Chapel tests
# PJ 2024-01-28 Now that we will have more than one test, let's tidy up.

test: vector3_test gas_test flow_test face_test
	./vector3_test
	./gas_test
	./flow_test
	./face_test

clean:
	- rm ./vector3_test ./gas_test ./flow_test ./face_test

vector3_test: vector3_test.chpl ../vector3.chpl
	chpl vector3_test.chpl ../vector3.chpl

gas_test: gas_test.chpl ../gas.chpl
	chpl gas_test.chpl ../gas.chpl

flow_test: flow_test.chpl ../vector3.chpl ../gas.chpl ../flow.chpl
	chpl flow_test.chpl ../vector3.chpl ../gas.chpl ../flow.chpl

face_test: face_test.chpl ../face.chpl ../flow.chpl ../gas.chpl ../vector3.chpl
	chpl face_test.chpl ../face.chpl ../flow.chpl ../gas.chpl ../vector3.chpl
