# makefile for the Chapel tests
# PJ 2024-01-28 Now that we will have more than one test, let's tidy up.

test: vector3_test
	./vector3_test

clean:
	- rm ./vector3_test

vector3_test: vector3_test.chpl ../vector3.chpl
	chpl vector3_test.chpl ../vector3.chpl

