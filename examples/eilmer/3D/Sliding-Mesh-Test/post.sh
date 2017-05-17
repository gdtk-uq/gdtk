#! /bin/sh
# run_Sliding_Mesh_test.sh
# Test Case for Developing a sliding mesh interface

# Postprocess and create .vtk
e4shared --post --job=Sliding_Mesh_Test --vtk-xml --tindx-plot=all

# Open paraview
paraview plot/Sliding_Mesh_Test.pvd




