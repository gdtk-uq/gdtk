# Run script for CUDA D Hello World example

# We need to point our linking path to this directory in order
# to find the cuda kernel shared object
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. ./hello
