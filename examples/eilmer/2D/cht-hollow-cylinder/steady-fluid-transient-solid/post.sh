# post-process simulation
e4shared --job=cylinder --post --tindx-plot=last --vtk-xml --add-vars="mach"
e4shared --job=cylinder --post --tindx-plot=last --vtk-xml --plotTag="residual"
