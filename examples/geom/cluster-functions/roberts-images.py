import pyvista as pv

rbl = pv.read("roberts-cluster-left.vtk")
rbr = pv.read("roberts-cluster-right.vtk")
rbb = pv.read("roberts-cluster-both.vtk")

plt = pv.Plotter(off_screen=True, window_size=[2000,768], shape=(3,1), border=False)
# plot left-end cluster
plt.subplot(0, 0)
plt.add_mesh(rbl, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("RobertsFunction:new{end0=true, end1=false, beta=1.05}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

# Plot right end clustering
plt.subplot(1, 0)
plt.add_mesh(rbr, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("RobertsFunction:new{end0=false, end1=true, beta=1.05}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

# Plot iboth ends clustering
plt.subplot(2, 0)
plt.add_mesh(rbb, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("RobertsFunction:new{end0=true, end1=true, beta=1.05}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

#plt.show()
plt.screenshot("roberts-cf-images.png")
