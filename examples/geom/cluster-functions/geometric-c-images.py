import pyvista as pv

geol = pv.read("geometric-cluster-left.vtk")
geol2 = pv.read("geometric-cluster-left2.vtk")
geor = pv.read("geometric-cluster-right.vtk")

plt = pv.Plotter(off_screen=True, window_size=[2000,768], shape=(3,1), border=False)
# plot left-end cluster
plt.subplot(0, 0)
plt.add_mesh(geol, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("GeometricFunction:new{a=0.002, r=1.1, N=60}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

# Plot right end clustering
plt.subplot(1, 0)
plt.add_mesh(geol2, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("GeometricFunction:new{a=0.004, r=1.2, N=60}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

# Plot iboth ends clustering
plt.subplot(2, 0)
plt.add_mesh(geor, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("GeometricFunction:new{a=0.004, r=1.2, N=60, reverse=true}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

#plt.show()
plt.screenshot("geometric-cf-images.png")
