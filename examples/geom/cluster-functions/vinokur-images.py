import pyvista as pv

vin0 = pv.read("vinokur-cluster-ex0.vtk")
vin1 = pv.read("vinokur-cluster-ex1.vtk")
vin2 = pv.read("vinokur-cluster-ex2.vtk")

plt = pv.Plotter(off_screen=True, window_size=[2000,768], shape=(3,1), border=False)
# plot left-end cluster
plt.subplot(0, 0)
plt.add_mesh(vin0, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("VinokurFunction:new{n=15, s0=0.01, s1=0.02}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

# Plot right end clustering
plt.subplot(1, 0)
plt.add_mesh(vin1, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("VinokurFunction:new{n=15, s0=0.005, s1=0.02}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

# Plot iboth ends clustering
plt.subplot(2, 0)
plt.add_mesh(vin2, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("VinokurFunction:new{n=15, s0=0.005, s1=0.03}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

#plt.show()
plt.screenshot("vinokur-cf-images.png")
