import pyvista as pv

gaussl = pv.read("gaussian-cluster-left.vtk")
gaussm = pv.read("gaussian-cluster-mid.vtk")
gaussr = pv.read("gaussian-cluster-right.vtk")

plt = pv.Plotter(off_screen=True, window_size=[2000,768], shape=(3,1), border=False)
# plot left-end cluster
plt.subplot(0, 0)
plt.add_mesh(gaussl, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("GaussianFunction:new{m=0.1, s=0.1, ratio=0.1}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

# Plot right end clustering
plt.subplot(1, 0)
plt.add_mesh(gaussm, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("GaussianFunction:new{m=0.5, s=0.2, ratio=0.1}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

# Plot iboth ends clustering
plt.subplot(2, 0)
plt.add_mesh(gaussr, line_width=4, color="black", point_size=12,  show_vertices=True, render_points_as_spheres=True)
plt.add_text("GaussianFunction:new{m=0.75, s=0.1, ratio=0.2}", position=(155,170), font='courier', color="black", font_size=18)
plt.add_text("0", position=(295,80), font_size=14, font='courier')
plt.add_text("1", position=(1690,80), font_size=14, font='courier')
plt.view_xy()
plt.camera.zoom(5.5)

#plt.show()
plt.screenshot("gaussian-cf-images.png")
