# gui_wdgets.py
"""
GUI widgets setup and interaction for the method of characteristics.
This is a regrowth of the old IMOC code that was implemented in C+Tcl/Tk.

Author(s):
  Peter J.
  Centre for Hypersonics,
  School of Mechanical Engineering, U of Q.

Version:
  2020-01-18: Rebuild of the old Tcl/Tk interface.
"""
# The following string will be displayed in the Help/About dialog.
versionString = \
"""Isentropic Method of Characteristics
Version: 0.1, 2020-Jan-18."""

import sys
import math
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox

import eilmer.imoc.kernel as kernel
import eilmer.imoc.unit_process as unit

cfg = {}
cfg["node_file"] = ""
cfg["plot_file"] = ""
# When displaying the characteristic mesh, we will give plotting
# commands to Tkinter in pixel units.
# A range of world units will be scaled to fit onto a patch of canvas.
# The scrolled-canvas widget will display a selected view-port
# of the full canvas.
cfg["canvas"] = None
cfg["dots_per_cm"] = 35; # guess for standard screens
cfg["vport_x_size"] = int(15.0 * cfg["dots_per_cm"])
cfg["vport_y_size"] = int(10.0 * cfg["dots_per_cm"])
cfg["canvas_x_size"] = None # will be set when the root window exists
cfg["canvas_y_size"] = None
cfg["canvas_x_offset"] = 40
cfg["canvas_y_offset"] = 40
cfg["same_scales"] = True
cfg["x_min"] = 0.0; cfg["x_max"] = 1.2
cfg["y_min"] = 0.0; cfg["y_max"] = 1.0
cfg['zoom_x1'] = None; cfg['zoom_y1'] = None

def set_xy_ranges(xmin=0.0, ymin=0.0, xmax=1.2, ymax=1.0):
    """
    Set the range of world units to display on the canvas.
    """
    if xmin > xmax: xmin, xmax = xmax, xmin
    if ymin > ymax: ymin, ymax = ymax, ymin
    cfg["x_min"] = xmin
    cfg["x_max"] = xmin+1.0 if xmax == xmin else xmax
    cfg["y_min"] = ymin
    cfg["y_max"] = ymin+1.0 if ymax == ymin else ymax
    set_xy_scales()
    return

def set_xy_scales():
    """
    Sets the world-to-canvas scales so that the full range will fit
    on the canvas with a few pixels to spare.
    """
    xscale = (cfg["canvas_x_size"]-cfg["canvas_x_offset"]-10)/(cfg["x_max"]-cfg["x_min"])
    yscale = (cfg["canvas_y_size"]-cfg["canvas_y_offset"]-10)/(cfg["y_max"]-cfg["y_min"])
    if cfg["same_scales"]:
        # Select the smaller scale, such that the full world-range fits.
        minscale = xscale if xscale < yscale else yscale
        cfg["x_scale"] = minscale; cfg["y_scale"] = minscale
    else:
        cfg["x_scale"] = xscale; cfg["y_scale"] = yscale
    return

def set_xy_tics(dx=0.2, dy=0.2):
    """
    Tic mark spacing in world units.
    """
    cfg["dx"] = dx; cfg["dy"] = dy
    return

def canvas_x(world_x):
    return (world_x - cfg["x_min"])*cfg["x_scale"] + cfg["canvas_x_offset"]

def canvas_y(world_y):
    return cfg["canvas_y_size"] - cfg["canvas_x_offset"] - \
        (world_y - cfg["y_min"])*cfg["y_scale"]

def world_x(canvas_x):
    return (canvas_x - cfg["canvas_x_offset"])/cfg["x_scale"] + cfg["x_min"]

def world_y(canvas_y):
    return (cfg["canvas_y_size"] - canvas_y - cfg["canvas_y_offset"])/cfg["y_scale"] \
        + cfg["y_min"]

cfg["show_node_numbers"] = False
cfg["show_char_mesh"] = True
cfg["show_stream_lines"] = True
cfg["display_dialog_for_coincident_nodes"] = True
cfg["pick_action"] = "select_node" # others select_corner_1, pick_corner_2
selected_nodes = []

def init_widgets():
    """
    By the end of this function, all of the GUI widgets will be in place.
    """
    # Some of our variables need to be visible to other functions.
    global root, cfg
    root = Tk()
    root.title("Isentropic Method of Characteristics")
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    cfg["canvas_x_size"] = int(0.7 * int(root.tk.eval("winfo screenwidth .")))
    cfg["canvas_y_size"] = int(0.7 * int(root.tk.eval("winfo screenheight .")))
    set_xy_ranges()
    set_xy_tics()
    root.geometry("+20+20") # place near top-left of screen
    #
    global var_g, var_p0, var_T0
    var_g = StringVar()
    var_p0 = StringVar()
    var_T0 = StringVar()
    global status_msg, coord_x, coord_y
    status_msg = StringVar(); status_msg.set("temporary message")
    coord_x = StringVar(); coord_x.set("0.0")
    coord_y = StringVar(); coord_y.set("0.0")
    #
    # mf (master frame) will resize when the user pulls the edges.
    mf = ttk.Frame(root, padding="3 3 3 3")
    mf.grid(column=0, row=0, sticky=(N, W, E, S))
    mf.columnconfigure(0, weight=1)
    mf.rowconfigure(0, weight=1)
    #
    # Scrolled Canvas for display of the MOC mesh.
    cf = ttk.Labelframe(mf, text="Mesh")
    h = ttk.Scrollbar(cf, orient=HORIZONTAL)
    v = ttk.Scrollbar(cf, orient=VERTICAL)
    c = Canvas(cf, width=cfg["vport_x_size"], height=cfg["vport_y_size"],
               scrollregion=(0, 0, cfg["canvas_x_size"], cfg["canvas_y_size"]),
               yscrollcommand=v.set, xscrollcommand=h.set)
    cfg["canvas"] = c # We will want access at event-processing times.
    c.yview_moveto(1.0-cfg["vport_y_size"]/cfg["canvas_x_size"]) # show lower region
    c.grid(column=0, row=0, sticky=(N,W,E,S))
    h.grid(column=0, row=1, sticky=(E,W))
    v.grid(column=1, row=0, sticky=(N,S))
    h['command'] = c.xview
    v['command'] = c.yview
    cf.columnconfigure(0, weight=1)
    cf.columnconfigure(1, weight=0)
    cf.rowconfigure(0, weight=1)
    cf.rowconfigure(1, weight=0)
    #
    # Status line.
    sf = ttk.Labelframe(mf, text="Status")
    sm = ttk.Label(sf, textvariable=status_msg, width=30)
    sm.grid(column=0, row=1, sticky=(W, E))
    xl = ttk.Label(sf, text="x:")
    xm = ttk.Label(sf, textvariable=coord_x, width=10)
    xl.grid(column=1, row=1, sticky=(W))
    xm.grid(column=2, row=1, sticky=(W))
    yl = ttk.Label(sf, text="y:")
    ym = ttk.Label(sf, textvariable=coord_y, width=10)
    yl.grid(column=3, row=1, sticky=(W))
    ym.grid(column=4, row=1, sticky=(W))
    sg = ttk.Sizegrip(sf)
    sg.grid(column=5, row=1, sticky=(S,E))
    sf.columnconfigure(0, weight=1)
    #
    cf.grid(column=0, row=0, sticky=(N,W,E,S))
    sf.grid(column=0, row=1, sticky=(W,E,S))
    #
    root.option_add('*tearOff', FALSE)
    # win = Toplevel(root)
    menubar = Menu(root)
    root['menu'] = menubar
    menu_file = Menu(menubar)
    menubar.add_cascade(menu=menu_file, label='File')
    menu_file.add_command(label='Source file...', command=sourceFile)
    menu_file.add_command(label='Quit', command=quitProgram)
    menu_edit = Menu(menubar)
    menubar.add_cascade(menu=menu_edit, label='Edit')
    menu_edit.add_command(label='Clear list of selected nodes', command=clearSelectedNodes)
    menu_edit.add_command(label='Delete selected nodes', command=deleteSelectedNodes)
    menu_compute = Menu(menubar)
    menubar.add_cascade(menu=menu_compute, label='Compute')
    menu_compute.add_command(label='Interior node', command=computeInteriorNode)
    menu_compute.add_command(label='C- to Wall 0', command=computeCminusWall0)
    menu_plot = Menu(menubar)
    menubar.add_cascade(menu=menu_plot, label='Plot')
    menu_plot.add_command(label='Refresh', command=refreshDisplay)
    menu_plot.add_command(label='Zoom', command=startZoom)
    menu_help = Menu(menubar)
    menubar.add_cascade(menu=menu_help, label='Help')
    menu_help.add_command(label='About', command=aboutProgram)
    #
    # Some bindings.
    c.bind("<Button-1>", pickSomething)
    c.bind("<Motion>", displayCoords)
    #
    # For the convenience of the user, leave the focus in the canvas widget
    c.focus()
    cfg["pick_action"] = "select_node"
    showStatusMsg("Ready to select node.")
    return

def sourceFile():
    filename = filedialog.askopenfilename()
    print("source file:", filename)
    f = open(filename, 'r')
    content = f.read()
    # We are going to trust the content of the file.
    exec(content, globals(), locals())
    return

def quitProgram():
    result = messagebox.askyesno(
        message="Quit program?",
        icon='question', title='IMOC',
        type='yesno')
    if result: sys.exit()
    return

def aboutProgram():
    result = messagebox.showinfo(
        message=versionString,
        icon='info', title='IMOC',
        type='ok')
    return

def showStatusMsg(text):
    status_msg.set(text)
    return

def pickSomething(event):
    """
    B1 event in the canvas picks one of:
        a node (this is the default pick_action)
        lower-left corner of zoom range
        upper-right corner of zoom-range
    """
    c = cfg['canvas']
    cx = c.canvasx(event.x); cy = c.canvasy(event.y)
    x = world_x(cx); y = world_y(cy)
    if cfg['pick_action'] == 'select_node':
        showStatusMsg("select node near x=%.3g y=%.3g" % (x, y))
        print("[TODO] finish pickSomething for selecting a node")
    elif cfg['pick_action'] == 'pick_corner_1':
        showStatusMsg("Lower-left corner for zoom x=%.3g y=%.3g" % (x, y))
        cfg['zoom_x1'] = x; cfg['zoom_y1'] = y
        cfg['pick_action'] = 'pick_corner_2'
    elif cfg['pick_action'] == 'pick_corner_2':
        showStatusMsg("Upper-right corner for zoom x=%.3g y=%.3g" % (x, y))
        set_xy_ranges(cfg['zoom_x1'], cfg['zoom_y1'], x, y)
        cfg['pick_action'] = 'select_node' # Return to default pick action.
        eraseCursors()
        refreshDisplay()
    else:
        print("Oops, should not have arrived here (in pickSomething).")
    return

def startZoom():
    cfg["pick_action"] = "pick_corner_1"
    return

def displayCoords(event):
    """
    On mouse movement in the canvas, update the displayed cursor position.
    Usually, this is just a text form in the status line, however,
    display cross-hair cursor if we are picking the corners of a zoom window.
    """
    c = cfg['canvas']
    cx = c.canvasx(event.x); cy = c.canvasy(event.y)
    x = world_x(cx); y = world_y(cy)
    coord_x.set("%.3g" % x); coord_y.set("%.3g" % y)

    c_x_min = c.canvasx(canvas_x(cfg['x_min']))
    c_y_min = c.canvasy(canvas_y(cfg['y_min']))
    c_x_max = c.canvasx(canvas_x(cfg['x_max']))
    c_y_max = c.canvasy(canvas_y(cfg['y_max']))
    if cfg['pick_action'] == 'pick_corner_1':
        c.delete("cursor1")
        c.create_line(cx, c_y_min, cx, c_y_max, fill='red', tags='cursor1')
        c.create_line(c_x_min, cy, c_x_max, cy, fill='red', tags='cursor1')
    if cfg['pick_action'] == 'pick_corner_2':
        c.delete("cursor2")
        c.create_line(cx, c_y_min, cx, c_y_max, fill='red', tags='cursor2')
        c.create_line(c_x_min, cy, c_x_max, cy, fill='red', tags='cursor2')
    return

def eraseCursors():
    c = cfg['canvas']
    c.delete('cursor1')
    c.delete('cursor2')
    return

def plotAxes():
    """
    Draw a set of axes and tic marks onto the canvas.
    """
    c = cfg['canvas']
    xmin = cfg['x_min']; xmax = cfg['x_max']
    ymin = cfg['y_min']; ymax = cfg['y_max']
    dx = cfg['dx']; dy = cfg['dy']
    #
    c.delete('axes')
    x1 = canvas_x(xmin); y1 = canvas_y(ymin)
    x2 = canvas_x(xmax); y2 = canvas_y(ymin)
    c.create_line(x1, y1, x2, y2, fill='black', tags='axes')
    x2 = canvas_x(xmin); y2 = canvas_y(ymax)
    c.create_line(x1, y1, x2, y2, fill='black', tags='axes')
    #
    eps = 1.0e-6 # something small
    x = xmin
    while x <= xmax+eps:
        x1 = canvas_x(x); y1 = canvas_y(ymin)
        c.create_line(x1, y1, x1, y1+5, fill='black', tags='axes')
        c.create_text(x1, y1+10, text="%.3g"%x, fill='black', tags='axes', anchor=N)
        x += cfg['dx']
    y = ymin
    while y <= ymax+eps:
        x1 = canvas_x(xmin); y1 = canvas_y(y)
        c.create_line(x1, y1, x1-5, y1, fill='black', tags='axes')
        c.create_text(x1-10, y1, text="%.3g"%y, fill='black', tags='axes', anchor=E)
        y += cfg['dy']
    return

def eraseAxes():
    cfg['canvas'].delete('axes')
    return

def plotWalls():
    c = cfg['canvas']
    c.delete('walls')
    xmin = cfg['x_min']; xmax = cfg['x_max']
    n = 100
    dx = (xmax - xmin)/n
    for wall in kernel.walls:
        x = xmin
        x1 = canvas_x(x); y1 = canvas_y(wall(x))
        for i in range(n):
            xnew = x+dx
            x2 = canvas_x(xnew); y2 = canvas_y(wall(xnew))
            c.create_line(x1, y1, x2, y2, fill='blue', tags='walls')
            x, x1, y1 = xnew, x2, y2
    return

def deleteSelectedNodes():
    print("[TODO] complete deleteSelectedNodes")
    return

def clearSelectedNodes():
    print("[TODO] clearSelectedNodes")
    return

def computeInteriorNode():
    print("[TODO] computeInteriorNode")
    return

def computeCminusWall0():
    print("[TODO] complete computeCminusWall0")
    return

def refreshDisplay():
    print("[TODO] complete refreshDisplay")
    plotAxes()
    plotWalls()
    return
