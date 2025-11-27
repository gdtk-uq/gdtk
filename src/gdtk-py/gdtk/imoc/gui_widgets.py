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
Version: 0.2, 2020-Apr-19."""

import sys
import math
import collections
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox

import gdtk.imoc.kernel as kernel
import gdtk.imoc.unit_process as unit

cfg = {}
cfg['node_file'] = ""
cfg['plot_file'] = ""
# When displaying the characteristic mesh, we will give plotting
# commands to Tkinter in pixel units.
# A range of world units will be scaled to fit onto a patch of canvas.
# The scrolled-canvas widget will display a selected view-port
# of the full canvas.
cfg['canvas'] = None
cfg['dots_per_cm'] = 35; # guess for standard screens
cfg['vport_x_size'] = int(15.0 * cfg['dots_per_cm'])
cfg['vport_y_size'] = int(10.0 * cfg['dots_per_cm'])
cfg['canvas_x_size'] = None # will be set when the root window exists
cfg['canvas_y_size'] = None
cfg['canvas_x_offset'] = 40
cfg['canvas_y_offset'] = 40
cfg['equal_scales'] = True
cfg['x_min'] = 0.0; cfg['x_max'] = 1.2
cfg['y_min'] = 0.0; cfg['y_max'] = 1.0
cfg['zoom_x1'] = None; cfg['zoom_y1'] = None

def set_xy_ranges(xmin=0.0, ymin=0.0, xmax=1.2, ymax=1.0):
    """
    Set the range of world units to display on the canvas.
    """
    if xmin > xmax: xmin, xmax = xmax, xmin
    if ymin > ymax: ymin, ymax = ymax, ymin
    cfg['x_min'] = xmin
    cfg['x_max'] = xmin+1.0 if xmax == xmin else xmax
    cfg['y_min'] = ymin
    cfg['y_max'] = ymin+1.0 if ymax == ymin else ymax
    set_xy_scales()
    return

def set_xy_scales():
    """
    Sets the world-to-canvas scales so that the full range will fit
    on the canvas with a few pixels to spare.
    """
    xscale = (cfg['canvas_x_size']-cfg['canvas_x_offset']-10)/(cfg['x_max']-cfg['x_min'])
    yscale = (cfg['canvas_y_size']-cfg['canvas_y_offset']-10)/(cfg['y_max']-cfg['y_min'])
    if cfg['equal_scales']:
        # Select the smaller scale, such that the full world-range fits.
        minscale = xscale if xscale < yscale else yscale
        cfg['x_scale'] = minscale; cfg['y_scale'] = minscale
    else:
        cfg['x_scale'] = xscale; cfg['y_scale'] = yscale
    return

def set_xy_tics(dx=0.2, dy=0.2):
    """
    Tic mark spacing in world units.
    """
    cfg['dx'] = dx; cfg['dy'] = dy
    return

def canvas_x(world_x):
    return (world_x - cfg['x_min'])*cfg['x_scale'] + cfg['canvas_x_offset']

def canvas_y(world_y):
    return cfg['canvas_y_size'] - cfg['canvas_x_offset'] - \
        (world_y - cfg['y_min'])*cfg['y_scale']

def world_x(canvas_x):
    return (canvas_x - cfg['canvas_x_offset'])/cfg['x_scale'] + cfg['x_min']

def world_y(canvas_y):
    return (cfg['canvas_y_size'] - canvas_y - cfg['canvas_y_offset'])/cfg['y_scale'] \
        + cfg['y_min']

cfg['node_indx'] = {}
cfg['object_id'] = {}
cfg['display_dialog_for_coincident_nodes'] = True
# Make something like an enum for use with the pick_action state variable.
PickAction = collections.namedtuple('_', 'node corner_1 corner_2')(*range(3))
cfg['pick_action'] = PickAction.node
cfg['selected_nodes'] = []

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
    cfg['canvas_x_size'] = int(0.7 * int(root.tk.eval("winfo screenwidth .")))
    cfg['canvas_y_size'] = int(0.7 * int(root.tk.eval("winfo screenheight .")))
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
    global var_show_node_numbers, var_show_char_mesh, var_show_streamlines
    global var_equal_scales
    var_show_node_numbers = StringVar(); var_show_node_numbers.set(0)
    var_show_char_mesh = StringVar(); var_show_char_mesh.set(1)
    var_show_streamlines = StringVar(); var_show_streamlines.set(1)
    var_equal_scales = StringVar(); var_equal_scales.set(cfg['equal_scales'])
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
    c = Canvas(cf, width=cfg['vport_x_size'], height=cfg['vport_y_size'],
               scrollregion=(0, 0, cfg['canvas_x_size'], cfg['canvas_y_size']),
               yscrollcommand=v.set, xscrollcommand=h.set)
    cfg['canvas'] = c # We will want access at event-processing times.
    c.yview_moveto(1.0-cfg['vport_y_size']/cfg['canvas_x_size']) # show lower region
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
    menu_file.add_command(label='Save node data...', command=saveDataFile)
    menu_file.add_command(label='Save view as Postscript...', command=savePostscriptFile)
    menu_file.add_command(label='Quit', command=quitProgram)
    menu_edit = Menu(menubar)
    menubar.add_cascade(menu=menu_edit, label='Edit')
    menu_edit.add_command(label='Clear list of selected nodes',
                          command=clearSelectedNodesList)
    menu_edit.add_command(label='Delete selected nodes from mesh',
                          command=deleteSelectedNodesFromMesh)
    menu_compute = Menu(menubar)
    menubar.add_cascade(menu=menu_compute, label='Compute')
    menu_compute.add_command(label='Interior node', command=computeInteriorNode)
    menu_compute.add_command(label='C- to Wall 0', command=computeCminusWall0)
    menu_plot = Menu(menubar)
    menubar.add_cascade(menu=menu_plot, label='Plot')
    menu_plot.add_command(label='Refresh', command=refreshDisplay)
    menu_plot.add_command(label='Zoom to cursor selection', command=startZoom)
    menu_plot.add_command(label='Zoom to include all', command=zoomToIncludeAll)
    menu_plot.add_checkbutton(label='Equal x,y-scales', variable=var_equal_scales,
                              onvalue=1, offvalue=0, command=adjustScalesRefreshDisplay)
    menu_plot.add_checkbutton(label='Show node numbers', variable=var_show_node_numbers,
                              onvalue=1, offvalue=0, command=refreshDisplay)
    menu_plot.add_checkbutton(label='Show char mesh', variable=var_show_char_mesh,
                              onvalue=1, offvalue=0, command=refreshDisplay)
    menu_plot.add_checkbutton(label='Show streamlines', variable=var_show_streamlines,
                              onvalue=1, offvalue=0, command=refreshDisplay)
    menu_help = Menu(menubar)
    menubar.add_cascade(menu=menu_help, label='Help')
    menu_help.add_command(label='About', command=aboutProgram)
    #
    # Some bindings.
    c.bind('<Button-1>', pickSomething)
    c.bind('<Motion>', displayCoords)
    # Pressing ESC key clears list of selected nodes.
    root.bind('<Escape>', clearSelectedNodesList)
    #
    # For the convenience of the user, leave the focus in the canvas widget
    c.focus()
    cfg['pick_action'] = PickAction.node
    showStatusMsg("Ready.")
    return

def sourceFile():
    filename = filedialog.askopenfilename(title="Open user script")
    print("source file:", filename)
    f = open(filename, 'r')
    content = f.read()
    # We are going to trust the content of the file.
    exec(content, globals(), locals())
    refreshDisplay()
    return

def saveDataFile():
    filename = filedialog.asksaveasfilename(title="Save node data as")
    print("data file:", filename)
    f = open(filename, 'w')
    for node in kernel.nodes:
        f.write("%s\n" % str(node))
    f.write("char_mesh=[\n")
    n = len(kernel.char_mesh)
    for i in range(n):
        f.write("%d" % kernel.char_mesh[i])
        if i+1 < n:
            f.write(",")
            ws = "\n" if (i+1) % 10 == 0 else " "
            f.write(ws)
    f.write("\n]\n")
    return

def savePostscriptFile():
    filename = filedialog.asksaveasfilename(title="Save view as Postscript")
    print("postscript file:", filename)
    c = cfg['canvas']
    c.postscript(file=filename)
    return

def quitProgram():
    result = messagebox.askyesno(
        message="Quit program?",
        icon='question', title='IMOC',
        type='yesno')
    if result: root.destroy()
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
    if cfg['pick_action'] == PickAction.node:
        showStatusMsg("Try to select node at x=%.3g y=%.3g" % (x, y))
        selectNode(c, cx, cy)
    elif cfg['pick_action'] == PickAction.corner_1:
        showStatusMsg("Lower-left corner for zoom x=%.3g y=%.3g" % (x, y))
        cfg['zoom_x1'] = x; cfg['zoom_y1'] = y
        cfg['pick_action'] = PickAction.corner_2
    elif cfg['pick_action'] == PickAction.corner_2:
        showStatusMsg("Upper-right corner for zoom x=%.3g y=%.3g" % (x, y))
        set_xy_ranges(cfg['zoom_x1'], cfg['zoom_y1'], x, y)
        cfg['pick_action'] = PickAction.node # Return to default pick action.
        eraseCursors()
        refreshDisplay()
        showStatusMsg("Zoom to selected range x=[%.3g, %.3g] y=[%.3g, %.3g]" %
                      (cfg['zoom_x1'], x, cfg['zoom_y1'], y))
    else:
        print("Oops, should not have arrived here (in pickSomething).")
    return

def selectNode(canvas, cx, cy):
    objectList = list(canvas.find_overlapping(cx, cy, cx, cy))
    objectList.reverse()
    print("DEBUG objectList=", objectList)
    foundNodeIndices = []
    for obj in objectList:
        if canvas.type(obj) == "oval": foundNodeIndices.append(cfg['node_indx'][obj])
    print("DEBUG foundNodeIndices=", foundNodeIndices)
    if len(foundNodeIndices) > 1 and cfg['display_dialog_for_coincident_nodes']:
        print("[TODO] Put up a list of the possible nodes so that the user can select one.")
    try:
        # Just take the top node for now.
        cfg['selected_nodes'].append(foundNodeIndices[0])
        showStatusMsg("Selected nodes: %s" % cfg['selected_nodes'])
    except IndexError:
        pass
    return

def clearSelectedNodesList(event=None):
    # We need to accept an event for the key-binding.
    cfg['selected_nodes'] = []
    showStatusMsg("Selected nodes: %s" % cfg['selected_nodes'])
    return

def startZoom():
    cfg['pick_action'] = PickAction.corner_1
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
    if cfg['pick_action'] == PickAction.corner_1:
        c.delete("cursor1")
        c.create_line(cx, c_y_min, cx, c_y_max, fill='red', tags='cursor1')
        c.create_line(c_x_min, cy, c_x_max, cy, fill='red', tags='cursor1')
    if cfg['pick_action'] == PickAction.corner_2:
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

    [TODO]: choose a good tic increment and display labelled tic marks.
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

def plotWalls():
    c = cfg['canvas']
    c.delete('walls')
    for wall in kernel.walls:
        xmin = max(cfg['x_min'], wall.x_min)
        xmax = min(cfg['x_max'], wall.x_max)
        n = 100
        dx = (xmax - xmin)/n
        x = xmin
        x1 = canvas_x(x); y1 = canvas_y(wall(x))
        for i in range(n):
            xnew = x+dx
            x2 = canvas_x(xnew); y2 = canvas_y(wall(xnew))
            c.create_line(x1, y1, x2, y2, fill='blue', width=3, tags='walls')
            x, x1, y1 = xnew, x2, y2
    return

def plotMesh():
    """
    Plot the nodes as circles with line segments indicating
    the characteristic lines and streamlines.

    [TODO]: Clip the drawing to the selected x,y-ranges.
    """
    c = cfg['canvas']
    c.delete('nodes'); c.delete('nodeids');
    c.delete('mesh'); c.delete('stream')
    show_node_numbers = var_show_node_numbers.get() == "1"
    show_char_mesh = var_show_char_mesh.get() == "1"
    show_streamlines = var_show_streamlines.get() == "1"
    r = 15 if show_node_numbers else 3
    node_indx = {}
    object_id = {}
    if show_char_mesh:
        for i in kernel.char_mesh:
            node = kernel.nodes[i]
            x = canvas_x(node.x); y = canvas_y(node.y)
            myid = c.create_oval(x-r, y-r, x+r, y+r, outline='black',
                                 fill='gray', tags='nodes')
            object_id[node.indx] = myid
            node_indx[myid] = node.indx
            if show_node_numbers:
                c.create_text(x, y, text="%d" % node.indx, anchor='center', tags='nodeids')
            if node.cplus_up is not None:
                nb = kernel.nodes[node.cplus_up]
                xb = canvas_x(nb.x); yb = canvas_y(nb.y)
                line_seg = c.create_line(x, y, xb, yb, fill='green', width=2, tags='mesh')
                c.lower(line_seg)
            if node.cminus_up is not None:
                nb = kernel.nodes[node.cminus_up]
                xb = canvas_x(nb.x); yb = canvas_y(nb.y)
                line_seg = c.create_line(x, y, xb, yb, fill='green', width=2, tags='mesh')
                c.lower(line_seg)
            if node.cplus_down is not None:
                nb = kernel.nodes[node.cplus_down]
                xb = canvas_x(nb.x); yb = canvas_y(nb.y)
                line_seg = c.create_line(x, y, xb, yb, fill='green', width=2, tags='mesh')
                c.lower(line_seg)
            if node.cminus_down is not None:
                nb = kernel.nodes[node.cminus_down]
                xb = canvas_x(nb.x); yb = canvas_y(nb.y)
                line_seg = c.create_line(x, y, xb, yb, fill='green', width=2, tags='mesh')
                c.lower(line_seg)
            # It may be that the node is also part of a streamline, for example,
            # a free-boundary node.
            if node.czero_up is not None:
                nb = kernel.nodes[node.czero_up]
                xb = canvas_x(nb.x); yb = canvas_y(nb.y)
                line_seg = c.create_line(x, y, xb, yb, fill='brown', width=2, tags='stream')
                c.lower(line_seg)
            if node.czero_down is not None:
                nb = kernel.nodes[node.czero_down]
                xb = canvas_x(nb.x); yb = canvas_y(nb.y)
                line_seg = c.create_line(x, y, xb, yb, fill='brown', width=2, tags='stream')
                c.lower(line_seg)
    if show_streamlines:
        # Collect the node indices for each streamline starting point.
        streamline_nodes = []
        for i in kernel.streamlines:
            streamline_nodes += kernel.get_streamline_nodes(i)
        for i in streamline_nodes:
            node = kernel.nodes[i]
            x = canvas_x(node.x); y = canvas_y(node.y)
            myid = c.create_oval(x-r, y-r, x+r, y+r, outline='black',
                                 fill='gray', tags='nodes')
            object_id[node.indx] = myid
            node_indx[myid] = node.indx
            if show_node_numbers:
                c.create_text(x, y, text="%d" % node.indx, anchor='center', tags='nodeids')
            if node.czero_up is not None:
                nb = kernel.nodes[node.czero_up]
                xb = canvas_x(nb.x); yb = canvas_y(nb.y)
                line_seg = c.create_line(x, y, xb, yb, fill='brown', width=2, tags='stream')
                c.lower(line_seg)
            if node.czero_down is not None:
                nb = kernel.nodes[node.czero_down]
                xb = canvas_x(nb.x); yb = canvas_y(nb.y)
                line_seg = c.create_line(x, y, xb, yb, fill='brown', width=2, tags='stream')
                c.lower(line_seg)
    # We want to retain the connection between node indx and canvas object id
    # for use with selecting nodes.
    cfg['node_indx'] = node_indx
    cfg['object_id'] = object_id
    return

def refreshDisplay():
    plotAxes()
    plotWalls()
    plotMesh()
    return

def adjustScalesRefreshDisplay():
    cfg['equal_scales'] = var_equal_scales.get() == "1"
    set_xy_scales()
    refreshDisplay()
    return

def zoomToIncludeAll():
    x_min = 0.0; x_max = 0.0;
    y_min = 0.0; y_max = 0.0;
    for w in kernel.walls:
        x = w.x_min; x_min = min(x_min, x); x_max = max(x_max, x);
        y = w(x); y_min = min(y_min, y); y_max = max(y_max, y);
        x = w.x_max; x_min = min(x_min, x); x_max = max(x_max, x);
        y = w(x); y_min = min(y_min, y); y_max = max(y_max, y);
    for n in kernel.nodes:
        x = n.x; x_min = min(x_min, x); x_max = max(x_max, x);
        y = n.y; y_min = min(y_min, y); y_max = max(y_max, y);
    if x_max <= x_min: x_max = x_min + 1.0
    if y_max <= y_min: y_max = y_min + 1.0
    showStatusMsg("Zoom to include all x=[%.3g, %.3g] y=[%.3g, %.3g]" %
                  (x_min, x_max, y_min, y_max))
    set_xy_ranges(x_min, y_min, x_max, y_max)
    eraseCursors()
    refreshDisplay()
    return

def deleteSelectedNodesFromMesh():
    print("[TODO] complete deleteSelectedNodesFromMesh")
    return

def computeInteriorNode():
    print("[TODO] computeInteriorNode")
    return

def computeCminusWall0():
    print("[TODO] complete computeCminusWall0")
    return
