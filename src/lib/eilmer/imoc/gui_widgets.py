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

def init_widgets():
    """
    By the end of this function, all of the GUI widgets will be in place.
    """
    # Some of our variables need to be visible to other functions.
    global root
    root = Tk()
    root.title("Isentropic Method of Characteristics")
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
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
    c = Canvas(cf, scrollregion=(0, 0, 1000, 1200),
               yscrollcommand=v.set, xscrollcommand=h.set)
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
    menu_help = Menu(menubar)
    menubar.add_cascade(menu=menu_file, label='File')
    menu_file.add_command(label='Source file...', command=sourceFile)
    menu_file.add_command(label='Quit', command=quitProgram)
    menubar.add_cascade(menu=menu_help, label='Help')
    menu_help.add_command(label='About', command=aboutProgram)
    #
    # Some bindings.
    c.bind("<Button-1>", selectNode)
    c.bind("<Motion>", displayCoords)
    #
    # For the convenience of the user, leave the focus in the canvas widget
    c.focus()
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

def selectNode(event):
    # [TODO] convert coordinates to physical units
    status_msg.set("select node x=%d y=%d" % (event.x, event.y))
    return

def displayCoords(event):
    # [TODO] convert coordinates to physical units
    coord_x.set(event.x)
    coord_y.set(event.y)
    return
