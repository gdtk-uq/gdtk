# lmr.py
# Tools for getting lmr configuration and simulation information.
#
# Peter J. 2025-05-14
#

import os
import json
import yaml
from typing import NamedTuple
import re
from enum import Enum
try:
    import pyvista as pv
    HAVE_PYVISTA = True
except ModuleNotFoundError:
    HAVE_PYVISTA = False

from gdtk.geom.sgrid import StructuredGrid


class DomainType(Enum):
    FLUID = 'fluid'
    SOLID = 'solid'

class LmrConfig:
    """
    A place to store the lmr code configuration.
    """
    __slots__ = ['data']

    def __init__(self, fname=None):
        """
        Reads the lmr config file from the specified file name
        or go looking in the installation directory.
        Input:
          fname: path/file name to the lmr configuration file
              If this is not supplied, the os environment is checked
              and a presumed path/file name is constructed.
        """
        if not fname:
            install_path = os.getenv('DGD')
            fname = os.path.join(install_path, 'etc', 'lmr.cfg')
        with open(fname, 'r') as fp:
            self.data = json.load(fp)
        return

    def __getitem__(self, indx):
        """
        The LmrConfig objects acts as a dictionary of config data.
        """
        return self.data[indx]


BlockInfo = NamedTuple("BlockInfo", [('id', int), ('type', str), ('label', str), ('ncells', int)])

class SimInfo:
    """
    A place to store the information about an lmr simulation.
    """
    __slots__ = ['lmr_cfg', 'sim_cfg', 'blocks', 'grids', 'grid_metadata',
                 'times', 'snapshots', 'fluid_variables']

    def __init__(self, lmr_cfg):
        """
        We need to know some things about the lmr code configuration
        so that we can find the relevant files, etc.
        """
        self.lmr_cfg = lmr_cfg
        sim_dir = self.lmr_cfg["simulation-directory"]
        fname = os.path.join(sim_dir, self.lmr_cfg["config-filename"])
        with open(fname, 'r') as fp:
            self.sim_cfg = json.load(fp)
        #
        # Read the metadata for the list of blocks.
        #
        fname = os.path.join(sim_dir, self.lmr_cfg["block-list-filename"])
        self.blocks = []
        if os.path.exists(fname):
            with open(fname, 'r') as fp:
                content = fp.readlines()
                for line in content:
                    if line.strip().startswith('#'): continue
                    items = line.strip().split()
                    self.blocks.append(BlockInfo(int(items[0]), items[1], items[2], int(items[3])))
        # Check that the blocks are in the expected order.
        ids_in_order = [i == self.blocks[i].id for i in range(len(self.blocks))]
        if not all(ids_in_order):
            print("Warning: not all block ids are in order.")
        #
        # Read the metadata for the time snapshots.
        #
        self.times = []
        if self.sim_cfg['solver_mode'] == 'transient':
            fname = os.path.join(sim_dir, self.lmr_cfg["times-filename"])
            if os.path.exists(fname):
                with open(fname, 'r') as fp:
                    content = fp.readlines()
                    # [TODO] Extract the data and put it into times_list.
        self.snapshots = []
        snaps_dirname = os.path.join(sim_dir, self.lmr_cfg["snapshot-directory"])
        names = os.listdir(snaps_dirname)
        names = [n for n in names if re.match(r'\d\d\d\d', n)]
        names.sort()
        names_in_order = [i == int(names[i]) for i in range(len(names))]
        if not all(names_in_order):
            print("Warning: not all snapshots in order.")
        self.snapshots = names
        #
        # List of names of the fluid variables.
        #
        self.fluid_variables = []
        fname = os.path.join(snaps_dirname, "fluid"+self.lmr_cfg["metadata-extension"])
        with open(fname, 'r') as fp:
            fluid_metadata = yaml.safe_load(fp)
            self.fluid_variables = fluid_metadata["variables"]
        #
        # Grid metadata
        #
        self.grids = []
        grid_dir = os.path.join(os.path.join(sim_dir, self.lmr_cfg["grid-directory"]))
        fname = os.path.join(os.path.join(grid_dir, "grid"+self.lmr_cfg["metadata-extension"]))
        with open(fname, 'r') as fp:
            self.grid_metadata = json.load(fp)
        ngrids = self.grid_metadata['ngrids']
        for i in range(ngrids):
            fname = ("grid-%04d" % i)+self.lmr_cfg["metadata-extension"]
            fname = os.path.join(os.path.join(grid_dir, fname))
            with open(fname, 'r') as fp:
                self.grids.append(json.load(fp))
        #
        return

    def read_grids(self):
        """
        Return a list containing all of the grids associated the simulation.
        """
        grids = []
        sim_dir = self.lmr_cfg["simulation-directory"]
        grid_dir = os.path.join(os.path.join(sim_dir, self.lmr_cfg["grid-directory"]))
        ngrids = self.grid_metadata['ngrids']
        grid_format = self.sim_cfg["grid_format"]
        for i in range(ngrids):
            grid_type = self.grids[i]["type"]
            if grid_type == "structured_grid":
                if grid_format == "gziptext":
                    fname = ("grid-%04d.gz" % i)
                    fname = os.path.join(os.path.join(grid_dir, fname))
                    grids.append(StructuredGrid(gzfile=fname))
                elif grid_format == "rawbinary":
                    fname = "grid-%04d.bin" % (i)
                    fname = os.path.join(os.path.join(grid_dir, fname))
                    grids.append(StructuredGrid(binaryfile=fname))
                else:
                    raise "Unknown format for grid[%d] %s" % (i, grid_format)
            elif grid_type == "unstructured_grid":
                raise "unstructured_grid reading not implemtented"
        return grids


# Service functions for picking up data from an Eilmer simulation

def load_pvd_into_pyvista(lmr_cfg, sim_info, domain_type=DomainType.FLUID, merged=False, as_point_data=False):
    """This loads the top-level PVD file for VTK output from Eilmer.

    The VTK output is ordinarily produced as the result of running 'snapshot2vtk'.

    Args:
        lmr_cfg (LmrConfig)      : config object for Eilmer program itself
        sim_info (SimInfo)       : information for specific simulation (in current working directory)
        domain_type (DomainType) : indicates if domain contains fluid or solid data
        merged (bool)            : set to True to combine blocks; default is to leave partitioned
        as_point_data (bool)     : set to True to convert to point data; default is to leave as cell data

    Returns:
        A pyvista.DataSet object.
    """
    if not HAVE_PYVISTA:
        raise "pyvista module is not loaded"
    # locate the PVD file
    vtk_dir = os.path.join(lmr_cfg['simulation-directory'], lmr_cfg['vtk-output-directory'])
    fname = os.path.join(vtk_dir, domain_type.value+".pvd")
    # now read
    reader = pv.get_reader(fname)
    pv_data = reader.read()
    if merged:
        pv_data = pv_data.combine()
    if as_point_data:
        pv_data = pv_data.cell_data_to_point_data()
    return pv_data
