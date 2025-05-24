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


DomainType = Enum('DomainType', ['FLUID', 'SOLID'])

GridType = Enum('GridType', ['STRUCTURED', 'UNSTRUCTURED'])

GridInfo = NamedTuple(
    'Gridinfo', [('tag', str),
                 ('type', GridType),
                 ('field_type', DomainType),
                 ('dimensions', int),
                 ('niv', int), ('njv', int), ('nkv', int),
                 ('nic', int), ('njc', int), ('nkc', int),
                 ('json', dict)])


BlockInfo = NamedTuple(
    'BlockInfo', [('id', int),
                  ('grid_type', GridType),
                  ('label', str),
                  ('ncells', int)])

class SimInfo:
    """
    A place to store the information about an lmr simulation.
    """
    __slots__ = ['lmr_cfg', 'sim_cfg',
                 'sim_dir', 'grid_dir', 'snaps_dir', 'vtk_dir',
                 'blocks', 'grids', 'snapshots', 'times',
                 'fluid_variables']

    def __init__(self, lmr_cfg):
        """
        We need to know some things about the lmr code configuration
        so that we can find the relevant files, etc.
        """
        self.lmr_cfg = lmr_cfg
        self.sim_dir = self.lmr_cfg["simulation-directory"]
        self.grid_dir = os.path.join(os.path.join(self.sim_dir, self.lmr_cfg["grid-directory"]))
        self.snaps_dir = os.path.join(self.sim_dir, self.lmr_cfg["snapshot-directory"])
        self.vtk_dir = os.path.join(self.sim_dir, self.lmr_cfg['vtk-output-directory'])
        #
        fname = os.path.join(self.sim_dir, self.lmr_cfg["config-filename"])
        with open(fname, 'r') as fp:
            self.sim_cfg = json.load(fp)
        #
        # Read the metadata for the list of blocks.
        #
        fname = os.path.join(self.sim_dir, self.lmr_cfg["block-list-filename"])
        self.blocks = []
        if os.path.exists(fname):
            with open(fname, 'r') as fp:
                content = fp.readlines()
                for line in content:
                    if line.strip().startswith('#'): continue
                    items = line.strip().split()
                    gt = GridType.STRUCTURED if items[1] == "structured_grid" else GridType.UNSTRUCTURED
                    bi = BlockInfo(id=int(items[0]), grid_type=gt, label=items[2], ncells=int(items[3]))
                    self.blocks.append(bi)
        # Check that the blocks are in the expected order.
        ids_in_order = [i == self.blocks[i].id for i in range(len(self.blocks))]
        if not all(ids_in_order):
            print("Warning: not all block ids are in order.")
        #
        # Read the metadata for the snapshots.
        # The primary source of information is the collection of directory names,
        # because both the steady-state and the transient code do the same thing
        # with these directories.  The transient simulation also writes metadata
        # for each snapshot that includes the time instants at which snapshots
        # are written.
        #
        self.snapshots = []
        names = os.listdir(self.snaps_dir)
        names = [n for n in names if re.match(r'\d\d\d\d', n)]
        names.sort()
        names_in_order = [i == int(names[i]) for i in range(len(names))]
        if not all(names_in_order):
            print("Warning: not all snapshots in order.")
        self.snapshots = names
        self.times = []
        if self.sim_cfg['solver_mode'] == 'transient':
            fname = os.path.join(self.snaps_dir, "snapshot-times-metadata")
            with open(fname, 'r') as fp:
                snapshots_metadata = yaml.safe_load(fp)
            snaps = list(snapshots_metadata.keys())
            snaps.sort()
            self.times = [float(snapshots_metadata[snap]['time']) for snap in snaps]
            found_snaps = [name in self.snapshots for name in snaps]
            if not all(found_snaps):
                print("Warning: not all snapshots listed in metadata are found.")
        #
        # List of names of the fluid variables.
        #
        self.fluid_variables = []
        fname = os.path.join(self.snaps_dir, "fluid"+self.lmr_cfg["metadata-extension"])
        with open(fname, 'r') as fp:
            fluid_metadata = yaml.safe_load(fp)
            self.fluid_variables = fluid_metadata["variables"]
        #
        # Grid metadata
        #
        self.grids = []
        fname = os.path.join(os.path.join(self.grid_dir, "grid"+self.lmr_cfg["metadata-extension"]))
        with open(fname, 'r') as fp:
            grid_metadata = json.load(fp)
        for i in range(grid_metadata['ngrids']):
            fname = ("grid-%04d" % i)+self.lmr_cfg["metadata-extension"]
            fname = os.path.join(os.path.join(self.grid_dir, fname))
            with open(fname, 'r') as fp:
                mdata = json.load(fp)
                gt = GridType.STRUCTURED if mdata['type'] == "structured_grid" else GridType.UNSTRUCTURED
                gi = GridInfo(tag=mdata['tag'], type=gt, field_type=mdata['fieldType'],
                              dimensions=mdata['dimensions'],
                              niv=mdata['niv'], njv=mdata['njv'], nkv=mdata['njv'],
                              nic=mdata['nic'], njc=mdata['njc'], nkc=mdata['nkc'],
                              json=mdata)
                self.grids.append(gi)
        #
        return

    def grid_metadata_as_dict(self):
        """Return the grid metadata a dict"""
        fname = os.path.join(os.path.join(self.grid_dir, "grid"+self.lmr_cfg["metadata-extension"]))
        with open(fname, 'r') as fp:
            grid_metadata = json.load(fp)
        return grid_metadata

    def read_grids(self):
        """
        Return a list containing all of the grids associated the simulation.
        """
        grids = []
        grid_format = self.sim_cfg["grid_format"]
        for i in range(len(self.grids)):
            grid_type = self.grids[i].type
            if grid_type == GridType.STRUCTURED:
                if grid_format == "gziptext":
                    fname = ("grid-%04d.gz" % i)
                    fname = os.path.join(os.path.join(self.grid_dir, fname))
                    grids.append(StructuredGrid(gzfile=fname))
                elif grid_format == "rawbinary":
                    fname = "grid-%04d.bin" % (i)
                    fname = os.path.join(os.path.join(self.grid_dir, fname))
                    grids.append(StructuredGrid(binaryfile=fname))
                else:
                    raise "Unknown format for grid[%d] %s" % (i, grid_format)
            elif grid_type == GridType.UNSTRUCTURED:
                raise "unstructured_grid reading not implemtented"
            else:
                raise "unknown grid type"
        return grids

    def read_snapshot(self, snap_name):
        """
        Read the solution data for a particular snapshot.
        """
        # [TODO] PJ 2025-05-23
        fields = []
        grids = []
        return fields, grids

    def load_pvd_into_pyvista(self, domain_type=DomainType.FLUID, snapshot=-1,
                              merged=False, as_point_data=False):
        """
        This loads the top-level PVD file for VTK output from Eilmer.

        The VTK output is ordinarily produced as the result of running 'snapshot2vtk'.

        Args:
            domain_type (DomainType) : indicates if domain contains fluid or solid data
            snapshot (int)           : select snapshot, default is last snapshot (-1 index into list)
            merged (bool)            : set to True to combine blocks; default is to leave partitioned
            as_point_data (bool)     : set to True to convert to point data; default is to leave as cell data

        Returns:
            A pyvista.DataSet object.
        """
        if not HAVE_PYVISTA:
            raise "pyvista module is not loaded"
        # locate the PVD file
        name = "fluid" if domain_type == DomainType.FLUID else "solid"
        fname = os.path.join(self.vtk_dir, name+".pvd")
        if not os.path.exists(fname):
            raise "Cannot find .pvd file: %s" % (fname,)
        reader = pv.get_reader(fname)
        reader.set_active_time_point(snapshot)
        pv_data = reader.read()
        if merged:
            pv_data = pv_data.combine()
        if as_point_data:
            pv_data = pv_data.cell_data_to_point_data()
        return pv_data

    # end of class SimInfo
