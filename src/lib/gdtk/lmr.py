# lmr.py
# Tools for getting lmr configuration and simulation information.
#
# Peter J. 2025-05-14, 2025-06-13
#

import os
import subprocess
import json
import yaml
from typing import NamedTuple
import re
from enum import Enum
import gzip
import numpy as np
try:
    import pyvista as pv
    HAVE_PYVISTA = True
except ModuleNotFoundError:
    HAVE_PYVISTA = False
try:
    import pandas as pd
    HAVE_PANDAS = True
except ModuleNotFoundError:
    HAVE_PANDAS = False

from gdtk.geom.sgrid import StructuredGrid
from gdtk.gas import GasModel


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

# We will use Enums and NamedTuples to label the available metadata
# and make it conveniently findable in the Python module.

Face = Enum('Face', ['west', 'east', 'south', 'north', 'bottom', 'top'])
FaceDict = {'west':0, 'east':1, 'south':2, 'north':3, 'bottom':4, 'top':5}

DomainType = Enum('DomainType', ['FLUID', 'SOLID'])

GridType = Enum('GridType', ['STRUCTURED', 'UNSTRUCTURED'])

GridInfo = NamedTuple(
    'GridInfo',
    [('tag', str),
     ('type', GridType),
     ('field_type', DomainType),
     ('dimensions', int),
     ('niv', int), ('njv', int), ('nkv', int),
     ('nic', int), ('njc', int), ('nkc', int),
     # ('json', dict) # The entire JSON dict could be made available.
     ])

GridArrayInfo = NamedTuple(
    'GridArrayInfo',
    [('tag', str),
     ('shock_fitting', bool),
     ('idarray', list),
     ('idflatlist', list),
     ('nib', int), ('njb', int), ('nkb', int),
     ('niv', int), ('njv', int), ('nkv', int),
     ('nics', list), ('njcs', list), ('nkcs', list),
     # ('json', dict)
     ])

BlockInfo = NamedTuple(
    'BlockInfo',
    [('id', int),
     ('grid_type', GridType),
     ('label', str),
     ('ncells', int)
     ])

class SimInfo:
    """
    A place to store the information about an lmr simulation.
    """
    __slots__ = ['lmr_cfg', 'sim_cfg',
                 'sim_dir', 'grid_dir', 'snaps_dir', 'vtk_dir', 'loads_dir',
                 'blocks', 'grids', 'gridarrays', 'moving_grid',
                 'snapshots', 'times', 'loads_indices',
                 'fluid_variables', 'gas_model']

    def __init__(self, lmr_cfg):
        """
        We need to know some things about the lmr code configuration
        so that we can find the relevant files, etc.
        """
        self.lmr_cfg = lmr_cfg
        self.sim_dir = self.lmr_cfg["simulation-directory"]
        self.grid_dir = os.path.join(self.sim_dir, self.lmr_cfg["grid-directory"])
        self.snaps_dir = os.path.join(self.sim_dir, self.lmr_cfg["snapshot-directory"])
        self.vtk_dir = os.path.join(self.sim_dir, self.lmr_cfg['vtk-output-directory'])
        self.loads_dir = os.path.join(self.sim_dir, self.lmr_cfg['loads-directory'])
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
        self.moving_grid = self.sim_cfg['grid_motion'] != "none"
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
        self.gridarrays = []
        fname = os.path.join(os.path.join(self.grid_dir, "grid"+self.lmr_cfg["metadata-extension"]))
        with open(fname, 'r') as fp:
            grid_metadata = json.load(fp)
        for i in range(grid_metadata['ngridarrays']):
            mdata = grid_metadata['gridarrays'][i]
            gai = GridArrayInfo(
                tag=mdata['tag'], shock_fitting=mdata['shock_fitting'],
                idarray=mdata['idarray'], idflatlist=mdata['idflatlist'],
                nib=mdata['nib'], njb=mdata['njb'], nkb=mdata['nkb'],
                niv=mdata['niv'], njv=mdata['njv'], nkv=mdata['nkv'],
                nics=mdata['nics'], njcs=mdata['njcs'], nkcs=mdata['nkcs'],
                # json=mdata
            )
            self.gridarrays.append(gai)
        for i in range(grid_metadata['ngrids']):
            fname = ("grid-%04d" % i)+self.lmr_cfg["metadata-extension"]
            fname = os.path.join(os.path.join(self.grid_dir, fname))
            with open(fname, 'r') as fp:
                mdata = json.load(fp)
                gt = GridType.STRUCTURED if mdata['type'] == "structured_grid" else GridType.UNSTRUCTURED
                gi = GridInfo(
                    tag=mdata['tag'], type=gt, field_type=mdata['fieldType'],
                    dimensions=mdata['dimensions'],
                    niv=mdata['niv'], njv=mdata['njv'], nkv=mdata['njv'],
                    nic=mdata['nic'], njc=mdata['njc'], nkc=mdata['nkc'],
                    # json=mdata
                )
                self.grids.append(gi)
        #
        # Loads that can be found.
        self.loads_indices = []
        names = os.listdir(self.loads_dir)
        indices = [int(n) for n in names if re.match(r'\d\d\d\d', n)]
        indices.sort()
        self.loads_indices = indices
        #
        # Finally, set up a GasModel object that may be handy.
        #
        self.gas_model = GasModel(self.sim_cfg['gas_model_file'])
        return

    def grid_metadata_as_dict(self):
        """
        Return the grid metadata a dict

        for use when the metadata writen by lmr is known to have the needed information
        but this module does not already provide it in a convenient form.
        """
        fname = os.path.join(os.path.join(self.grid_dir, "grid"+self.lmr_cfg["metadata-extension"]))
        with open(fname, 'r') as fp:
            grid_metadata = json.load(fp)
        return grid_metadata

    def read_grids(self, grid_dir=None):
        """
        Return a list containing all of the grids associated the simulation.

        By default, look in the grid_dir associated with the preparation of the simulation.
        We can override this to pick the '0000' snapshot or any later snapshot
        for a transient simulation with moving-grid.
        """
        if not grid_dir: grid_dir = self.grid_dir
        grids = []
        grid_format = self.sim_cfg["grid_format"]
        for i in range(len(self.grids)):
            grid_type = self.grids[i].type
            if grid_type == GridType.STRUCTURED:
                if grid_format == "gziptext":
                    fname = ("grid-%04d.gz" % i)
                    fname = os.path.join(os.path.join(grid_dir, fname))
                    grids.append(StructuredGrid(gzfile=fname))
                elif grid_format == "rawbinary":
                    fname = "grid-%04d.bin" % (i)
                    fname = os.path.join(os.path.join(grid_dir, fname))
                    grids.append(StructuredGrid(binaryfile=fname))
                else:
                    raise RuntimeError(f"Unknown format for grid[{i}] {grid_format}")
            elif grid_type == GridType.UNSTRUCTURED:
                raise RuntimeError("unstructured_grid reading not implemtented")
            else:
                raise RuntimeError("unknown grid type")
        return grids

    def read_snapshot(self, snap_name):
        """
        Read the solution data for a particular snapshot.
        """
        return Snapshot(self, snap_name)

    def load_pvd_into_pyvista(self, domain_type=DomainType.FLUID, snapshot=-1,
                              merged=False, as_point_data=False, make_vtk_files=True):
        """
        Returns a pyvista.DataSet object containing the simulated flowfield.

        This is done by loading the top-level PVD file for VTK output from Eilmer
        which is ordinarily produced as the result of running 'snapshot2vtk'.

        Args:
        domain_type (DomainType) : indicates if domain contains fluid or solid data
        snapshot (int)           : select snapshot, default is last snapshot (-1 index into list)
        merged (bool)            : set to True to combine blocks; default is to leave partitioned
        as_point_data (bool)     : set to True to convert to point data; default is to leave as cell data
        make_vtk_files (bool)    : make the VTK files if they are not already present
        """
        if not HAVE_PYVISTA:
            raise RuntimeError("pyvista module is not loaded")
        # locate the PVD file
        name = "fluid" if domain_type == DomainType.FLUID else "solid"
        fname = os.path.join(self.vtk_dir, name+".pvd")
        if not os.path.exists(fname):
            if make_vtk_files:
                proc = subprocess.run(['lmr', 'snapshot2vtk', '--all'], capture_output=True)
                if not proc.returncode == 0:
                    raise RuntimeError(f"Error when running lmr snapshot2vtk: {proc.stdout}, {proc.stderr}")
            else:
                raise RuntimeError(f"Cannot find .pvd file: {fname}")
        #
        reader = pv.get_reader(fname)
        reader.set_active_time_point(snapshot)
        pv_data = reader.read()
        if merged:
            pv_data = pv_data.combine()
        if as_point_data:
            pv_data = pv_data.cell_data_to_point_data()
        return pv_data

    def read_loads(self, indx, group, bf_list=None):
        """
        Returns a Pandas DataFrame with the surface loads for the specified group.

        Args:
        indx: int index of the loads files to be read
        group: string name of the group associated with the required loads
        bf_list: list of tuples that specify, in order, particular block+face pairs to load

        The order in which the loads files are appended to the DateFrame
        is determined by the directory listing, so it may make sense to
        sort the resulting DataFrame on one of the columns, say 'pos.x'.
        """
        if not HAVE_PANDAS:
            raise RuntimeError("pandas module is not loaded")
        d = os.path.join(self.loads_dir, ('%04d' % indx))
        if bf_list:
            # We have been given a list of block+face tuples.
            names = []
            for blk,face in bf_list:
                if isinstance(face, str): face = FaceDict[face.lower()]
                names.append('blk-%04d-bndry-%d-%s.dat' % (blk, face, group))
        else:
            # Locate the files by name, filtering for the group name.
            names = os.listdir(d)
            target_str = '-' + group + '.'
            names = [n for n in names if n.find(target_str) > 0]
        blk_dfs = []
        for n in names:
            f = os.path.join(d, n)
            blk_dfs.append(pd.read_table(f, sep=r'\s+'))
        df = pd.concat(blk_dfs)
        return df

    # end of class SimInfo


class Snapshot:
    """
    A place to store the field data and a set of grids.

    And some functions to provide convenient access.
    """
    __slots__ = ['sim', 'fields', 'grids']

    def __init__(self, sim, snap_name):
        """
        We pick up a specific snapshot for a simulation.
        """
        self.sim = sim
        snap_dir = os.path.join(sim.snaps_dir, snap_name)
        if not os.path.exists(snap_dir):
            raise RuntimeError(f"Could not find snapshot {snap_dir}")
        grid_dir = snap_dir if sim.moving_grid else os.path.join(sim.snaps_dir, sim.snapshots[0])
        self.grids = sim.read_grids(grid_dir)
        n = len(sim.grids)
        assert n == len(self.grids), f"Expected number of grids ({n}) and actually loaded grids ({len(grids)})."
        field_format = sim.sim_cfg["field_format"]
        binary_data = field_format == "rawbinary"
        self.fields = []
        for i in range(n):
            grid_meta = sim.grids[i]
            fname = grid_meta.field_type + ("-%04d" % (i)) + (".bin" if binary_data else ".gz")
            fname = os.path.join(os.path.join(snap_dir, fname))
            if not os.path.exists(fname):
                raise RuntimeError(f"Could not find field file {fname}")
            f = open(fname, 'rb') if binary_data else gzip.open(fname, 'rt')
            combined_data = np.frombuffer(f.read(), dtype=float) if binary_data else np.loadtxt(f)
            f.close()
            ntotal = combined_data.shape[0]
            # [FIX-ME] Need to handle solid blocks appropriately
            nvars = len(sim.fluid_variables)
            ncells = ntotal // nvars
            combined_data = combined_data.reshape((nvars,ncells))
            field_data = {}
            for j,var in enumerate(sim.fluid_variables):
                vdata = combined_data[j,:]
                if grid_meta.type == GridType.STRUCTURED:
                    nic = grid_meta.nic; njc = grid_meta.njc; nkc = grid_meta.nkc
                    if grid_meta.dimensions == 2:
                        vdata = vdata.reshape((njc,nic)).transpose()
                    else:
                        vdata = vdata.reshape((nkc,njc,nic)).transpose()
                field_data[var] = vdata
            self.fields.append(field_data)
        return

    def get_slice(self, var, blocks=None, i=None, j=None, k=None):
        """
        Returns four arrays of data for x, y, z and var for all cells along one index.

        var    : str, name of particular field variable
        blocks : int, index, or range of indices, given in a list
        i      : int or None
        j      : int or None
        k      : int or None

        In 2D, we supply either a specific i-index or j-index and
        select all cells along the index that is None.
        In 3D, two of the three need to be integers and all cells
        will be selected along the final index direction.
        """
        # [TODO] Check for valid arguments
        #
        dims = self.sim.grids[0].dimensions
        x = np.array([], dtype=float)
        y = np.array([], dtype=float)
        if dims == 3: z = np.array([], dtype=float)
        values = np.array([], dtype=float)
        #
        # Make sure that we have a list of block ids.
        typ = type(blocks)
        if typ == type([1]):
            pass
        elif typ == type(1) or typ == type(range(1)):
            blocks = list(blocks)
        elif typ == type(None):
            blocks = list(range(len(self.fields)))
        else:
            raise RuntimeError("Don't know how to handle blocks="+str(blocks))
        #
        # We select a single line of cells along one index direction,
        # the one that is None.
        #
        for ib in blocks:
            if self.sim.grids[ib].type == GridType.UNSTRUCTURED: continue
            # The following works for a structured grid only.
            if dims == 3:
                if i != None and j != None:
                    x = np.concat((x, self.fields[ib]['pos.x'][i,j,:]), axis=None)
                    y = np.concat((y, self.fields[ib]['pos.y'][i,j,:]), axis=None)
                    z = np.concat((z, self.fields[ib]['pos.z'][i,j,:]), axis=None)
                    values = np.concat((values, self.fields[ib][var][i,j,:]), axis=None)
                elif j != None and k != None:
                    x = np.concat((x, self.fields[ib]['pos.x'][:,j,k]), axis=None)
                    y = np.concat((y, self.fields[ib]['pos.y'][:,j,k]), axis=None)
                    z = np.concat((z, self.fields[ib]['pos.z'][:,j,k]), axis=None)
                    values = np.concat((values, self.fields[ib][var][:,j,k]), axis=None)
                elif i != None and k != None:
                    x = np.concat((x, self.fields[ib]['pos.x'][i,:,k]), axis=None)
                    y = np.concat((y, self.fields[ib]['pos.y'][i,:,k]), axis=None)
                    z = np.concat((z, self.fields[ib]['pos.z'][i,:,k]), axis=None)
                    values = np.concat((values, self.fields[ib][var][i,:,k]), axis=None)
                else:
                    raise RuntimeError("Did not correctly select i, j or k index.")
            else:
                # Presumably 2D
                if i != None:
                    x = np.concat((x, self.fields[ib]['pos.x'][i,:]), axis=None)
                    y = np.concat((y, self.fields[ib]['pos.y'][i,:]), axis=None)
                    values = np.concat((values, self.fields[ib][var][i,:]), axis=None)
                elif j != None:
                    x = np.concat((x, self.fields[ib]['pos.x'][:,j]), axis=None)
                    y = np.concat((y, self.fields[ib]['pos.y'][:,j]), axis=None)
                    values = np.concat((values, self.fields[ib][var][:,j]), axis=None)
                else:
                    raise RuntimeError("Did not correctly select i or j index.")
        if dims == 3: return x, y, z, values
        return x, y, values

    def add_vars(self, names):
        """
        Add a derived variable to the field dictionary for all blocks in the Snapshot.

        names: list(str) names may contain mach, pitot, total_p, total_h, total_T
        """
        current_vars = list(self.fields[0].keys())
        names = [name.lower() for name in names]
        for field in self.fields:
            p = field['p']
            rho = field['rho']
            a = field['a']
            vx = field['vel.x']; vy = field['vel.y']
            if 'vel.z' in current_vars:
                vz = field['vel.z']
                v = np.sqrt(vx*vx + vy*vy + vz*vz)
            else:
                v = np.sqrt(vx*vx + vy*vy)
            M = v/a
            g = a*a*rho/p # approximation for gamma
            #
            if 'mach' in names:
                field['mach'] = M
            if 'pitot' in names:
                # Go through the shock and isentropic compression.
                t1 = (g+1)*M*M/2;
                t2 = (g+1)/(2*g*M*M - (g-1));
                pitot_p_supersonic = p * np.pow(t1,(g/(g-1))) * np.pow(t2,(1/(g-1)))
                # Isentropic compression only.
                t1 = 1 + 0.5*(g-1)*M*M
                pitot_p_subsonic = p * np.pow(t1,(g/(g-1)))
                # Select one to save.
                field['pitot'] = np.where(M > 1, pitot_p_supersonic, pitot_p_subsonic)
            if 'total_p' in names:
                # Isentropic process only.
                t1 = 1 + 0.5*(g-1)*M*M;
                field['total_p'] = p * np.pow(t1,(g/(g-1)))
            if 'total_h' in names:
                emodes_names = [var for var in current_vars if var.find('e_modes') >= 0]
                e_int = np.zeros(p.shape, dtype=float)
                for var in emodes_names: e_int += field[var]
                total_h = p/rho + field['e'] + 0.5*v*v
                if 'tke' in current_vars: total_h += field['tke']
                field['total_h'] = total_h
            if 'total_T' in names:
                field['total_T'] = T * (1.0 + 0.5*(g - 1.0)*M*M)
        return
