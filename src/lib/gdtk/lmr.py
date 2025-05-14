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
    __slots__ = ['lmr_cfg', 'sim_cfg', 'blocks', 'grids',
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
        names = [n for n in names if re.match('\d\d\d\d', n)]
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
            grid_metadata = json.load(fp)
        print("grid_metadata=", grid_metadata)
        ngrids = grid_metadata['ngrids']
        for i in range(ngrids):
            fname = ("grid-%04d" % i)+self.lmr_cfg["metadata-extension"]
            fname = os.path.join(os.path.join(grid_dir, fname))
            with open(fname, 'r') as fp:
                self.grids.append(json.load(fp))
        #
        return
