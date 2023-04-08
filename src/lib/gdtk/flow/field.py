# field.py
"""
Field class for holding Eilmer field data for a FluidBlock.

PJ, 2023-04-08
"""
import os
import re
import zipfile
import gzip
import json
import numpy

class Field():
    """
    Pick up a block of flow-field data from an Eilmer simulation.

    Once in memory, the data are available in a dictionary of numpy arrays.
    """
    __slots__ = ['new_format', 'ziparchive', 'gzfile',
                 'sim_time', 'variables', 'data', 'lengths',
                 'nic', 'njc', 'nkc', 'structured']

    def __init__(self, ziparchive=None, gzfile=None):
        """
        Read a block's field data from either a zip-archives or a gzipped file.

        ziparchive: file name of Daryl's zip archive.
        gzfile: file name of the "classic" Eilmer flow file.
        """
        if ziparchive and not gzfile:
            self.new_format = True
        elif not ziparchive and gzfile:
            self.new_format = False
        else:
            raise RuntimeError("Request either a ziparchive or a gzfile as input.")
        #
        if self.new_format:
            # Daryl's zip archive format.
            if not zipfile.is_zipfile(ziparchive):
                raise RuntimeError(f"Not a zip archive: {ziparchive}")
            self.ziparchive = ziparchive
            za = zipfile.ZipFile(ziparchive, mode='r')
            text = za.read('header.txt')
            header = json.loads(text)
            self.sim_time = header['sim_time']
            self.variables = list(header['variables'].keys())
            self.nic = header['nic']
            self.njc = header['njc']
            self.nkc = header['nkc']
            self.data = {}
            self.lengths = {} # a length of 1 indicates a scalar quantity
            for var in self.variables:
                metaData = header['variables'][var]
                dataFile = metaData['data']
                n_cells = metaData['dimension'][0]
                length = metaData['dimension'][1]
                with za.open(dataFile, 'r') as fp:
                    dataBlob = numpy.loadtxt(fp)
                if length == 1:
                    # Scalar quantity
                    if self.nkc > 1:
                        dataBlob = dataBlob.reshape(self.nkc, self.njc, self.nic)
                    else:
                        dataBlob = dataBlob.reshape(self.njc, self.nic)
                    self.data[var] = dataBlob.transpose() # index as [i,j] or [i,j,k]
                else:
                    # Array quantity
                    if self.nkc > 1:
                        dataBlob = dataBlob.reshape(length, self.nkc, self.njc, self.nic)
                    else:
                        dataBlob = dataBlob.reshape(length, self.njc, self.nic)
                    self.data[var] = dataBlob.transpose() # index as [i,j,n] or [i,j,k,n]
                self.lengths[var] = length
        else:
            # Classic (gzipped text) format for Eilmer flow file.
            if not os.path.exists(gzfile):
                raise RuntimeError(f"Could not find gzfile: {gzfile}")
            self.gzfile = gzfile
            fp = gzip.open(gzfile, 'rt')
            text = fp.readline()
            if text.strip() != "structured_grid_flow 1.0":
                raise RuntimeError(f"Did not see correct first line: {text}")
            text = fp.readline() # label
            text = fp.readline() # sim_time
            items = re.split(":", text)
            self.sim_time = float(items[1])
            text = fp.readline() # variables:
            items = re.split(":", text)
            nvariables = int(items[1])
            print(f"nvariables={nvariables}")
            text = fp.readline().strip() # The variable names
            # Assume no spaces within each name.
            items = text.split(' ')
            self.variables = [item.replace('"', '') for item in items]
            text = fp.readline() # dimensions
            items = re.split(":", text)
            dimensions = int(items[1])
            text = fp.readline() # nicell
            items = re.split(":", text)
            self.nic = int(items[1])
            text = fp.readline() # njcell
            items = re.split(":", text)
            self.njc = int(items[1])
            text = fp.readline() # nkcell
            items = re.split(":", text)
            self.nkc = int(items[1])
            self.data = {}
            if dimensions==2:
                if not self.nkc==1:
                    raise RuntimeError(f"For a 2D grid, expect nkc==1 but got {self.nkc}")
                for var in self.variables:
                    self.data[var] = numpy.zeros((self.nic,self.njc), float)
                for j in range(self.njc):
                    for i in range(self.nic):
                        text = fp.readline().strip()
                        items = text.split(' ')
                        assert nvariables == len(items), "Wrong number of items on data line."
                        for n in range(nvariables):
                            self.data[self.variables[n]][i,j] = float(items[n])
            else:
                # 3D grid
                for var in self.variables:
                    self.data[var] = numpy.zeros((self.nic,self.njc,self.nkc), float)
                for k in range(self.nkc):
                    for j in range(self.njc):
                        for i in range(self.nic):
                            text = fp.readline().strip()
                            items = text.split(' ')
                            assert nvariables == len(items), "Wrong number of items on data line."
                            for n in range(nvariables):
                                self.data[self.variables[n]][i,j,k] = float(items[n])
            fp.close()
        return

    def __repr__(self):
        str = "Field("
        str += f"ziparchive={self.ziparchive}" if self.new_format else f"gzfile={self.gzfile}"
        str += f", sim_time={self.sim_time}"
        str += f", variables={self.variables}"
        str += f", nic={self.nic}, njc={self.njc}, nkc={self.nkc}"
        str += ")"
        return str

