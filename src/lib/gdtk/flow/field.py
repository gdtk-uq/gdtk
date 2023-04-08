# field.py
"""
Field class for holding Eilmer field data for a FluidBlock.

PJ, 2023-04-08
"""

import zipfile
import json
import numpy

class Field():
    """
    Pick up a block of field-data from an Eilmer simulation.

    Once in memory, the data are available in a dictionary of numpy arrays.
    """
    __slots__ = ['ziparchive',
                 'variables', 'data', 'lengths',
                 'nic', 'njc', 'nkc', 'structured']

    def __init__(self, ziparchive):
        """
        Read a block's field data from one of Daryl's zip archives.

        ziparchive: file name of the zip archive.
        """
        if not zipfile.is_zipfile(ziparchive):
            raise RuntimeError(f"Not a zip arvhive: {ziparchive}")
        self.ziparchive = ziparchive
        zf = zipfile.ZipFile(ziparchive, mode='r')
        with zf.open('header.txt', 'r') as fp:
            text = fp.read()
        header = json.loads(text)
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
            with zf.open(dataFile, 'r') as fp:
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
        return

    def __repr__(self):
        str = f"Field(ziparchive={self.ziparchive}"
        str += f", variables={self.variables}"
        str += f", nic={self.nic}, njc={self.njc}, nkc={self.nkc}"
        str += ")"
        return str

