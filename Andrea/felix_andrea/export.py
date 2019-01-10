import cPickle
import numpy.lib.format as npfmt
import os

class GrowingArrayFile(object):

    def __init__(self, filename, elementShape, dtype=float):
        self.fp = open(filename, 'wb')
        self.elementShape = elementShape
        self.numElements = 0
        self.dtype = dtype
        self.__writeHeader()


    def __writeHeader(self):
        self.fp.write(npfmt.magic(1, 0))
        npfmt.write_array_header_1_0(self.fp, {
            'shape': (self.numElements,) + self.elementShape,
            'fortran_order': False,
            'descr': npfmt.dtype_to_descr(self.dtype) })


    def __updateHeader(self):
        self.fp.seek(0)
        self.__writeHeader()
        self.fp.seek(0, os.SEEK_END)


    def write(self, array):
        assert array.shape[1:] == self.elementShape
        assert array.dtype == self.dtype

        # See for reference:
        # https://github.com/numpy/numpy/blob/maintenance/1.7.x/numpy/lib/format.py#L398
        if array.dtype.hasobject:
            # We contain Python objects so we cannot write out the data
            # directly. Instead, we will pickle it out with version 2 of the
            # pickle protocol.
            cPickle.dump(array, self.fp, protocol=2)
        else:
            array.tofile(self.fp)

        self.numElements += 1


    def __del__(self):
        self.__updateHeader()