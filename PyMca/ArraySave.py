#/*##########################################################################
# Copyright (C) 2004-2012 European Synchrotron Radiation Facility
#
# This file is part of the PyMca X-ray Fluorescence Toolkit developed at
# the ESRF by the Software group.
#
# This file is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# This file is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
#############################################################################*/
__author__ = "V.A. Sole - ESRF Data Analysis"
import os
import numpy
import time

try:
    from PyMca import EdfFile
    from PyMca import TiffIO
except ImportError:
    print("ArraySave.py is importing EdfFile and TiffIO from local directory")
    import EdfFile
    import TiffIO

HDF5 = True
try:
    import h5py
except ImportError:
    HDF5 = False


DEBUG = 0


def getDate():
    localtime = time.localtime()
    gtime = time.gmtime()
    #year, month, day, hour, minute, second,\
    #      week_day, year_day, delta = time.localtime()
    year = localtime[0]
    month = localtime[1]
    day = localtime[2]
    hour = localtime[3]
    minute = localtime[4]
    second = localtime[5]
    #get the difference against Greenwich
    delta = hour - gtime[3]
    return "%4d-%02d-%02dT%02d:%02d:%02d%+02d:00" % (year, month, day, hour,
                                                     minute, second, delta)


def save2DArrayListAsASCII(datalist, filename,
                           labels=None, csv=False, csvseparator=";"):
    if type(datalist) != type([]):
        datalist = [datalist]
    r, c = datalist[0].shape
    ndata = len(datalist)
    if os.path.exists(filename):
        try:
            os.remove(filename)
        except OSError:
            pass
    if labels is None:
        labels = []
        for i in range(len(datalist)):
            labels.append("Array_%d" % i)
    if len(labels) != len(datalist):
        raise ValueError("Incorrect number of labels")
    if csv:
        header = '"row"%s"column"' % csvseparator
        for label in labels:
            header += '%s"%s"' % (csvseparator, label)
    else:
        header = "row  column"
        for label in labels:
            header += "  %s" % label
    filehandle = open(filename, 'w+')
    filehandle.write('%s\n' % header)
    fileline = ""
    if csv:
        for row in range(r):
            for col in range(c):
                fileline += "%d" % row
                fileline += "%s%d" % (csvseparator, col)
                for i in range(ndata):
                    fileline += "%s%g" % (csvseparator, datalist[i][row, col])
                fileline += "\n"
                filehandle.write("%s" % fileline)
                fileline = ""
    else:
        for row in range(r):
            for col in range(c):
                fileline += "%d" % row
                fileline += "  %d" % col
                for i in range(ndata):
                    fileline += "  %g" % datalist[i][row, col]
                fileline += "\n"
                filehandle.write("%s" % fileline)
                fileline = ""
    filehandle.write("\n")
    filehandle.close()


def save2DArrayListAsEDF(datalist, filename, labels=None, dtype=None):
    if type(datalist) != type([]):
        datalist = [datalist]
    ndata = len(datalist)
    if os.path.exists(filename):
        try:
            os.remove(filename)
        except OSError:
            pass
    if labels is None:
        labels = []
        for i in range(ndata):
            labels.append("Array_%d" % i)
    if len(labels) != ndata:
        raise ValueError("Incorrect number of labels")
    edfout = EdfFile.EdfFile(filename, access="ab")
    for i in range(ndata):
        if dtype is None:
            edfout.WriteImage({'Title': labels[i]},
                                datalist[i], Append=1)
        else:
            edfout.WriteImage({'Title': labels[i]},
                               datalist[i].astype(dtype),
                               Append=1)
    del edfout  # force file close


def save2DArrayListAsMonochromaticTiff(datalist, filename,
                                       labels=None, dtype=None):
    if type(datalist) != type([]):
        datalist = [datalist]
    ndata = len(datalist)
    if dtype is None:
        dtype = datalist[0].dtype
        for i in range(len(datalist)):
            dtypeI = datalist[i].dtype
            if dtypeI in [numpy.float32, numpy.float64] or\
               dtypeI.str[-2] == 'f':
                dtype = numpy.float32
                break
            elif dtypeI != dtype:
                dtype = numpy.float32
                break
    if os.path.exists(filename):
        try:
            os.remove(filename)
        except OSError:
            pass
    if labels is None:
        labels = []
        for i in range(ndata):
            labels.append("Array_%d" % i)
    if len(labels) != ndata:
        raise ValueError("Incorrect number of labels")
    outfileInstance = TiffIO.TiffIO(filename, mode="wb+")
    for i in range(ndata):
        if i == 1:
            outfileInstance = TiffIO.TiffIO(filename, mode="rb+")
        if dtype is None:
            data = datalist[i]
        else:
            data = datalist[i].astype(dtype)
        outfileInstance.writeImage(data, info={'Title': labels[i]})
    outfileInstance.close()  # force file close


def openHDF5File(name, mode='a', **kwargs):
    """
    Open an HDF5 file.

    Valid modes (like Python's file() modes) are:
    - r   Readonly, file must exist
    - r+  Read/write, file must exist
    - w   Create file, truncate if exists
    - w-  Create file, fail if exists
    - a   Read/write if exists, create otherwise (default)

    sorted_with is a callable function like python's builtin sorted, or
    None.
    """

    h5file = h5py.File(name, mode, **kwargs)
    if h5file.mode != 'r' and len(h5file) == 0:
        if 'file_name' not in h5file.attrs:
            attr = 'file_name'
            txt = "%s" % name
            dtype = '<S%d' % len(txt)
            h5file.attrs.create(attr, txt, dtype=dtype)
        if 'file_time' not in h5file.attrs:
            attr = 'file_time'
            txt = "%s" % getDate()
            dtype = '<S%d' % len(txt)
            h5file.attrs.create(attr, txt, dtype=dtype)
        if 'HDF5_version' not in h5file.attrs:
            attr = 'HDF5_version'
            txt = "%s" % h5py.version.hdf5_version
            dtype = '<S%d' % len(txt)
            h5file.attrs.create(attr, txt, dtype=dtype)
        if 'HDF5_API_version' not in h5file.attrs:
            attr = 'HDF5_API_version'
            txt = "%s" % h5py.version.api_version
            dtype = '<S%d' % len(txt)
            h5file.attrs.create(attr, txt, dtype=dtype)
        if 'h5py_version' not in h5file.attrs:
            attr = 'h5py_version'
            txt = "%s" % h5py.version.version
            dtype = '<S%d' % len(txt)
            h5file.attrs.create(attr, txt, dtype=dtype)
        if 'creator' not in h5file.attrs:
            attr = 'creator'
            txt = "%s" % 'PyMca'
            dtype = '<S%d' % len(txt)
            h5file.attrs.create(attr, txt, dtype=dtype)
        #if 'format_version' not in self.attrs and len(h5file) == 0:
        #    h5file.attrs['format_version'] = __format_version__

    return h5file


def getHDF5FileInstanceAndBuffer(filename, shape,
                                 buffername="data",
                                 dtype=numpy.float32,
                                 interpretation=None,
                                 compression=None):
    if not HDF5:
        raise IOError('h5py does not seem to be installed in your system')

    if os.path.exists(filename):
        try:
            os.remove(filename)
        except:
            raise IOError("Cannot overwrite existing file!")
    hdf = openHDF5File(filename, 'a')
    entryName = "data"

    #entry
    nxEntry = hdf.require_group(entryName)
    if 'NX_class' not in nxEntry.attrs:
        nxEntry.attrs['NX_class'] = 'NXentry'.encode('utf-8')
    elif nxEntry.attrs['NX_class'] != 'NXentry'.encode('utf-8'):
        #should I raise an error?
        pass
    nxEntry['title'] = "PyMca saved 3D Array".encode('utf-8')
    nxEntry['start_time'] = getDate().encode('utf-8')
    nxData = nxEntry.require_group('NXdata')
    if 'NX_class' not in nxData.attrs:
        nxData.attrs['NX_class'] = 'NXdata'.encode('utf-8')
    elif nxData.attrs['NX_class'] == 'NXdata'.encode('utf-8'):
        #should I raise an error?
        pass
    if compression:
        if DEBUG:
            print("Saving compressed and chunked dataset")
        chunk1 = int(shape[1] / 10)
        if chunk1 == 0:
            chunk1 = shape[1]
        for i in [11, 10, 8, 7, 5, 4]:
            if (shape[1] % i) == 0:
                chunk1 = int(shape[1] / i)
                break
        chunk2 = int(shape[2] / 10)
        if chunk2 == 0:
            chunk2 = shape[2]
        for i in [11, 10, 8, 7, 5, 4]:
            if (shape[2] % i) == 0:
                chunk2 = int(shape[2] / i)
                break
        data = nxData.require_dataset(buffername,
                           shape=shape,
                           dtype=dtype,
                           chunks=(1, chunk1, chunk2),
                           compression=compression)
    else:
        #no chunking
        if DEBUG:
            print("Saving not compressed and not chunked dataset")
        data = nxData.require_dataset(buffername,
                           shape=shape,
                           dtype=dtype,
                           compression=None)
    data.attrs['signal'] = numpy.int32(1)
    if interpretation is not None:
        data.attrs['interpretation'] = interpretation.encode('utf-8')
    for i in range(len(shape)):
        dim = numpy.arange(shape[i]).astype(numpy.float32)
        dset = nxData.require_dataset('dim_%d' % i,
                               dim.shape,
                               dim.dtype,
                               dim,
                               chunks=dim.shape)
        dset.attrs['axis'] = numpy.int32(i + 1)
    nxEntry['end_time'] = getDate().encode('utf-8')
    return hdf, data


def save3DArrayAsMonochromaticTiff(data, filename,
                                   labels=None, dtype=None, mcaindex=-1):
    ndata = data.shape[mcaindex]
    if dtype is None:
        dtype = numpy.float32
    if os.path.exists(filename):
        try:
            os.remove(filename)
        except OSError:
            pass
    if labels is None:
        labels = []
        for i in range(ndata):
            labels.append("Array_%d" % i)
    if len(labels) != ndata:
        raise ValueError("Incorrect number of labels")
    outfileInstance = TiffIO.TiffIO(filename, mode="wb+")
    if mcaindex in [2, -1]:
        for i in range(ndata):
            if i == 1:
                outfileInstance = TiffIO.TiffIO(filename, mode="rb+")
            if dtype is None:
                tmpData = data[:, :, i]
            else:
                tmpData = data[:, :, i].astype(dtype)
            outfileInstance.writeImage(tmpData, info={'Title': labels[i]})
            if (ndata > 10):
                print("Saved image %d of %d" % (i + 1, ndata))
    elif mcaindex == 1:
        for i in range(ndata):
            if i == 1:
                outfileInstance = TiffIO.TiffIO(filename, mode="rb+")
            if dtype is None:
                tmpData = data[:, i, :]
            else:
                tmpData = data[:, i, :].astype(dtype)
            outfileInstance.writeImage(tmpData, info={'Title': labels[i]})
            if (ndata > 10):
                print("Saved image %d of %d" % (i + 1, ndata))
    else:
        for i in range(ndata):
            if i == 1:
                outfileInstance = TiffIO.TiffIO(filename, mode="rb+")
            if dtype is None:
                tmpData = data[i]
            else:
                tmpData = data[i].astype(dtype)
            outfileInstance.writeImage(tmpData, info={'Title': labels[i]})
            if (ndata > 10):
                print("Saved image %d of %d" % (i + 1, ndata))
    outfileInstance.close()  # force file close

# it should be used to name the data that for the time being is named 'data'.
def save3DArrayAsHDF5(data, filename, axes=None, labels=None, dtype=None, mode='nexus',
                      mcaindex=-1, interpretation=None, compression=None):
    if not HDF5:
        raise IOError('h5py does not seem to be installed in your system')
    if (mcaindex == 0) and (interpretation in ["spectrum", None]):
        #stack of images to be saved as stack of spectra
        modify = True
        shape = [data.shape[1], data.shape[2], data.shape[0]]
    elif (mcaindex != 0) and (interpretation in ["image"]):
        #stack of spectra to be saved as stack of images
        modify = True
        shape = [data.shape[2], data.shape[0], data.shape[1]]
    else:
        modify = False
        shape = data.shape
    if dtype is None:
        dtype = data.dtype
    if mode.lower() in ['nexus', 'nexus+']:
        #raise IOError, 'NeXus data saving not implemented yet'
        if os.path.exists(filename):
            try:
                os.remove(filename)
            except:
                raise IOError("Cannot overwrite existing file!")
        hdf = openHDF5File(filename, 'a')
        entryName = "data"
        #entry
        nxEntry = hdf.require_group(entryName)
        if 'NX_class' not in nxEntry.attrs:
            nxEntry.attrs['NX_class'] = 'NXentry'.encode('utf-8')
        elif nxEntry.attrs['NX_class'] != 'NXentry'.encode('utf-8'):
            #should I raise an error?
            pass

        nxEntry['title'] = "PyMca saved 3D Array".encode('utf-8')
        nxEntry['start_time'] = getDate().encode('utf-8')
        nxData = nxEntry.require_group('NXdata')
        if ('NX_class' not in nxData.attrs):
            nxData.attrs['NX_class'] = 'NXdata'.encode('utf-8')
        elif nxData.attrs['NX_class'] != 'NXdata'.encode('utf-8'):
            #should I raise an error?
            pass
        if modify:
            if interpretation in ["image", "image".encode('utf-8')]:
                if compression:
                    if DEBUG:
                        print("Saving compressed and chunked dataset")
                    #risk of taking a 10 % more space in disk
                    chunk1 = int(shape[1] / 10)
                    if chunk1 == 0:
                        chunk1 = shape[1]
                    for i in [11, 10, 8, 7, 5, 4]:
                        if (shape[1] % i) == 0:
                            chunk1 = int(shape[1] / i)
                            break
                    chunk2 = int(shape[2] / 10)
                    for i in [11, 10, 8, 7, 5, 4]:
                        if (shape[2] % i) == 0:
                            chunk2 = int(shape[2] / i)
                            break
                    dset = nxData.require_dataset('data',
                                       shape=shape,
                                       dtype=dtype,
                                       chunks=(1, chunk1, chunk2),
                                       compression=compression)
                else:
                    if DEBUG:
                        print("Saving not compressed and not chunked dataset")
                    #print not compressed -> Not chunked
                    dset = nxData.require_dataset('data',
                                       shape=shape,
                                       dtype=dtype,
                                       compression=None)
                for i in range(data.shape[-1]):
                    tmp = data[:, :, i:i + 1]
                    tmp.shape = 1, shape[1], shape[2]
                    dset[i, 0:shape[1], :] = tmp
                    print("Saved item %d of %d" % (i + 1, data.shape[-1]))
            elif 0:
                #if I do not match the input and output shapes it takes ages
                #to save the images as spectra. However, it is much faster
                #when performing spectra operations.
                dset = nxData.require_dataset('data',
                               shape=shape,
                               dtype=dtype,
                               chunks=(1, shape[1], shape[2]))
                for i in range(data.shape[1]):  # shape[0]
                    chunk = numpy.zeros((1, data.shape[2], data.shape[0]),
                                        dtype)
                    for k in range(data.shape[0]):  # shape[2]
                        if 0:
                            tmpData = data[k:k + 1]
                            for j in range(data.shape[2]):  # shape[1]
                                tmpData.shape = data.shape[1], data.shape[2]
                                chunk[0, j, k] = tmpData[i, j]
                        else:
                            tmpData = data[k:k + 1, i, :]
                            tmpData.shape = -1
                            chunk[0, :, k] = tmpData
                    print("Saving item %d of %d" % (i, data.shape[1]))
                    dset[i, :, :] = chunk
            else:
                #if I do not match the input and output shapes it takes ages
                #to save the images as spectra. This is a very fast saving, but
                #the performance is awful when reading.
                if compression:
                    if DEBUG:
                        print("Saving compressed and chunked dataset")
                    dset = nxData.require_dataset('data',
                               shape=shape,
                               dtype=dtype,
                               chunks=(shape[0], shape[1], 1),
                               compression=compression)
                else:
                    if DEBUG:
                        print("Saving not compressed and not chunked dataset")
                    dset = nxData.require_dataset('data',
                               shape=shape,
                               dtype=dtype,
                               compression=None)
                for i in range(data.shape[0]):
                    tmp = data[i:i + 1, :, :]
                    tmp.shape = shape[0], shape[1], 1
                    dset[:, :, i:i + 1] = tmp
        else:
            if compression:
                if DEBUG:
                    print("Saving compressed and chunked dataset")
                chunk1 = int(shape[1] / 10)
                if chunk1 == 0:
                    chunk1 = shape[1]
                for i in [11, 10, 8, 7, 5, 4]:
                    if (shape[1] % i) == 0:
                        chunk1 = int(shape[1] / i)
                        break
                chunk2 = int(shape[2] / 10)
                if chunk2 == 0:
                    chunk2 = shape[2]
                for i in [11, 10, 8, 7, 5, 4]:
                    if (shape[2] % i) == 0:
                        chunk2 = int(shape[2] / i)
                        break
                if DEBUG:
                    print("Used chunk size = (1, %d, %d)" % (chunk1, chunk2))
                dset = nxData.require_dataset('data',
                               shape=shape,
                               dtype=dtype,
                               chunks=(1, chunk1, chunk2),
                               compression=compression)
            else:
                if DEBUG:
                    print("Saving not compressed and notchunked dataset")
                dset = nxData.require_dataset('data',
                               shape=shape,
                               dtype=dtype,
                               compression=None)
            tmpData = numpy.zeros((1, data.shape[1], data.shape[2]),
                                  data.dtype)
            for i in range(data.shape[0]):
                tmpData[0:1] = data[i:i + 1]
                dset[i:i + 1] = tmpData[0:1]
                print("Saved item %d of %d" % (i + 1, data.shape[0]))

        dset.attrs['signal'] = "1".encode('utf-8')
        if interpretation is not None:
            dset.attrs['interpretation'] = interpretation.encode('utf-8')
        axesAttribute = []
        for i in range(len(shape)):
            if axes is None:
                dim = numpy.arange(shape[i]).astype(numpy.float32)
                dimlabel = 'dim_%d' % i
            elif axes[i] is not None:
                dim = axes[i]
                try:
                    dimlabel = "%s" % labels[i]
                except:
                    dimlabel = 'dim_%d' % i
            else:
                dim = numpy.arange(shape[i]).astype(numpy.float32)
                dimlabel = 'dim_%d' % i
            axesAttribute.append(dimlabel)
            adset = nxData.require_dataset(dimlabel,
                                   dim.shape,
                                   dim.dtype,
                                   compression=None)
            adset[:] = dim[:]
            adset.attrs['axis'] = i + 1
        dset.attrs['axes'] = (":".join(axesAttribute)).encode('utf-8')
        nxEntry['end_time'] = getDate().encode('utf-8')
        if mode.lower() == 'nexus+':
            #create link
            g = h5py.h5g.open(hdf.fid, '/'.encode('utf-8'))
            g.link('/data/NXdata/data'.encode('utf-8'),
                   '/data/data'.encode('utf-8'),
                   h5py.h5g.LINK_HARD)

    elif mode.lower() == 'simplest':
        if os.path.exists(filename):
            try:
                os.remove(filename)
            except:
                raise IOError("Cannot overwrite existing file!")
        hdf = h5py.File(filename, 'a')
        if compression:
            hdf.require_dataset('data',
                           shape=shape,
                           dtype=dtype,
                           data=data,
                           chunks=(1, shape[1], shape[2]),
                           compression=compression)
        else:
            hdf.require_dataset('data',
                           shape=shape,
                           data=data,
                           dtype=dtype,
                           compression=None)
    else:
        if os.path.exists(filename):
            try:
                os.remove(filename)
            except:
                raise IOError("Cannot overwrite existing file!")
        shape = data.shape
        dtype = data.dtype
        hdf = h5py.File(filename, 'a')
        dataGroup = hdf.require_group('data')
        dataGroup.require_dataset('data',
                           shape=shape,
                           dtype=dtype,
                           data=data,
                           chunks=(1, shape[1], shape[2]))
    hdf.flush()
    hdf.close()

def main():
    a = numpy.arange(1000000.)
    a.shape = 20, 50, 1000
    save3DArrayAsHDF5(a, '/test.h5', mode='nexus+', interpretation='image')
    getHDF5FileInstanceAndBuffer('/test2.h5', (100, 100, 100))
    print("Date String = ", getDate())

if __name__ == "__main__":
    main()

