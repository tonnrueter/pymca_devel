#/*##########################################################################
# Copyright (C) 2004-2012 European Synchrotron Radiation Facility
#
# This file is part of the PyMca X-ray Fluorescence Toolkit developed at
# the ESRF by the Software group.
#
# This toolkit is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# PyMca is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# PyMca; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# PyMca follows the dual licensing model of Riverbank's PyQt and cannot be
# used as a free plugin for a non-free program.
#
# Please contact the ESRF industrial unit (industry@esrf.fr) if this license
# is a problem for you.
#############################################################################*/
import sys
import os
import numpy
from PyMca import DataObject
from PyMca import specfilewrapper as specfile
from PyMca import SpecFileDataSource

HDF5 = False
try:
    import h5py
    HDF5 = True
except:
    pass
SOURCE_TYPE = "SpecFileStack"
DEBUG = 0

X_AXIS=0
Y_AXIS=1
Z_AXIS=2

class SpecFileStack(DataObject.DataObject):
    def __init__(self, filelist = None):
        DataObject.DataObject.__init__(self)
        self.incrProgressBar=0
        self.__keyList = []
        if filelist is not None:
            if type(filelist) != type([]):
                filelist = [filelist]
            if len(filelist) == 1:
                self.loadIndexedStack(filelist)
            else:
                self.loadFileList(filelist)

    def loadFileList(self, filelist, fileindex=0, shape=None):
        if type(filelist) == type(''):
            filelist = [filelist]
        self.__keyList = []
        self.sourceName = filelist
        self.__indexedStack = True
        self.sourceType = SOURCE_TYPE
        self.info = {}
        self.nbFiles=len(filelist)

        #read first file
        #get information
        tempInstance=SpecFileDataSource.SpecFileDataSource(filelist[0])
        keylist = tempInstance.getSourceInfo()['KeyList']
        nscans = len(keylist)        #that is the number of scans
        nmca = 0
        numberofdetectors = 0
        for key in keylist:
            info = tempInstance.getKeyInfo(key)
            numberofmca       = info['NbMca']
            if numberofmca > 0:
                numberofdetectors = info['NbMcaDet']
            scantype          = info["ScanType"]
            if numberofmca:
                nmca += numberofmca
        if numberofdetectors == 0:
            raise ValueError("No MCA found in file %s" % filelist[0])

        if (nscans > 1) and ((nmca/numberofdetectors) == nscans):
            SLOW_METHOD = True
        else:
            SLOW_METHOD = False
        #get last mca of first point
        key = "%s.1.%s" % (keylist[-1], numberofmca)
        dataObject = tempInstance._getMcaData(key)
        self.info.update(dataObject.info)
        arrRet = dataObject.data
        self.onBegin(self.nbFiles*nmca/numberofdetectors)

        self.incrProgressBar= 0
        if info['NbMcaDet'] > 1:
            #Should I generate a map for each mca and not just for the last one as I am doing?
            iterlist = range(info['NbMcaDet'],info['NbMca']+1,info['NbMcaDet'])
        else:
            iterlist = [1]
        if SLOW_METHOD and shape is None:
            self.data = numpy.zeros((self.nbFiles,
                                   nmca/numberofdetectors,
                                   arrRet.shape[0]),
                                   arrRet.dtype.char)
            filecounter         = 0
            for tempFileName in filelist:
                tempInstance=SpecFileDataSource.SpecFileDataSource(tempFileName)
                mca_number = -1
                for keyindex in keylist:
                    info = tempInstance.getKeyInfo(keyindex)
                    numberofmca       = info['NbMca']                    
                    if numberofmca <= 0:
                        continue
                    key = "%s.1.%s" % (keyindex, numberofmca)
                    dataObject = tempInstance._getMcaData(key)
                    arrRet = dataObject.data
                    mca_number += 1
                    for i in iterlist:
                        #mcadata = scan_obj.mca(i)
                        self.data[filecounter,
                                  mca_number,
                                  :] = arrRet[:]
                        self.incrProgressBar += 1
                        self.onProgress(self.incrProgressBar)
                filecounter += 1
        elif shape is None and (self.nbFiles == 1) and (iterlist == [1]):
            #it can only be here if there is one file
            #it can only be here if there is only one scan
            #it can only be here if there is only one detector
            self.data = numpy.zeros((1,
                               numberofmca,
                               arrRet.shape[0]),
                               arrRet.dtype.char)
            for tempFileName in filelist:
                tempInstance=specfile.Specfile(tempFileName)
                #it can only be here if there is one scan per file
                #prevent problems if the scan number is different
                #scan = tempInstance.select(keylist[-1])
                scan = tempInstance[-1]
                iterationList = range(scan.nbmca()) 
                for i in iterationList:
                    #mcadata = scan_obj.mca(i)
                    self.data[0,
                              i,
                              :] = scan.mca(i+1)[:]
                    self.incrProgressBar += 1
                    self.onProgress(self.incrProgressBar)
                filecounter = 1
        elif shape is None:
            #it can only be here if there is one scan per file
            try:
                self.data = numpy.zeros((self.nbFiles,
                                   numberofmca/numberofdetectors,
                                   arrRet.shape[0]),
                                   arrRet.dtype.char)
                filecounter         = 0
                for tempFileName in filelist:
                    tempInstance=specfile.Specfile(tempFileName)
                    #it can only be here if there is one scan per file
                    #prevent problems if the scan number is different
                    #scan = tempInstance.select(keylist[-1])
                    scan = tempInstance[-1]
                    for i in iterlist:
                        #mcadata = scan_obj.mca(i)
                        self.data[filecounter,
                                  0,
                                  :] = scan.mca(i)[:]
                        self.incrProgressBar += 1
                        self.onProgress(self.incrProgressBar)
                    filecounter += 1
            except MemoryError:
                qtflag = False
                if ('PyQt4.QtCore' in sys.modules) or\
                   ('PyMcaQt' in sys.modules) or\
                   ('PyMca.PyMcaQt' in sys.modules):
                    qtflag = True
                hdf5done = False
                if HDF5 and qtflag:
                    import PyMcaQt as qt
                    import ArraySave
                    msg=qt.QMessageBox.information( None,
                      "Memory error\n",
                      "Do you want to convert your data to HDF5?\n",
                      qt.QMessageBox.Yes,qt.QMessageBox.No)
                    if msg != qt.QMessageBox.No:
                        hdf5file = qt.QFileDialog.getSaveFileName(None,
                                    "Please select output file name",
                                    os.path.dirname(filelist[0]),
                                    "HDF5 files *.h5")
                        if not len(hdf5file):
                            raise IOError("Invalid output file")
                        hdf5file = qt.safe_str(hdf5file)
                        if not hdf5file.endswith(".h5"):
                            hdf5file += ".h5"

                        #get the final shape
                        from PyMca.RGBCorrelatorWidget import ImageShapeDialog
                        stackImageShape = self.nbFiles,\
                                     int(numberofmca/numberofdetectors)
                        dialog = ImageShapeDialog(None, shape =stackImageShape)
                        dialog.setModal(True)
                        ret = dialog.exec_()
                        if ret:
                            stackImageShape= dialog.getImageShape()
                            dialog.close()
                            del dialog
                        hdf, self.data =  ArraySave.getHDF5FileInstanceAndBuffer(hdf5file,
                                       (stackImageShape[0],
                                        stackImageShape[1],
                                        arrRet.shape[0]),
                                       compression=None,
                                       interpretation="spectrum")
                        nRow = 0
                        nCol = 0
                        for tempFileName in filelist:
                            tempInstance=specfile.Specfile(tempFileName)
                            #it can only be here if there is one scan per file
                            #prevent problems if the scan number is different
                            #scan = tempInstance.select(keylist[-1])
                            scan = tempInstance[-1]
                            nRow = int(self.incrProgressBar/stackImageShape[1])
                            nCol = self.incrProgressBar%stackImageShape[1]
                            for i in iterlist:
                                #mcadata = scan_obj.mca(i)
                                self.data[nRow,
                                          nCol,
                                          :] = scan.mca(i)[:]
                                self.incrProgressBar += 1
                                self.onProgress(self.incrProgressBar)
                        hdf5done = True
                        hdf.flush()
                    self.onEnd()
                    self.info["SourceType"] = "HDF5Stack1D"
                    self.info["McaIndex"] = 2
                    self.info["FileIndex"] = 0
                    self.info["SourceName"] = [hdf5file]
                    self.info["NumberOfFiles"] = 1
                    self.info["Size"]       = 1
                    return
                else:
                    raise
        else:
            sampling_order = 1
            s0 = shape[0]
            s1 = shape[1]
            MEMORY_ERROR = False
            try:
                self.data = numpy.zeros((shape[0],
                                   shape[1],
                                   arrRet.shape[0]),
                                   arrRet.dtype.char)
            except MemoryError:
                try:
                    self.data = numpy.zeros((shape[0],
                                   shape[1],
                                   arrRet.shape[0]),
                                   numpy.float32)
                except MemoryError:                    
                    MEMORY_ERROR = True
            while MEMORY_ERROR:
                try:
                    for i in range(5):
                        print("\7")
                    sampling_order += 1
                    print("**************************************************")
                    print(" Memory error!, attempting %dx%d sub-sampling " %\
                          (sampling_order, sampling_order))
                    print("**************************************************")
                    s0 = int(shape[0]/sampling_order)
                    s1 = int(shape[1]/sampling_order)
                    #if shape[0] % sampling_order:
                    #    s0 = s0 + 1
                    #if shape[1] % sampling_order:
                    #    s1 = s1 + 1
                    self.data = numpy.zeros((s0, s1,
                                             arrRet.shape[0]),
                                             numpy.float32)
                    MEMORY_ERROR = False
                except MemoryError:
                    pass
            filecounter         = 0
            for j in range(s0):
                filecounter = (j * sampling_order) * shape[1]
                for k in range(s1):
                    tempFileName = filelist[filecounter]
                    tempInstance=specfile.Specfile(tempFileName)
                    if tempInstance is None:
                        if not os.path.exists(tempFileName):
                            print("File %s does not exists"  % tempFileName)
                            raise IOError(\
                                "File %s does not exists"  % tempFileName)
                    scan = tempInstance.select(keylist[-1])
                    for i in iterlist:
                        #sum the present mcas
                        self.data[j,
                                  k,
                                  :] += scan.mca(i)[:]
                        self.incrProgressBar += 1
                        self.onProgress(self.incrProgressBar)
                    filecounter += sampling_order
            self.nbFiles = s0 * s1
        self.onEnd()

        """
        # Scan types
        # ----------    
        #SF_EMPTY       = 0        # empty scan
        #SF_SCAN        = 1        # non-empty scan
        #SF_MESH        = 2        # mesh scan
        #SF_MCA         = 4        # single mca
        #SF_NMCA        = 8        # multi mca (more than 1 mca per acq)

        case = None
        if scantype == (SpecFileDataSource.SF_MESH + \
                        SpecFileDataSource.SF_MCA):
            # SINGLE MESH + SINGLE MCA
            # nfiles  = 1
            # nscans  = 1
            # nmca    = 1
            # there is a danger if it can be considered an indexed file ...
            pass

        elif scantype == (SpecFileDataSource.SF_MESH + \
                        SpecFileDataSource.SF_NMCA):
            # SINGLE MESH + MULTIPLE MCA
            # nfiles  = 1
            # nscans  = 1
            # nmca    > 1
            # there is a danger if it can be considered an indexed file ...
            #for the time being I take last mca
            pass

        elif scantype == (SpecFileDataSource.SF_SCAN+ \
                          SpecFileDataSource.SF_MCA):
            #Assumed scans containing always 1 detector
            pass
        
        elif scantype == (SpecFileDataSource.SF_MCA):
            #Assumed scans containing always 1 detector
            pass

        elif scantype == (SpecFileDataSource.SF_SCAN+ \
                          SpecFileDataSource.SF_NMCA):
            #Assumed scans containing the same number of detectors
            #for the time being I take last mca
            pass
        
        elif scantype == (SpecFileDataSource.SF_NMCA):
            #Assumed scans containing the same number of detectors
            #for the time being I take last mca
            pass

        else:
            raise ValueError, "Unhandled scan type = %s" % scantype

        """

        self.__nFiles         = self.nbFiles
        self.__nImagesPerFile = 1
        shape = self.data.shape
        for i in range(len(shape)):
            key = 'Dim_%d' % (i+1,)
            self.info[key] = shape[i]
        self.info["SourceType"] = SOURCE_TYPE
        self.info["SourceName"] = self.sourceName
        self.info["Size"]       = self.__nFiles * self.__nImagesPerFile
        self.info["NumberOfFiles"] = self.__nFiles * 1
        self.info["FileIndex"] = fileindex

    def onBegin(self, n):
        pass

    def onProgress(self, n):
        pass

    def onEnd(self):
        pass

    def loadIndexedStack(self,filename,begin=None,end=None, skip = None, fileindex=0):
        #if begin is None: begin = 0
        if type(filename) == type([]):
            filename = filename[0]
        if not os.path.exists(filename):
            raise IOError("File %s does not exists" % filename)
        name = os.path.basename(filename)
        n = len(name)
        i = 1
        numbers = ['0', '1', '2', '3', '4', '5',
                   '6', '7', '8','9']
        while (i <= n):
            c = name[n-i:n-i+1]
            if c in ['0', '1', '2',
                                '3', '4', '5',
                                '6', '7', '8',
                                '9']:
                break
            i += 1
        suffix = name[n-i+1:]
        if len(name) == len(suffix):
            #just one file, one should use standard widget
            #and not this one.
            self.loadFileList(filename, fileindex=fileindex)
        else:
            nchain = []
            while (i<=n):
                c = name[n-i:n-i+1]
                if c not in ['0', '1', '2',
                                    '3', '4', '5',
                                    '6', '7', '8',
                                    '9']:
                    break
                else:
                    nchain.append(c)
                i += 1
            number = ""
            nchain.reverse()
            for c in nchain:
                number += c
            fformat = "%" + "0%dd" % len(number)
            if (len(number) + len(suffix)) == len(name):
                prefix = ""
            else:
                prefix = name[0:n-i+1]
            prefix = os.path.join(os.path.dirname(filename),prefix)
            if not os.path.exists(prefix + number + suffix):
                print("Internal error in EDFStack")
                print("file should exist: %s" % (prefix + number + suffix))
                return
            i = 0
            if begin is None:
                begin = 0
                testname = prefix+ fformat % begin+suffix
                while not os.path.exists(prefix + fformat % begin+suffix):
                    begin += 1
                    testname = prefix + fformat % begin+suffix
                    if len(testname) > len(filename):break
                i = begin
            else:
                i = begin
            if not os.path.exists(prefix+ fformat % i+suffix):
                raise ValueError("Invalid start index file = %s" % \
                      (prefix+ fformat % i+suffix))
            f = prefix + fformat % i+suffix
            filelist = []
            while os.path.exists(f):
                filelist.append(f)
                i += 1
                if end is not None:
                    if i > end:
                        break
                f = prefix + fformat % i+suffix
            self.loadFileList(filelist, fileindex=fileindex)

    def getSourceInfo(self):
        sourceInfo = {}
        sourceInfo["SourceType"]=SOURCE_TYPE
        if self.__keyList == []:
            for i in range(1, self.__nFiles + 1):
                for j in range(1, self.__nImages + 1):
                    self.__keyList.append("%d.%d" % (i,j))
        sourceInfo["KeyList"]= self.__keyList

    def getKeyInfo(self, key):
        print("Not implemented")
        return {}

    def isIndexedStack(self):
        return self.__indexedStack
    
    def getZSelectionArray(self,z=0):
        return (self.data[:,:,z]).astype(numpy.float)
        
    def getXYSelectionArray(self,coord=(0,0)):
        x,y=coord    
        return (self.data[y,x,:]).astype(numpy.float)

if __name__ == "__main__":
    import time
    import sys
    t0= time.time()
    stack = SpecFileStack()
    #stack.loadIndexedStack("Z:\COTTE\ch09\ch09__mca_0005_0000_0070.edf")
    if len(sys.argv) > 1:
        stack.loadIndexedStack(sys.argv[1])
    else:
        stack.loadIndexedStack("..\..\mca3\c449b01_001.mca")
    shape = stack.data.shape
    print("elapsed = %f" % (time.time() - t0))
    #guess the MCA
    imax = 0
    for i in range(len(shape)):
        if shape[i] > shape[imax]:
            imax = i

    print("selections ")
    print("getZSelectionArray  shape = ", stack.getZSelectionArray().shape)
    print("getXYSelectionArray shape = ", stack.getXYSelectionArray().shape)

    from PyMca import PyMcaQt as qt
    app = qt.QApplication([])
    qt.QObject.connect(app, qt.SIGNAL("lastWindowClosed()"),
                       app, qt.SLOT("quit()"))
    if 1:
        from PyMca import RGBCorrelatorGraph
        w = RGBCorrelatorGraph.RGBCorrelatorGraph()
        graph = w.graph
    else:
        from PyMca import QtBlissGraph
        w = QtBlissGraph.QtBlissGraph()
        graph = w
    print("shape sum 0 = ", numpy.sum(stack.data, 0).shape)
    print("shape sum 1 = ", numpy.sum(stack.data, 1).shape)
    print("shape sum 2 = ", numpy.sum(stack.data, 2).shape)
    a = numpy.sum(stack.data, imax)
    print(a.shape)
    graph.setX1AxisLimits(0, a.shape[0])
    if 0:
        w.setY1AxisLimits(0, a.shape[1])
        w.setY1AxisInverted(True)
    else:
        graph.setY1AxisInverted(True)
        graph.setY1AxisLimits(0, a.shape[1])
    graph.imagePlot(a, ymirror=0)
    if imax == 0:
        graph.x1Label('Column Number')
    else:
        graph.x1Label('Row Number')
    graph.ylabel('File Number')
    w.show()

    if imax == 0:
        mcaData0 = numpy.sum(numpy.sum(stack.data, 2),1)
    else:
        mcaData0 = numpy.sum(numpy.sum(stack.data, 2),0)

    from PyMca import McaWindow
    mca = McaWindow.McaWidget()
    sel = {}
    sel['SourceName'] = "Specfile Stack"
    sel['Key']        = "SUM"
    sel['legend']     = "EDF Stack SUM"
    mcaData = DataObject.DataObject()
    mcaData.info = {'McaCalib': [0 , 2.0 ,0],
                    "selectiontype":"1D",
                    "SourceName":"Specfile Stack",
                    "Key":"SUM"}
    mcaData.x = [numpy.arange(len(mcaData0)).astype(numpy.float)]
    mcaData.y = [mcaData0]
    sel['dataobject'] = mcaData
    mca.show()
    mca._addSelection([sel])
    graph.enableSelection(True)
    def graphSlot(ddict):
        if ddict['event'] == "MouseSelection":
            ix1 = int(ddict['xmin'])
            ix2 = int(ddict['xmax'])+1
            iy1 = int(ddict['xmin'])
            iy2 = int(ddict['xmax'])+1
            if imax == 0:
                selectedData = numpy.sum(numpy.sum(stack.data[:,ix1:ix2, iy1:iy2], 2),1)
            else:
                selectedData = numpy.sum(numpy.sum(stack.data[ix1:ix2,:, iy1:iy2], 2),0)
            sel = {}
            sel['SourceName'] = "Specfile Stack"
            sel['Key'] = "Selection"
            sel["selectiontype"] = "1D"
            sel['legend'] = "EDF Stack Selection"
            selDataObject = DataObject.DataObject()
            selDataObject.info={'McaCalib': [100 , 2.0 ,0],
                                "selectiontype":"1D",
                                "SourceName":"EDF Stack Selection",
                                "Key":"Selection"}
            selDataObject.x = [numpy.arange(len(mcaData0)).astype(numpy.float)]
            selDataObject.y = [selectedData]
            sel['dataobject'] = selDataObject
            mca._addSelection([sel])
    qt.QObject.connect(graph, qt.SIGNAL('QtBlissGraphSignal'),
                       graphSlot)
    #w.replot()
    if qt.qVersion() < '4.0.0':
        app.exec_loop()
    else:
        app.exec_()
