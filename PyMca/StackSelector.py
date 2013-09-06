#!/usr/bin/env python
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
import copy
from PyMca import PyMcaQt as qt
from PyMca import PyMcaDirs
from PyMca import DataObject
from PyMca import OmnicMap
from PyMca import OpusDPTMap
from PyMca import LuciaMap
from PyMca import SupaVisioMap
from PyMca import AifiraMap
from PyMca import TextImageStack
from PyMca import TiffStack
from PyMca.QEDFStackWidget import QStack, QSpecFileStack
try:
    from PyMca import QHDF5Stack1D
    import h5py
    HDF5 = True
except ImportError:
    HDF5 = False
    pass

QTVERSION = qt.qVersion()
DEBUG = 0


class StackSelector(object):
    def __init__(self, parent=None):
        self.parent = parent

    def getStack(self, filelist=None, imagestack=None):
        if filelist in [None, []]:
            filelist, filefilter = self._getStackOfFiles(getfilter=True)
        else:
            filefilter = ""

        if not len(filelist):
            return None

        if filefilter == "":
            if HDF5:
                if h5py.is_hdf5(filelist[0]):
                    filefilter = "HDF5"

        fileindex = 0
        begin = None
        end = None
        aifirafile = False
        if len(filelist):
            PyMcaDirs.inputDir = os.path.dirname(filelist[0])
            #if we are dealing with HDF5, no more tests needed
            if not filefilter.upper().startswith('HDF5'):
                f = open(filelist[0], 'rb')
                #read 10 characters
                if sys.version < '3.0':
                    line = f.read(10)
                else:
                    line = str(f.read(10).decode())
                f.close()
            omnicfile = False
            if filefilter.upper().startswith('HDF5'):
                stack = QHDF5Stack1D.QHDF5Stack1D(filelist)
                omnicfile = True
            elif filefilter.upper().startswith('OPUS-DPT'):
                stack = OpusDPTMap.OpusDPTMap(filelist[0])
                omnicfile = True
            elif filefilter.upper().startswith("AIFIRA"):
                stack = AifiraMap.AifiraMap(filelist[0])
                omnicfile = True
                aifirafile = True
            elif filefilter.upper().startswith("SUPAVISIO"):
                stack = SupaVisioMap.SupaVisioMap(filelist[0])
                omnicfile = True
            elif filefilter.upper().startswith("TEXTIMAGE"):
                imagestack = True
                fileindex = 0
                stack = TextImageStack.TextImageStack(imagestack=True)
            elif filefilter.upper().startswith("IMAGE") and\
                 (filelist[0].upper().endswith("TIF") or\
                  filelist[0].upper().endswith("TIFF")):
                stack = TiffStack.TiffStack(imagestack=True)
            elif filefilter == "" and\
                 (filelist[0].upper().endswith("TIF") or\
                  filelist[0].upper().endswith("TIFF")):
                stack = TiffStack.TiffStack(imagestack=True)
            elif filefilter.upper().startswith("IMAGE"):
                if imagestack is None:
                    imagestack = True
                fileindex = 0
                stack = QStack(imagestack=imagestack)
            elif line[0] == "{":
                if filelist[0].upper().endswith("RAW"):
                    if imagestack is None:
                        imagestack=True
                stack = QStack(imagestack=imagestack)
            elif line[0:2] in ["II", "MM"]:
                if imagestack is None:
                    imagestack = True
                stack = QStack(imagestack=imagestack)
            elif line.startswith('Spectral'):
                stack = OmnicMap.OmnicMap(filelist[0])
                omnicfile = True
            elif line.startswith('#\tDate'):
                stack = LuciaMap.LuciaMap(filelist[0])
                omnicfile = True
            elif filelist[0].upper().endswith("RAW.GZ")or\
                 filelist[0].upper().endswith("EDF.GZ")or\
                 filelist[0].upper().endswith("CCD.GZ")or\
                 filelist[0].upper().endswith("RAW.BZ2")or\
                 filelist[0].upper().endswith("EDF.BZ2")or\
                 filelist[0].upper().endswith("CCD.BZ2"):
                if imagestack is None:
                    imagestack = True
                stack = QStack(imagestack=imagestack)
            elif filelist[0][-4:].upper() in ["PIGE", "PIGE"]:
                stack = SupaVisioMap.SupaVisioMap(filelist[0])
                omnicfile = True
            elif filelist[0][-3:].upper() in ["RBS"]:
                stack = SupaVisioMap.SupaVisioMap(filelist[0])
                omnicfile = True
            elif filelist[0][-3:].upper() in ["SPE"] and\
                 (line[0] not in ['$', '#']):
                #Roper Scientific format
                #handle it as MarCCD stack
                stack = QStack(imagestack=True)
            else:
                stack = QSpecFileStack()

        if len(filelist) == 1:
            if not omnicfile:
                try:
                    stack.loadIndexedStack(filelist[0], begin, end,
                                           fileindex=fileindex)
                except:
                    msg = qt.QMessageBox()
                    msg.setIcon(qt.QMessageBox.Critical)
                    msg.setText("%s" % sys.exc_info()[1])
                    if QTVERSION < '4.0.0':
                        msg.exec_loop()
                    else:
                        msg.exec_()
                    if DEBUG:
                        raise
        elif len(filelist):
            if not omnicfile:
                try:
                    stack.loadFileList(filelist, fileindex=fileindex)
                except:
                    msg = qt.QMessageBox()
                    msg.setIcon(qt.QMessageBox.Critical)
                    msg.setText("%s" % sys.exc_info()[1])
                    if QTVERSION < '4.0.0':
                        msg.exec_loop()
                    else:
                        msg.exec_()
                    if DEBUG:
                        raise
        if aifirafile:
            masterStack = DataObject.DataObject()
            masterStack.info = copy.deepcopy(stack.info)
            masterStack.data = stack.data[:, :, 0:1024]
            masterStack.info['Dim_2'] = int(masterStack.info['Dim_2'] / 2)

            slaveStack = DataObject.DataObject()
            slaveStack.info = copy.deepcopy(stack.info)
            slaveStack.data = stack.data[:, :, 1024:]
            slaveStack.info['Dim_2'] = int(slaveStack.info['Dim_2'] / 2)
            return [masterStack, slaveStack]
        else:
            return stack

    def getStackFromPattern(self, filepattern, begin, end, increment=None,
                            imagestack=None, fileindex=0):
        #get the first filename
        filename = filepattern % tuple(begin)
        if not os.path.exists(filename):
            raise IOError("Filename %s does not exist." % filename)
        #get the file list
        args = self.getFileListFromPattern(filepattern, begin, end,
                                           increment=increment)

        #get the file type
        f = open(args[0], 'rb')
        #read 10 characters
        line = f.read(10)
        f.close()

        specfile = False
        marCCD = False
        if line.startswith("II") or line.startswith("MM"):
            marCCD = True
        if line[0] == "\n":
            line = line[1:]
        if (line[0] == "{") or marCCD:
            if imagestack is None:
                if marCCD:
                    imagestack = True
            if imagestack:
                #prevent any modification
                fileindex = 0
            if filepattern is not None:
                #this dows not seem to put any trouble
                #(because of no redimensioning attempt)
                if False and (len(begin) != 1):
                    raise IOError("EDF stack redimensioning not supported yet")
            stack = QStack(imagestack=imagestack)
        elif line.startswith('Spectral'):
            stack = OmnicMap.OmnicMap(args[0])
        elif line.startswith('#\tDate:'):
            stack = LuciaMap.LuciaMap(args[0])
        elif args[0][-4:].upper() in ["PIGE", "PIXE"]:
            stack = SupaVisioMap.SupaVisioMap(args[0])
        elif args[0][-3:].upper() in ["RBS"]:
            stack = SupaVisioMap.SupaVisioMap(args[0])
        elif args[0][-3:].lower() in [".h5", "nxs", "hdf"]:
            if not HDF5:
                raise IOError(\
                    "No HDF5 support while trying to read an HDF5 file")
            stack = QHDF5Stack1D.QHDF5Stack1D(args)
        elif args[0].upper().endswith("RAW.GZ")or\
             args[0].upper().endswith("EDF.GZ")or\
             args[0].upper().endswith("CCD.GZ")or\
             args[0].upper().endswith("RAW.BZ2")or\
             args[0].upper().endswith("EDF.BZ2")or\
             args[0].upper().endswith("CCD.BZ2"):
            if imagestack is None:
                imagestack = True
            stack = QStack(imagestack=imagestack)            
        else:
            if HDF5:
                if h5py.is_hdf5(args[0]):
                    stack = QHDF5Stack1D.QHDF5Stack1D(args)
                else:
                    stack = QSpecFileStack()
                    specfile = True
            else:
                stack = QSpecFileStack()
                specfile = True

        if specfile and (len(begin) == 2):
            if increment is None:
                increment = [1] * len(begin)
            shape = (len(range(begin[0], end[0] + 1, increment[0])),
                     len(range(begin[1], end[1] + 1, increment[1])))
            stack.loadFileList(args, fileindex=fileindex, shape=shape)
        else:
            stack.loadFileList(args, fileindex=fileindex)
        return stack

    def _getFileList(self, fileTypeList, message=None, getfilter=None):
        if message is None:
            message = "Please select a file"
        if getfilter is None:
            getfilter = False
        wdir = PyMcaDirs.inputDir
        filterused = None
        if QTVERSION < '4.0.0':
            if sys.platform != 'darwin':
                filetypes = ""
                for filetype in fileTypeList:
                    filetypes += filetype + "\n"
                filelist = qt.QFileDialog.getOpenFileNames(filetypes,
                            wdir,
                            self.parent,
                            message,
                            message)
                if not len(filelist):
                    if getfilter:
                        return [], filterused
                    else:
                        return []
        else:
            #if (QTVERSION < '4.3.0') and (sys.platform != 'darwin'):
            if (PyMcaDirs.nativeFileDialogs) and (not getfilter):
                filetypes = ""
                for filetype in fileTypeList:
                    filetypes += filetype + "\n"
                filelist = qt.QFileDialog.getOpenFileNames(self.parent,
                        message,
                        wdir,
                        filetypes)
                if not len(filelist):
                    if getfilter:
                        return [], filterused
                    else:
                        return []
                else:
                    sample = qt.safe_str(filelist[0])
                    for filetype in fileTypeList:
                        ftype = filetype.replace("(", "")
                        ftype = ftype.replace(")", "")
                        extensions = ftype.split()[2:]
                        for extension in extensions:
                            if sample.endswith(extension[-3:]):
                                filterused = filetype
                                break
            else:
                fdialog = qt.QFileDialog(self.parent)
                fdialog.setModal(True)
                fdialog.setWindowTitle(message)
                if hasattr(qt, "QStringList"):
                    strlist = qt.QStringList()
                else:
                    strlist = []
                for filetype in fileTypeList:
                    strlist.append(filetype)
                fdialog.setFilters(strlist)
                fdialog.setFileMode(fdialog.ExistingFiles)
                fdialog.setDirectory(wdir)
                if QTVERSION > '4.3.0':
                    history = fdialog.history()
                    if len(history) > 6:
                        fdialog.setHistory(history[-6:])
                ret = fdialog.exec_()
                if ret == qt.QDialog.Accepted:
                    filelist = fdialog.selectedFiles()
                    if getfilter:
                        filterused = qt.safe_str(fdialog.selectedFilter())
                    fdialog.close()
                    del fdialog
                else:
                    fdialog.close()
                    del fdialog
                    if getfilter:
                        return [], filterused
                    else:
                        return []
        if sys.version < '3.0':
            try:
                filelist = [qt.safe_str(x) for x in filelist]
            except UnicodeEncodeError:
                filelist = [unicode(x) for x in filelist]
        if not(len(filelist)):
            return []
        PyMcaDirs.inputDir = os.path.dirname(filelist[0])
        if PyMcaDirs.outputDir is None:
            PyMcaDirs.outputDir = os.path.dirname(filelist[0])

        #This should not be there in order to allow the user other choices
        #filelist.sort()
        if getfilter:
            return filelist, filterused
        else:
            return filelist

    def _getStackOfFiles(self, getfilter=None):
        if getfilter is None:
            getfilter = False
        fileTypeList = ["EDF Files (*edf)",
                        "EDF Files (*ccd)",
                        "Specfile Files (*mca)",
                        "Specfile Files (*dat)",
                        "OMNIC Files (*map)",
                        "OPUS-DPT Files (*.DPT *.dpt)",
                        "HDF5 Files (*.nxs *.hdf *.h5)",
                        "AIFIRA Files (*DAT)",
                        "SupaVisio Files (*pige *pixe *rbs)",
                        "Image Files (*edf *ccd *raw *edf.gz *ccd.gz *raw.gz)",
                        "Image Files (*tif *tiff *TIF *TIFF)",
                        "TextImage Files (*txt)",
                        "All Files (*)"]
        if not HDF5:
            idx = fileTypeList.index("HDF5 Files (*.nxs *.hdf *.h5)")
            del fileTypeList[idx]
        message = "Open ONE indexed stack or SEVERAL files"
        return self._getFileList(fileTypeList, message=message,
                                 getfilter=getfilter)

    def getFileListFromPattern(self, pattern, begin, end, increment=None):
        if type(begin) == type(1):
            begin = [begin]
        if type(end) == type(1):
            end = [end]
        if len(begin) != len(end):
            raise ValueError(\
                "Begin list and end list do not have same length")
        if increment is None:
            increment = [1] * len(begin)
        elif type(increment) == type(1):
            increment = [increment]
        if len(increment) != len(begin):
            raise ValueError(\
                "Increment list and begin list do not have same length")
        fileList = []
        if len(begin) == 1:
            for j in range(begin[0], end[0] + increment[0], increment[0]):
                fileList.append(pattern % (j))
        elif len(begin) == 2:
            for j in range(begin[0], end[0] + increment[0], increment[0]):
                for k in range(begin[1], end[1] + increment[1], increment[1]):
                    fileList.append(pattern % (j, k))
        elif len(begin) == 3:
            raise ValueError("Cannot handle three indices yet.")
            for j in range(begin[0], end[0] + increment[0], increment[0]):
                for k in range(begin[1], end[1] + increment[1], increment[1]):
                    for l in range(begin[2], end[2] + increment[2], increment[2]):
                        fileList.append(pattern % (j, k, l))
        else:
            raise ValueError("Cannot handle more than three indices.")
        return fileList


if __name__ == "__main__":
    from PyMca import QStackWidget
    import getopt
    options = ''
    longoptions = ["fileindex=",
                   "filepattern=", "begin=", "end=", "increment=",
                   "nativefiledialogs=", "imagestack="]
    try:
        opts, args = getopt.getopt(
                     sys.argv[1:],
                     options,
                     longoptions)
    except:
        print(sys.exc_info()[1])
        sys.exit(1)
    fileindex = 0
    filepattern = None
    begin = None
    end = None
    imagestack = False
    increment = None
    for opt, arg in opts:
        if opt in '--begin':
            if "," in arg:
                begin = [int(x) for x in arg.split(",")]
            else:
                begin = [int(arg)]
        elif opt in '--end':
            if "," in arg:
                end = [int(x) for x in arg.split(",")]
            else:
                end = int(arg)
        elif opt in '--increment':
            if "," in arg:
                increment = [int(x) for x in arg.split(",")]
            else:
                increment = int(arg)
        elif opt in '--filepattern':
            filepattern = arg.replace('"', '')
            filepattern = filepattern.replace("'", "")
        elif opt in '--fileindex':
            fileindex = int(arg)
        elif opt in '--imagestack':
            imagestack = int(arg)
        elif opt in '--nativefiledialogs':
            if int(arg):
                PyMcaDirs.nativeFileDialogs = True
            else:
                PyMcaDirs.nativeFileDialogs = False
    if filepattern is not None:
        if (begin is None) or (end is None):
            raise ValueError(\
                "A file pattern needs at least a set of begin and end indices")
    app = qt.QApplication([])
    widget = QStackWidget.QStackWidget()
    w = StackSelector(widget)
    if filepattern is not None:
        #ignore the args even if present
        stack = w.getStackFromPattern(filepattern, begin, end,
                                      increment=increment,
                                      imagestack=imagestack)
    else:
        stack = w.getStack(args, imagestack=imagestack)
    if type(stack) == type([]):
        #aifira like, two stacks
        widget.setStack(stack[0])
        slave = QStackWidget.QStackWidget(master=False,
                                          rgbwidget=widget.rgbWidget)
        slave.setStack(stack[1])
        widget.setSlave(slave)
        stack = None
    else:
        widget.setStack(stack)
    widget.show()
    app.exec_()
