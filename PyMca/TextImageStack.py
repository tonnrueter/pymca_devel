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
import numpy
import sys
import os
from PyMca import DataObject
        
SOURCE_TYPE = "EdfFileStack"
DEBUG = 0

class TextImageStack(DataObject.DataObject):
    def __init__(self, filelist = None, imagestack=None, dtype=None):
        DataObject.DataObject.__init__(self)
        self.incrProgressBar=0
        self.__keyList = []
        if imagestack is None:
            self.__imageStack = True
        else:
            self.__imageStack = imagestack
        self.__dtype = dtype
        if filelist is not None:
            if type(filelist) != type([]):
                filelist = [filelist]
            if len(filelist) == 1:
                self.loadIndexedStack(filelist)
            else:
                self.loadFileList(filelist)

    def loadFileList(self, filelist, fileindex=0):
        if type(filelist) == type(''):
            filelist = [filelist]
        self.__keyList = []
        self.sourceName = filelist
        self.__indexedStack = True
        self.sourceType = SOURCE_TYPE
        self.info = {}
        self.__nFiles=len(filelist)
        #read first file
        arrRet = numpy.loadtxt(filelist[0])
        if self.__dtype is None:
            self.__dtype = arrRet.dtype
        self.__nImagesPerFile = 1
            
        #try to allocate the memory
        shape = self.__nFiles, arrRet.shape[0], arrRet.shape[1]
        samplingStep = 1
        try:
            self.data = numpy.zeros(shape, self.__dtype)
        except (MemoryError, ValueError):
            for i in range(3):
                print("\7")
            samplingStep = None
            i = 2
            while samplingStep is None:
                print("**************************************************")
                print(" Memory error!, attempting %dx%d sampling reduction ") % (i,i)
                print("**************************************************")
                s1, s2 = arrRet[::i, ::i].shape
                try:
                    self.data = numpy.zeros((self.__nFiles, s1, s2),
                                         self.__dtype)
                    samplingStep = i
                except (MemoryError, ValueError):
                    i += 1

        #fill the array
        self.onBegin(self.__nFiles)
        self.__imageStack = True
        self.incrProgressBar=0
        if samplingStep == 1:
            for tempFileName in filelist:
                self.data[self.incrProgressBar]=numpy.loadtxt(tempFileName,
                                                        dtype=self.__dtype)
                self.incrProgressBar += 1
                self.onProgress(self.incrProgressBar)
        else:
            for tempFileName in filelist:
                pieceOfStack=numpy.loadtxt(tempFileName, dtype=self.__dtype)
                self.data[self.incrProgressBar] = pieceOfStack[::samplingStep,
                                                               ::samplingStep]
                self.incrProgressBar += 1
                self.onProgress(self.incrProgressBar)
        self.onEnd()
        if self.__imageStack:
            self.info["McaIndex"] = 0
            self.info["FileIndex"] = 1
        else:
            self.info["McaIndex"] = 2
            self.info["FileIndex"] = 0
        self.info["SourceType"] = SOURCE_TYPE
        self.info["SourceName"] = self.sourceName
        self.info["NumberOfFiles"] = self.__nFiles * 1
        self.info["Size"] = self.__nFiles * self.__nImagesPerFile

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
                print("file should exist: %s " % (prefix + number + suffix))
                return
            i = 0
            if begin is None:
                begin = 0
                testname = prefix+fformat % begin+suffix
                while not os.path.exists(prefix+fformat % begin+suffix):
                    begin += 1
                    testname = prefix+fformat % begin+suffix
                    if len(testname) > len(filename):break
                i = begin
            else:
                i = begin
            if not os.path.exists(prefix+fformat % i+suffix):
                raise ValueError("Invalid start index file = %s" % \
                      (prefix+fformat % i+suffix))
            f = prefix+fformat % i+suffix
            filelist = []
            while os.path.exists(f):
                filelist.append(f)
                i += 1
                if end is not None:
                    if i > end:
                        break
                f = prefix+fformat % i+suffix
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

