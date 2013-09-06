#/*##########################################################################
# Copyright (C) 2012 European Synchrotron Radiation Facility
#
# This file is part of the PyMca X-ray Fluorescence Toolkit developed at
# the ESRF by the Software group.
#
# This file is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# This file is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
#############################################################################*/
__author__ = "V.A. Sole - ESRF Data Analysis"
"""

A Stack plugin is a module that will be automatically added to the PyMca stack windows
in order to perform user defined operations on the data stack.

These plugins will be compatible with any stack window that provides the functions:
    #data related
    getStackDataObject
    getStackData
    getStackInfo
    setStack

    #images related
    addImage
    removeImage
    replaceImage

    #mask related
    setSelectionMask
    getSelectionMask

    #displayed curves
    getActiveCurve
    getGraphXLimits
    getGraphYLimits

    #information method
    stackUpdated
    selectionMaskUpdated
"""
import numpy
try:
    from PyMca import StackPluginBase
except ImportError:
    from . import StackPluginBase

DEBUG = 0

class ReverseStackPlugin(StackPluginBase.StackPluginBase):
    def __init__(self, stackWindow, **kw):
        StackPluginBase.DEBUG = DEBUG
        StackPluginBase.StackPluginBase.__init__(self, stackWindow, **kw)
        self.methodDict = {}
        text  = "Replace current stack by one\n"
        text += "with odd rows reversed."
        function = self.reverseOddRows
        info = text
        icon = None
        self.methodDict["Reverse Odd Rows"] =[function,
                                              info,
                                              icon]
        text  = "Replace current stack by one\n"
        text += "with even rows reversed."
        function = self.reverseEvenRows
        info = text
        icon = None
        self.methodDict["Reverse Even Rows"] =[function,
                                              info,
                                              icon]
        text  = "Replace current stack by one\n"
        text += "with odd columns reversed."
        function = self.reverseOddColumns
        info = text
        icon = None
        self.methodDict["Reverse Odd Columns"] =[function,
                                                 info,
                                                 icon]
        text  = "Replace current stack by one\n"
        text += "with odd columns reversed."
        function = self.reverseEvenColumns
        info = text
        icon = None
        self.methodDict["Reverse Even Columns"] =[function,
                                              info,
                                              icon]


        self.__methodKeys = ["Reverse Odd Rows",
                             "Reverse Even Rows",
                             "Reverse Odd Columns",
                             "Reverse Even Columns"]
        
    #Methods implemented by the plugin
    def getMethods(self):
        return self.__methodKeys

    def getMethodToolTip(self, name):
        return self.methodDict[name][1]

    def getMethodPixmap(self, name):
        return self.methodDict[name][2]

    def applyMethod(self, name):
        return self.methodDict[name][0]()

    def reverseOddRows(self):
        self.reverseRows(offset=1)

    def reverseEvenRows(self):
        self.reverseRows(offset=0)

    def reverseOddColumns(self):
        self.reverseColumns(offset=1)

    def reverseEvenColumns(self):
        self.reverseColumns(offset=0)

    def reverseRows(self, offset=1):
        stack = self.getStackDataObject()
        if not isinstance(stack.data, numpy.ndarray):
            text = "This method does not work with dynamically loaded stacks"
            raise TypeError(text)
        mcaIndex = stack.info.get('McaIndex', -1)
        if mcaIndex in [-1, 2]:
            ndata = stack.data.shape[1]
            limit = 0.5 * ndata
            for i in range(offset, stack.data.shape[0], 2):
                j = 0
                while j < limit:
                    tmp = stack.data[i, j, :] * 1
                    stack.data[i, j, :] = stack.data[i,(ndata-j-1),:] * 1
                    stack.data[i,(ndata-j-1),:] = tmp
                    j += 1
        elif mcaIndex == 0:
            ndata = stack.data.shape[2]
            limit = 0.5 * ndata
            for i in range(offset, stack.data.shape[1], 2):
                j = 0
                while j < limit:
                    tmp = stack.data[:, i, j] * 1
                    stack.data[:, i, j] = stack.data[:, i,(ndata-j-1)] * 1
                    stack.data[:, i,(ndata-j-1)] = tmp
                    j += 1            
        else:
            raise ValueError("Invalid 1D index %d" % mcaIndex)
        self.setStack(stack) 

    def reverseColumns(self, offset=1):
        stack = self.getStackDataObject()
        if not isinstance(stack.data, numpy.ndarray):
            text = "This method does not work with dynamically loaded stacks"
            raise TypeError(text)
        mcaIndex = stack.info.get('McaIndex', -1)
        if mcaIndex in [-1, 2]:
            ndata = stack.data.shape[0]
            limit = 0.5 * ndata
            for i in range(offset, stack.data.shape[1], 2):
                j = 0
                while j < limit:
                    tmp = stack.data[j, i, :] * 1
                    stack.data[j, i, :] = stack.data[(ndata-j-1), i,:] * 1
                    stack.data[(ndata-j-1), i,:] = tmp
                    j += 1
        elif mcaIndex == 0:
            ndata = stack.data.shape[1]
            limit = 0.5 * ndata
            for i in range(offset, stack.data.shape[2], 2):
                j = 0
                while j < limit:
                    tmp = stack.data[:, j, i] * 1
                    stack.data[:, j, i] = stack.data[:,(ndata-j-1), i] * 1
                    stack.data[:, (ndata-j-1), i] = tmp
                    j += 1            
        else:
            raise ValueError("Invalid 1D index %d" % mcaIndex)
        self.setStack(stack) 

MENU_TEXT = "Stack Row or Column Reversing"
def getStackPluginInstance(stackWindow, **kw):
    ob = ReverseStackPlugin(stackWindow)
    return ob
