#/*##########################################################################
# Copyright (C) 2004-2013 European Synchrotron Radiation Facility
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
__author__ = "V.A. Sole - ESRF Data Analysis"
import numpy
from numpy import cos, sin
import sys
try:
    from PyMca import Plugin1DBase
except ImportError:
    from . import Plugin1DBase

import os
try:
    from PyMca import PyMcaQt as qt
except ImportError:
    print("WARNING: RegularMeshPlugin Using huge PyQt4 import")
    import PyQt4.Qt as qt

from PyMca import MaskImageWidget

DEBUG = 0

class RegularMeshPlugins(Plugin1DBase.Plugin1DBase):
    def __init__(self, plotWindow, **kw):
        Plugin1DBase.Plugin1DBase.__init__(self, plotWindow, **kw)
        self.methodDict = {}
        self.methodDict['Show Image'] = [self._convert,
                                             "Show mesh as image",
                                             None]
                           
        self.imageWidget = None
        
    #Methods to be implemented by the plugin
    def getMethods(self, plottype=None):
        """
        A list with the NAMES  associated to the callable methods
        that are applicable to the specified plot.

        Plot type can be "SCAN", "MCA", None, ...        
        """
        names = list(self.methodDict.keys())
        names.sort()
        return names

    def getMethodToolTip(self, name):
        """
        Returns the help associated to the particular method name or None.
        """
        return self.methodDict[name][1]

    def getMethodPixmap(self, name):
        """
        Returns the pixmap associated to the particular method name or None.
        """
        return self.methodDict[name][2]

    def applyMethod(self, name):
        """
        The plugin is asked to apply the method associated to name.
        """
        if DEBUG:
                self.methodDict[name][0]()
        else:
            try:
                self.methodDict[name][0]()
            except:
                print(sys.exc_info())
                raise

    def _convert(self):
        x, y, legend, info = self.getActiveCurve()
        self._x = x
        self._y = y
        if 'Header' not in info:
            raise ValueError("Active curve does not seem to be a mesh scan")
            
        header = info['Header'][0]

        item = header.split()
        if item[2] not in ['mesh', 'hklmesh']:
            raise ValueError("Active curve does not seem to be a mesh scan")

        self._xLabel = self.getGraphXTitle()
        self._yLabel = self.getGraphYTitle()

        self._motor0Mne = item[3]
        self._motor1Mne = item[7]

        #print("Scanned motors are %s and %s" % (motor0Mne, motor1Mne))
        
        #Assume an EXACTLY regular mesh for both motors
        self._motor0 = numpy.linspace(float(item[4]), float(item[5]), int(item[6])+1)
        self._motor1 = numpy.linspace(float(item[8]), float(item[9]), int(item[10])+1)

        try:
            if xLabel.upper() == motor0Mne.upper():
                self._motor0 = self._x
                self._motor0Mne = self._xLabel
            elif xLabel.upper() == motor1Mne.upper():
                self._motor1 = self._x
                self._motor1Mne = self._xLabel
            elif xLabel == info['selection']['cntlist'][0]:
                self._motor0 = self._x
                self._motor0Mne = self._xLabel
            elif xLabel == info['selection']['cntlist'][1]:
                self._motor1 = self._x
                self._motor1Mne = self._xLabel
        except:
            if DEBUG:
                print("XLabel should be one of the scanned motors")

        self._legend = legend
        self._info = info
        y.shape = len(self._motor1), len(self._motor0)
        if self.imageWidget is None:
            self.imageWidget = MaskImageWidget.MaskImageWidget(\
                                        imageicons=False,
                                        selection=False,
                                        profileselection=True,
                                        scanwindow=self)
        self.imageWidget.setImageData(y,
                                      xScale=(self._motor0[0],
                                              self._motor0[-1]),                                              
                                      yScale=(self._motor1[0],
                                              self._motor1[-1]))
        self.imageWidget.setXLabel(self._motor0Mne)
        self.imageWidget.setYLabel(self._motor1Mne)
        self.imageWidget.show()

MENU_TEXT = "RegularMeshPlugins"
def getPlugin1DInstance(plotWindow, **kw):
    ob = RegularMeshPlugins(plotWindow)
    return ob

if __name__ == "__main__":
    from PyMca import Plot1D
    app = qt.QApplication([])
    DEBUG = 1
    x = numpy.arange(100.)
    y = x * x
    plot = Plot1D.Plot1D()
    plot.addCurve(x, y, "dummy")
    plot.addCurve(x+100, -x*x)
    plugin = getPlugin1DInstance(plot)
    for method in plugin.getMethods():
        print(method, ":", plugin.getMethodToolTip(method))
    plugin.applyMethod(plugin.getMethods()[1])
    curves = plugin.getAllCurves()
    for curve in curves:
        print(curve[2])
    print("LIMITS = ", plugin.getGraphYLimits())
