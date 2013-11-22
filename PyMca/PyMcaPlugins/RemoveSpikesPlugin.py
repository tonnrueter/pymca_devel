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
__author__ = "T. Rueter - ESRF Data Analysis Unit"

try:
    from PyMca import Plugin1DBase
except ImportError:
    from . import Plugin1DBase

try:
    from PyMca.PyMcaSciPy.signal import medfilt1d
except ImportError:
    print("RemoveSpikesPlugin importing directly")
    from PyMca.PyMcaSciPy.signal import medfilt1d

from PyMca import PyMcaQt as qt
import numpy
    
DEBUG = 0
class RemoveSpikes(Plugin1DBase.Plugin1DBase):
    def __init__(self,  plotWindow,  **kw):
        Plugin1DBase.Plugin1DBase.__init__(self,  plotWindow,  **kw)
        self.methodDict = {
            'Apply to active curve':
                [self.removeSpikesActive,
                 'Apply sliding median filter to active curve',
                 None],
            'Apply to all curves':
                [self.removeSpikesAll,
                 'Apply sliding median filter to all curves',
                 None],
            'Configure median filter':
                [self.configureFilter,
                 'Set threshold and width of the filter',
                 None]
        }
        self.threshold = 0.66
        self.width = 9
    
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
        self.methodDict[name][0]()
        return

    def configureFilter(self):
        msg = qt.QDialog()
        msgLayout = qt.QGridLayout()
        buttonLayout = qt.QHBoxLayout()
        
        inpThreshold = qt.QDoubleSpinBox()
        inpThreshold.setRange(0.,10.)
        inpThreshold.setSingleStep(.1)
        inpThreshold.setValue(2.0)
        inpThreshold.setToolTip('Increase width for broad spikes')

        inpWidth = qt.QSpinBox()
        inpWidth.setRange(1,101)
        inpWidth.setSingleStep(2)
        inpWidth.setValue(self.width)
        inpWidth.setToolTip('Set low threshold for multiple spikes of different markedness')
        
        labelWidth = qt.QLabel('Width (must be odd)')
        labelThreshold = qt.QLabel('Threshold (multiple of deviation)')
        buttonOK = qt.QPushButton('Ok')
        buttonOK.clicked.connect(msg.accept)
        buttonCancel = qt.QPushButton('Cancel')
        buttonCancel.clicked.connect(msg.reject)
        
        allActiveBG = qt.QButtonGroup()
        buttonAll = qt.QCheckBox('Apply to All')
        buttonActive = qt.QCheckBox('Apply to Active')
        allActiveBG.addButton(buttonAll, 0)
        allActiveBG.addButton(buttonActive, 1)
        
        buttonLayout.addWidget(qt.HorizontalSpacer())
        buttonLayout.addWidget(buttonOK)
        buttonLayout.addWidget(buttonCancel)
        
        msgLayout.addWidget(labelWidth,0,0)
        msgLayout.addWidget(inpWidth,0,1)
        msgLayout.addWidget(labelThreshold,1,0)
        msgLayout.addWidget(inpThreshold,1,1)
        msgLayout.addWidget(buttonActive,2,0)
        msgLayout.addWidget(buttonAll,2,1)
        msgLayout.addLayout(buttonLayout,3,0,1,2)
        msg.setLayout(msgLayout)
        if msg.exec_():
            try:
                self.threshold = float(inpThreshold.value())/100.
                self.threshold = float(inpThreshold.value())
                self.width = int(inpWidth.value())
            except:
                self.threshold = 0.66
                self.width = 9
            if not (self.width%2):
                self.width += 1
            if buttonActive.isChecked():
                if DEBUG:
                    print('ActiveChecked')
                self.removeSpikesActive()
            if buttonAll.isChecked():
                if DEBUG:
                    print('AllChecked')
                self.removeSpikesAll()

    def removeSpikesAll(self):
        self.medianThresholdFilter(False, self.threshold, self.width)

    def removeSpikesActive(self):
        self.medianThresholdFilter(True, self.threshold, self.width)

    def medianThresholdFilter(self, activeOnly, threshold, length):
        if activeOnly:
            active = self._plotWindow.getActiveCurve()
            if not active:
                return
            else:
                x, y, legend, info = active
                self.removeCurve(legend)
                spectra = [active]
        else:
            spectra = self._plotWindow.getAllCurves()
        for (idx, spec) in enumerate(spectra):
            x, y, legend, info = spec
            filtered = medfilt1d(y, length)
            diff = filtered-y
            mean = diff.mean()
            sigma = (x-mean)**2
            sigma = numpy.sqrt(sigma.sum()/len(sigma))
            ynew = numpy.where(abs(diff) > threshold * sigma, filtered, y)
            legend = info.get('selectionlegend','') + ' SR'
            if (idx==0) and (len(spectra)!=1):
                self.addCurve(x,ynew,legend,info,replace=True) 
            else:
                self.addCurve(x,ynew,legend,info)
        

MENU_TEXT = "Remove spikes"
def getPlugin1DInstance(plotWindow,  **kw):
    ob = RemoveSpikes(plotWindow)
    return ob
    
if __name__ == "__main__":
    from PyMca import ScanWindow
    from PyMca import PyMcaQt as qt
    import numpy
    app = qt.QApplication([])
    
    sw = ScanWindow.ScanWindow()

    x = numpy.arange(1000.)
    y0 =  10 * x + 10000. * numpy.exp(-0.5*(x-500)*(x-500)/400) + 1500 * numpy.random.random(1000.)
    y1 =  10 * x + 10000. * numpy.exp(-0.5*(x-600)*(x-600)/400) + 1500 * numpy.random.random(1000.)
    y2 =  10 * x + 10000. * numpy.exp(-0.5*(x-400)*(x-400)/400) + 1500 * numpy.random.random(1000.)
    y3 =  10 * x + 10000. * numpy.exp(-0.5*(x-700)*(x-700)/400) + 1500 * numpy.random.random(1000.)
    y0[320:322] = 50000.
    y2[400:405] = 12300.
    y2[200:205] = 10400.
    y1[620:623] = 16300.
    y1[800:803] = 50000.
    y3[664:666] = 16950.
    y3[699:701] = 20000.
    y3[730:733] = 20000.
    
    
    sw.addCurve(x, y0, legend="Curve0")
    sw.addCurve(x, y1, legend="Curve1")
    sw.addCurve(x, y2, legend="Curve2")
    sw.addCurve(x, y3, legend="Curve3")
    
    plugin = getPlugin1DInstance(sw)    
    plugin.configureFilter()
    
    sw.show()
    app.exec_()
