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
__author__ = "T. Rueter - ESRF Data Analysis"

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
    
DEBUG = True
class RemoveSpikes(Plugin1DBase.Plugin1DBase):
    def __init__(self,  plotWindow,  **kw):
        Plugin1DBase.Plugin1DBase.__init__(self,  plotWindow,  **kw)
        self.methodDict = {
            'Apply to active curve':
                [self.removeSpikesActive,
                 'Remove spikes from active curve',
                 None],
            'Apply to all curves':
                [self.removeSpikesAll,
                 'Remove spikes from all curves',
                 None],
            'Configure median filter':
                [self.configureFilter,
                 'Set threshold and width of the filter',
                 None]
        }
        # TODO: Set threshold and width default values
        print 'Here come the mD:\n',self.methodDict
    
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
        # TODO: Spawn QDialog, ask for threshold % window width
        msg = qt.QDialog()
        msgLayout = qt.QGridLayout()
        buttonLayout = qt.QHBoxLayout()
        threshold = qt.QSpinBox()
        threshold.setRange(1,99)
        threshold.setSingleStep(2)
        threshold.setValue(66)
##        threshold.setValidator(qt.QDoubleValidator(0.,0.99, 2, threshold))
        width = qt.QSpinBox()
        width.setRange(1,100)
        width.setSingleStep(2)
        width.setValue(9)
        labelWidth = qt.QLabel('Width (must be odd)')
        labelThreshold = qt.QLabel('Threshold (Percent of Maximum)')
        buttonOK = qt.QPushButton('Ok')
        buttonOK.clicked.connect(msg.accept)
        buttonCancel = qt.QPushButton('Cancel')
        buttonCancel.clicked.connect(msg.reject)
        buttonLayout.addWidget(qt.HorizontalSpacer())
        buttonLayout.addWidget(buttonOK)
        buttonLayout.addWidget(buttonCancel)
        msgLayout.addWidget(labelWidth,0,0)
        msgLayout.addWidget(width,0,1)
        msgLayout.addWidget(labelThreshold,1,0)
        msgLayout.addWidget(threshold,1,1)
        msgLayout.addLayout(buttonLayout,2,0,1,2)
        msg.setLayout(msgLayout)
        val = msg.exec_()
        if val:
            try:
                thr = float(str(threshold.text()))
                wid = float(str(width.value()))
            except:
                thr = 0.66
                wid = 9
        print thr, wid

    def removeSpikesAll(self):
        self.medianThresholdFilter(False)

    def removeSpikesActive(self):
        self.medianThresholdFilter(True)

    def medianThresholdFilter(self, activeOnly, threshold=None, length=9):
        if activeOnly:
            active = self._plotWindow.getActiveCurve()
            if not active:
                return
            else:
                x, y, legend, info = spec
                self.removeCurve(legend)
                spectra = [active]
        else:
            spectra = self._plotWindow.getAllCurves()
        for (i,spec) in enumerate(spectra):
            x, y, legend, info = spec
            filtered = medfilt1d(y, length)
            diff = abs(filtered-y)
            if not threshold:
                threshold = .66 * diff.max() # OR: diff.mean()
            ynew = numpy.where(diff>threshold, filtered, y)
            legend = info.get('selectionlegend','') + ' SR'
            if (i==0) and (len(spectra)!=1):
                self.addCurve(x,ynew,legend,info, replace=True)
            else:
                self.addCurve(x,ynew,legend,info)
        

MENU_TEXT = "Spike-X Plugin" #Remove spikes
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
    y2[320:322] = 50000.
    info0 = {'FileHeader':['#F /data/id08/inhouse/somerandomname'],'xlabel': 'foo', 'ylabel': 'arb', 'MotorNames': 'oxPS Motor11 Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3',  'MotorValues': '1 8.69271399699 21.9836418539 0.198068826612 0.484475455792 0.350252217264 0.663925270933 0.813033264421 0.221149410218 0.593188258866 0.678010392881 0.267389247833 0.677890617858'}
    info1 = {'MotorNames': 'PhaseD oxPS Motor16 Motor15 Motor14 Motor13 Motor12 Motor11 Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3', 'MotorValues': '0.470746882688 0.695816070299 0.825780811755 0.25876374531 0.739264467436 0.090842892619 2 0.213445659833 0.823400550314 0.020278096857 0.568744021322 0.85378115537 0.696730386891 0.269196313956 0.793293334395 0.769216567757 0.959092709527 0.0109264683697 0.538264972553'}
    info2 = {'MotorNames': 'PhaseD oxPS Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3',  'MotorValues': '2 0.44400576644 0.613870067852 0.901968648111 0.319768771085 0.571432278628 0.278675836163 0.154436774878 0.416231999332 0.294201017231 0.813913587748 0.577572903105 0.869045182568'}
    
    sw.addCurve(x, y0, legend="Curve0", info=info0, replot=False, replace=False)
    sw.addCurve(x, y1, legend="Curve1", info=info1, replot=False, replace=False)
    sw.addCurve(x, y2, legend="Curve2", info=info2, replot=False, replace=False)

    plugin = getPlugin1DInstance(sw)
    plugin.configureFilter()
##    print plugin.methodDict
##    plugin.applyMethod(plugin.getMethods()[0])
    
##    sw.show()
    
    app.exec_()
    
