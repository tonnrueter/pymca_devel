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
__author__ = "V.A. Sole - ESRF Data Analysis"
"""

A Stack plugin is a module that will be automatically added to the PyMca
stack windows in order to perform user defined operations on the data stack.

These plugins will be compatible with any stack window that provides the
functions:
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
    from PyMca import CalculationThread
except ImportError:
    #python 3 might need this
    from . import StackPluginBase
    from . import CalculationThread
try:
    from PyMca.NNMAWindow import NNMAParametersDialog
    from PyMca import StackPluginResultsWindow
    import PyMca.PyMca_Icons as PyMca_Icons
except ImportError:
    print("NNMAStackPlugin importing from somewhere else")
    from NNMAWindow import NNMAParametersDialog
    import StackPluginResultsWindow
    import PyMca_Icons

qt = StackPluginResultsWindow.qt
DEBUG = 0


class NNMAStackPlugin(StackPluginBase.StackPluginBase):
    def __init__(self, stackWindow, **kw):
        StackPluginBase.DEBUG = DEBUG
        StackPluginBase.StackPluginBase.__init__(self, stackWindow, **kw)
        self.methodDict = {'Calculate': [self.calculate,
                                         "Perform NNMA",
                                         None],
                           'Show': [self._showWidget,
                                    "Show last results",
                                    PyMca_Icons.brushselect]}
        self.__methodKeys = ['Calculate', 'Show']
        self.configurationWidget = None
        self.widget = None
        self.thread = None

    def stackUpdated(self):
        if DEBUG:
            print("NNMAStackPlugin.stackUpdated() called")
        self.configurationWidget = None
        self.widget = None

    def selectionMaskUpdated(self):
        if self.widget is None:
            return
        if self.widget.isHidden():
            return
        mask = self.getStackSelectionMask()
        self.widget.setSelectionMask(mask)

    def mySlot(self, ddict):
        if DEBUG:
            print("mySlot ", ddict['event'], ddict.keys())
        if ddict['event'] == "selectionMaskChanged":
            self.setStackSelectionMask(ddict['current'])
        elif ddict['event'] == "addImageClicked":
            self.addImage(ddict['image'], ddict['title'])
        elif ddict['event'] == "removeImageClicked":
            self.removeImage(ddict['title'])
        elif ddict['event'] == "replaceImageClicked":
            self.replaceImage(ddict['image'], ddict['title'])
        elif ddict['event'] == "resetSelection":
            self.setStackSelectionMask(None)

    #Methods implemented by the plugin
    def getMethods(self):
        if self.widget is None:
            return [self.__methodKeys[0]]
        else:
            return self.__methodKeys

    def getMethodToolTip(self, name):
        return self.methodDict[name][1]

    def getMethodPixmap(self, name):
        return self.methodDict[name][2]

    def applyMethod(self, name):
        return self.methodDict[name][0]()

    #The specific part
    def calculate(self):
        if self.configurationWidget is None:
            self.configurationWidget = NNMAParametersDialog(None)
            self._status = qt.QLabel(self.configurationWidget)
            self._status.setAlignment(qt.Qt.AlignHCenter)
            font = qt.QFont(self._status.font())
            font.setBold(True)
            self._status.setFont(font)
            self._status.setText("Ready")
            self.configurationWidget.layout().addWidget(self._status)
        activeCurve = self.getActiveCurve()
        if activeCurve is None:
            #I could get some defaults from the stack itslef
            raise ValueError("Please select an active curve")
            return
        x, spectrum, legend, info = activeCurve
        spectrumLength = max(spectrum.shape)
        oldValue = self.configurationWidget.nPC.value()
        self.configurationWidget.nPC.setMaximum(spectrumLength)
        self.configurationWidget.nPC.setValue(min(oldValue, spectrumLength))
        binningOptions = [1]
        for number in [2, 3, 4, 5, 7, 9, 10, 11, 13, 15, 17, 19]:
            if (spectrumLength % number) == 0:
                binningOptions.append(number)
        ddict = {'options': binningOptions,
                 'binning': 1,
                 'method': 0}
        self.configurationWidget.setParameters(ddict)
        #y = spectrum
        #self.configurationWidget.setSpectrum(x, y)
        self.configurationWidget.show()
        self.configurationWidget.raise_()
        ret = self.configurationWidget.exec_()
        if ret:
            self._executeFunctionAndParameters()

    def _executeFunctionAndParameters(self):
        if DEBUG:
            print("NNMAStackPlugin _executeFunctionAndParameters")
        self.widget = None
        self.thread = CalculationThread.CalculationThread(\
                            calculation_method=self.actualCalculation)
        qt.QObject.connect(self.thread,
                     qt.SIGNAL('finished()'),
                     self.threadFinished)
        self.configurationWidget.show()
        message = "Please wait. NNMA Calculation going on."
        if DEBUG:
            print("NNMAStackPlugin starting thread")
        self.thread.start()
        if DEBUG:
            print("NNMAStackPlugin waitingMessageDialog")
        CalculationThread.waitingMessageDialog(self.thread,
                                message=message,
                                parent=self.configurationWidget)
        if DEBUG:
            print("NNMAStackPlugin waitingMessageDialog passed")

    def actualCalculation(self):
        if DEBUG:
            print("NNMAStackPlugin actualCalculation")
        nnmaParameters = self.configurationWidget.getParameters()
        self._status.setText("Calculation going on")
        self.configurationWidget.setEnabled(False)
        #self.configurationWidget.close()
        #At some point I should make sure I get directly the
        #function and the parameters from the configuration widget
        function = nnmaParameters['function']
        ddict = {}
        ddict.update(nnmaParameters['kw'])
        ddict['ncomponents'] = nnmaParameters['npc']
        ddict['binning'] = nnmaParameters['binning']
        #ddict['kmeans'] = False
        del nnmaParameters
        stack = self.getStackDataObject()
        if isinstance(stack, numpy.ndarray):
            if stack.data.dtype not in [numpy.float, numpy.float32]:
                print("WARNING: Non floating point data")
                text = "Calculation going on."
                text += " WARNING: Non floating point data."
                self._status.setText(text)

        oldShape = stack.data.shape
        result = function(stack, **ddict)
        if stack.data.shape != oldShape:
            stack.data.shape = oldShape
        return result

    def threadFinished(self):
        if DEBUG:
            print("NNMAStackPlugin threadFinished")
        result = self.thread.result
        self.thread = None
        if type(result) == type((1,)):
            #if we receive a tuple there was an error
            if len(result):
                if result[0] == "Exception":
                    self._status.setText("Ready after calculation error")
                    self.configurationWidget.setEnabled(True)
                    raise Exception(result[1], result[2])
                    return
        self._status.setText("Ready")
        self.configurationWidget.setEnabled(True)
        self.configurationWidget.close()

        images, eigenValues, eigenVectors = result
        imageNames = None
        vectorNames = None
        nimages = images.shape[0]
        imageNames = []
        vectorNames = []
        vectorTitles = []
        for i in range(nimages):
            imageNames.append("NNMA Image %02d" % i)
            vectorNames.append("NNMA Component %02d" % i)
            vectorTitles.append("%g %% explained intensity" %\
                                               eigenValues[i])
        if DEBUG:
            print("NNMAStackPlugin threadFinished. Create widget")
        self.widget = StackPluginResultsWindow.StackPluginResultsWindow(\
                                        usetab=True)
        if DEBUG:
            print("NNMAStackPlugin threadFinished. Widget created")
        self.widget.buildAndConnectImageButtonBox()
        qt = StackPluginResultsWindow.qt
        qt.QObject.connect(self.widget,
                           qt.SIGNAL('MaskImageWidgetSignal'),
                           self.mySlot)

        self.widget.setStackPluginResults(images,
                                          spectra=eigenVectors,
                                          image_names=imageNames,
                                          spectra_names=vectorNames,
                                          spectra_titles=vectorTitles)
        self._showWidget()

    def _showWidget(self):
        if self.widget is None:
            return
        #Show
        self.widget.show()
        self.widget.raise_()

        #update
        self.selectionMaskUpdated()

MENU_TEXT = "PyMca NNMA"


def getStackPluginInstance(stackWindow, **kw):
    ob = NNMAStackPlugin(stackWindow)
    return ob
