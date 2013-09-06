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
__author__ = "V.A. Sole - ESRF Software Group"
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
    from PyMca import CalculationThread
except ImportError:
    print("XASStackNormalizationPlugin importing bases from somewhere else")
    from . import StackPluginBase
    from . import CalculationThread

try:
    from PyMca import PyMcaQt as qt
    from PyMca import StackPluginResultsWindow
    from PyMca import XASNormalization
    from PyMca import XASNormalizationWindow
except ImportError:
    print("XASStackNormalizationPlugin importing from somewhere else")
    import PyMcaQt as qt
    import StackPluginResultsWindow
    import XASNormalization
    import XASNormalizationWindow

DEBUG = 0

class XASStackNormalizationPlugin(StackPluginBase.StackPluginBase):
    def __init__(self, stackWindow, **kw):
        StackPluginBase.DEBUG = DEBUG
        StackPluginBase.StackPluginBase.__init__(self, stackWindow, **kw)
        self.methodDict = {}
        text = "Replace current stack by a normalized one."
        function = self.XASNormalize
        info = text
        icon = None
        self.methodDict["XANES Normalization"] =[function,
                                                 info,
                                                 icon]

        self.__methodKeys = ["XANES Normalization"]
        self.widget = None
        self.imageWidget = None
        
    #Methods implemented by the plugin
    def stackUpdated(self):
        if self.widget is not None:
            self.widget.close()
        self.widget = None

    def selectionMaskUpdated(self):
        if self.imageWidget is None:
            return
        if self.imageWidget.isHidden():
            return
        mask = self.getStackSelectionMask()
        self.imageWidget.setSelectionMask(mask)

    def getMethods(self):
        return self.__methodKeys

    def getMethodToolTip(self, name):
        return self.methodDict[name][1]

    def getMethodPixmap(self, name):
        return self.methodDict[name][2]

    def applyMethod(self, name):
        return self.methodDict[name][0]()

    
    # own stuff
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
    
    def XASNormalize(self):
        stack = self.getStackDataObject()
        if not isinstance(stack.data, numpy.ndarray):
            text = "This method does not work with dynamically loaded stacks"
            raise TypeError(text)
        activeCurve = self.getActiveCurve()
        if activeCurve is None:
            return
        x, spectrum, legend, info = activeCurve
        if self.widget is None:
            self.widget = XASNormalizationWindow.XASNormalizationDialog(None,
                                                spectrum, energy=x)
        else:
            oldParameters = self.widget.getParameters()
            oldEnergy = self.widget.parametersWidget.energy
            oldEMin = oldEnergy.min()
            oldEMax = oldEnergy.max()
            self.widget.setData(spectrum, energy=x)
            if abs(oldEMin - x.min()) < 1:
                if abs(oldEMax - x.max()) < 1:
                    self.widget.setParameters(oldParameters)
        ret = self.widget.exec_()
        if ret:
            parameters = self.widget.getParameters()
            # TODO: this dictionary adaptation should be made
            #       by the configuration
            if parameters['auto_edge']:
                edge = None
            else:
                edge = parameters['edge_energy']
            energy = x
            pre_edge_regions = parameters['pre_edge']['regions']
            post_edge_regions = parameters['post_edge']['regions']
            algorithm ='polynomial'
            algorithm_parameters = {}
            algorithm_parameters['pre_edge_order'] = parameters['pre_edge']\
                                                             ['polynomial']
            algorithm_parameters['post_edge_order'] = parameters['post_edge']\
                                                             ['polynomial']

            result  = self.__replaceStackByXASNormalizedData(stack,
                                            energy=energy,
                                            edge=edge,
                                            pre_edge_regions=pre_edge_regions,
                                            post_edge_regions=post_edge_regions,
                                            algorithm=algorithm,
                                            algorithm_parameters=algorithm_parameters)
            if result[0] == 'Exception':
                # exception occurred
                raise Exception(result[1], result[2], result[3])
            else:
                edges, jumps, errors = result
            images, names = self.getStackROIImagesAndNames()
            edges.shape = images[0].shape
            jumps.shape = images[0].shape
            errors.shape = images[0].shape
            self.setStack(stack)
            if self.imageWidget is None:
                self.imageWidget = StackPluginResultsWindow.StackPluginResultsWindow(\
                                        usetab=False)
                self.imageWidget.buildAndConnectImageButtonBox()
                qt = StackPluginResultsWindow.qt
                qt.QObject.connect(self.imageWidget,
                           qt.SIGNAL('MaskImageWidgetSignal'),
                           self.mySlot)
                self.methodDict["Show Images"] =[self._showImageWidget,
                                                 "Show calculated jump and edge position images",
                                                 None]
                self.__methodKeys.append("Show Images")
                self.imageWidget.setStackPluginResults([jumps, errors, edges],
                                                        image_names=['Jump',
                                                                     'Errors',
                                                                     'Edge Position'])
            self._showImageWidget()

    def __replaceStackByXASNormalizedData(self, *var, **kw):
        self._progress = 0.0
        thread = CalculationThread.CalculationThread(\
                calculation_method=self.replaceStackByXASNormalizedData,
                calculation_vars=var,
                calculation_kw=kw)
        thread.start()
        CalculationThread.waitingMessageDialog(thread,
                                               message="Please wait. Calculation going on.",
                                               parent=self.widget,
                                               modal=True,
                                               update_callback=self._waitingCallback)
        return thread.result

    def _waitingCallback(self):
        ddict = {}
        ddict['message'] = "Calculation Progress = %d %%" % self._progress
        return ddict

    def _showImageWidget(self):
        if self.imageWidget is None:
            return
        #Show
        self.imageWidget.show()
        self.imageWidget.raise_()

        #update
        self.selectionMaskUpdated()

    def replaceStackByXASNormalizedData(self,
                                        stack,
                                        energy=None,
                                        edge=None,
                                        pre_edge_regions=None,
                                        post_edge_regions=None,
                                        algorithm='polynomial',
                                        algorithm_parameters=None):
        """
        Performs an in place replacement of a set of spectra by a set of
        normalized spectra.
        """
        mcaIndex = -1
        if hasattr(stack, "info") and hasattr(stack, "data"):
            actualData = stack.data
            mcaIndex = stack.info.get('McaIndex', -1)
        else:
            actualData = stack
        if not isinstance(actualData, numpy.ndarray):
            raise TypeError("Currently this method only supports numpy arrays")

        # Take a data view
        oldShape = actualData.shape
        data = actualData[:]
        DONE = 0
        if mcaIndex in [-1, len(data.shape)-1]:
            data.shape = -1, oldShape[-1]
            edges = numpy.zeros(data.shape[0], numpy.float32)
            jumps = numpy.zeros(data.shape[0], numpy.float32)
            errors = numpy.zeros(data.shape[0], numpy.float32)
            total = 0.01 * data.shape[0]
            for i in range(data.shape[0]):
                self._progress = i / total
                try:
                    ene, spe, ed, jmp = XASNormalization.XASNormalization(data[i,:],
                                energy=energy,
                                edge=edge,
                                pre_edge_regions=pre_edge_regions,
                                post_edge_regions=post_edge_regions,
                                algorithm=algorithm,
                                algorithm_parameters=algorithm_parameters)[0:4]
                except:
                    # what to do?
                    # for the data is clear (set to 0)
                    # for the jump 0 is also a good compromise
                    # for the edge?
                    data[i, :] = 0
                    errors[i] = 1
                    jumps[i] = 0
                    edges[i] = 0
                    continue
                if not DONE:
                    c0 = (numpy.nonzero(energy >= (ed + pre_edge_regions[0][0]))[0]).min()
                    c1 = (numpy.nonzero(energy <= (ed + post_edge_regions[-1][1]))[-1]).max()
                    DONE = True
                if ((spe.max()-spe.min()) > 10.) or (jmp < 0):
                    data[i, :] = 0.0
                    # should I give some useless values?
                    edges[i] = 0.0
                    jumps[i] = 0.0
                elif 0:
                    # this approach removed
                    data[i,:c0] = spe[c0]
                    data[i, c0:c1] = spe[c0:c1]
                    data[i, c1:] = spe[c1]
                    edges[i] = ed
                    jumps[i] = jmp
                else:
                    # it seems more appropriate to set the channels below and
                    # above limits to 0 than to the corresponding limits of the region
                    data[i,:c0] = 0.0
                    data[i, c0:c1] = spe[c0:c1]
                    data[i, c1:] = 0.0
                    edges[i] = ed
                    jumps[i] = jmp
            self._progress = 100
            data.shape = oldShape
        elif mcaIndex == 0:
            data.shape = oldShape[0], -1
            edges = numpy.zeros(data.shape[-1], numpy.float32)
            jumps = numpy.zeros(data.shape[-1], numpy.float32)
            errors = numpy.zeros(data.shape[-1], numpy.float32)
            total = 0.01 * data.shape[-1]
            for i in range(data.shape[-1]):
                self._progress = i / total 
                try:
                    ene, spe, ed, jmp = XASNormalization.XASNormalization(data[:, i],
                              energy=energy,
                              edge=edge,
                              pre_edge_regions=pre_edge_regions,
                              post_edge_regions=post_edge_regions,
                              algorithm=algorithm,
                              algorithm_parameters=algorithm_parameters)[0:4]
                except:
                    # what to do?
                    # for the data is clear (set to 0)
                    # for the jump 0 is also a good compromise
                    # for the edge?
                    data[:, i] = 0
                    jumps[i] = 0
                    edges[i] = 0
                    errors[i] = 1
                    continue
                if not DONE:
                    c0 = (numpy.nonzero(energy >= (ed + pre_edge_regions[0][0]))[0]).min()
                    c1 = (numpy.nonzero(energy <= (ed + post_edge_regions[-1][1]))[-1]).max()
                    DONE = True
                if ((spe.max()-spe.min()) > 10.) or (jmp < 0):
                    data[:, i] = 0.0
                    # should I give some useless values?
                    edges[i] = 0.0
                    jumps[i] = 0.0
                else:
                    # it seems more appropriate to set the channels below and
                    # above limits to 0 than to the corresponding limits of the region
                    data[:c0, i] = 0.0
                    data[c0:c1, i] = spe[c0:c1]
                    data[c1:, i] = 0.0
                    edges[i] = ed
                    jumps[i] = jmp
            self._progress = 100
            data.shape = oldShape
        else:
            raise ValueError("Invalid 1D index %d" % mcaIndex)
        return edges, jumps, errors


MENU_TEXT = "XAS Stack Normalization"
def getStackPluginInstance(stackWindow, **kw):
    ob = XASStackNormalizationPlugin(stackWindow)
    return ob
