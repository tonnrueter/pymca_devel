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
from PyMca import PyMcaQt as qt
from PyMca import PyMcaDirs, PyMcaFileDialogs
from PyMca import ConfigDict

try:
    from PyMca import Plugin1DBase
except ImportError:
    print("WARNING:AlignmentScanPlugin import from somewhere else")
    from . import Plugin1DBase

from PyMca import SpecfitFuns

#class TableItem(qt.QTableWidgetItem):
#    def __init__(self, parent, content):
        
DEBUG = True
class AlignmentWidget(qt.QDialog):
    def __init__(self, parent, ddict, plugin):
#    def __init__(self, parent, plugin):
        qt.QDialog.__init__(self, parent)
        self.setWindowTitle('FFT Shifts')
        
        nCols = 2
        nRows = len(ddict)
        self.plugin = plugin
        
        # Buttons
        buttonSave   = qt.QPushButton('Save')
        buttonLoad   = qt.QPushButton('Load')
        buttonStore   = qt.QPushButton('Store')
        buttonApply  = qt.QPushButton('Apply')
        buttonCancel = qt.QPushButton('Cancel')
        buttonCalc = qt.QPushButton('Cancel')

        # Table
        self.shiftTab = qt.QTableWidget(nRows, nCols)
        self.shiftTab.verticalHeader().hide()
        self.shiftTab.setHorizontalHeaderLabels(['Legend','Shift'])
        self.shiftTab.horizontalHeader().setStretchLastSection(True)
        
        # Shift Method selector
        self.shiftMethodComboBox = qt.QComboBox()
        self.shiftMethodComboBox.addItems(
            ['FFT',
             'X'])
             
        # Alignment Method selector
        self.alignmentMethodComboBox = qt.QComboBox()
        self.alignmentMethodComboBox.addItems(
            ['FFT',
             'MAX'])
        
        # Fill table with data
#        self.setDict(plugin.shiftDict)
        self.setDict(ddict)
        
        #Layouts
        topLayout = qt.QHBoxLayout()
        topLayout.addWidget(buttonCalc)
        topLayout.addWidget(qt.HorizontalSpacer())
        topLayout.addWidget(qt.QLabel('Alignment method:'))
        topLayout.addWidget(self.alignmentMethodComboBox)
        topLayout.addWidget(qt.QLabel('Shift method:'))
        topLayout.addWidget(self.shiftMethodComboBox)
        
        buttonLayout = qt.QHBoxLayout()
        buttonLayout.addWidget(buttonSave)
        buttonLayout.addWidget(buttonLoad)
        buttonLayout.addWidget(qt.HorizontalSpacer())
        buttonLayout.addWidget(buttonApply)
        buttonLayout.addWidget(buttonStore)
        buttonLayout.addWidget(buttonCancel)
        
        mainLayout = qt.QVBoxLayout()
        mainLayout.addLayout(topLayout)
        mainLayout.addWidget(self.shiftTab)
        mainLayout.addLayout(buttonLayout)
        self.setLayout(mainLayout)
        
        # Connects
        self.shiftTab.cellChanged.connect(self.validateInput)
        buttonApply.clicked.connect(self.accept)
        buttonCancel.clicked.connect(self.reject)
        buttonStore.clicked.connect(self.store)
        buttonSave.clicked.connect(self.saveDict)
        buttonLoad.clicked.connect(self.loadDict)
        # ..to Plugin instance
        buttonCalc.clicked.connect(self.recalculate)
        self.alignmentMethodComboBox.currentIndexChanged['QString'].connect(self.plugin.setAlignmentMethod)
        self.shiftMethodComboBox.currentIndexChanged['QString'].connect(self.plugin.setShiftMethod)

    def recalculate(self):
#        self.setDict(self.plugin.calculateShifts())
        ddict = self.plugin.calculateShifts()
        print ddict
        self.setDict(ddict)

    def store(self):
        self.done(2)

    def loadDict(self):
        openDir = PyMcaDirs.outputDir
        filter = 'PyMca (*.shift)'
        filename = qt.QFileDialog.\
                    getOpenFileName(self,
                                    'Load Shifts obtained from FFTAlignment',
                                    openDir,
                                    filter)
        inDict = ConfigDict.ConfigDict()
        try:
            inDict.read(filename)
        except IOError:
            msg = qt.QMessageBox()
            msg.setTitle('FFTAlignment Load Error')
            msg.setText('Unable to read shifts form file \'%s\''%filename)
            return
        if 'Shifts' not in inDict:
            return
        try:
            self.setDict(inDict['Shifts'])
        except:
            msg = qt.QMessageBox()
            msg.setWindowTitle('FFTAlignment Load Error')
            msg.setText('Configuration file \'%s\' corruted' % filename)
            msg.exec_()

    def saveDict(self):
        saveDir = PyMcaDirs.outputDir
        filter = ['PyMca (*.shift)']
        try:
            filename = PyMcaFileDialogs.\
                        getFileList(parent=self,
                            filetypelist=filter,
                            message='Safe FFT Alignment shifts',
                            mode='SAVE',
                            single=True)[0]
        except IndexError:
            # Returned list is empty
            return
        if DEBUG:
            print('saveOptions -- Filename: "%s"' % filename)
        if len(filename) == 0:
            return False
        if not str(filename).endswith('.shift'):
            filename += '.shift'
        outDict = ConfigDict.ConfigDict()
        outDict['Shifts'] = self.getDict()
        try:
            outDict.write(filename)
        except IOError:
            msg = qt.QMessageBox()
            msg.setWindowTitle('FFTAlignment Save Error')
            msg.setText('Unable to write configuration to \'%s\''%filename)
            msg.exec_()
        return True

    def getAlignmentMethodName(self):
        return self.alignmentMethodComboBox.currentText()

    def getShiftMethodName(self):
        return self.shiftMethodComboBox.currentText()

    def validateInput(self, row, col):
        if col == 0:
            return
        elif col == 1:
            item  = self.shiftTab.item(row, 1)
            try:
                floatValue = float(item.text())
                item.setText('%.4f'%floatValue)
            except:
                floatValue = float('NaN')
                item.setText(str(floatValue))

    def getDict(self):
        ddict = {}
        for idx in range(self.shiftTab.rowCount()):
            legend = self.shiftTab.item(idx, 0)
            value  = self.shiftTab.item(idx, 1)
            try:
                floatValue = float(value.text())
            except:
                floatValue = float('NaN')
            ddict[str(legend.text())] = floatValue
#        print '-- getDict : ', ddict
        return ddict
    
    def setDict(self, ddict):
        self.shiftTab.clear()
        dkeys = sorted(ddict.keys())
        dvals = ['%.4f'%ddict[k] for k in dkeys]
        for j, dlist in enumerate([dkeys, dvals]):
            for i in range(len(dlist)):
                if j == 0:
                    elem = qt.QTableWidgetItem(dlist[i])
                    elem.setFlags(qt.Qt.ItemIsEnabled)
                elif j == 1:
                    elem = qt.QTableWidgetItem(str(dlist[i]))
                    elem.setTextAlignment(qt.Qt.AlignRight)
                    elem.setTextAlignment(qt.Qt.AlignRight + qt.Qt.AlignVCenter)
                    elem.setFlags(qt.Qt.ItemIsEditable | qt.Qt.ItemIsEnabled)
                else:
                    elem = qt.QTableWidgetItem('')
                self.shiftTab.setItem(i,j, elem)
        self.shiftTab.resizeColumnToContents(0)    
        self.shiftTab.resizeRowsToContents()

DEBUG = False
class AlignmentScanPlugin(Plugin1DBase.Plugin1DBase):
    def __init__(self, plotWindow, **kw):
        Plugin1DBase.Plugin1DBase.__init__(self, plotWindow, **kw)
        self.__randomization = True
        self.__methodKeys = []
        self.methodDict = {}

        function = self.calculateAndApplyShifts
        method = "Perform FFT Alignment"
        text = "Perform FFT based alignment and shift"
        info = text
        icon = None
        self.methodDict[method] = [function,
                                   info,
                                   icon]
        self.__methodKeys.append(method)

        function = self.applyShifts
        method = "Apply FFT Alignment"
        text = "Apply calculated shifts to curves"
        info = text
        icon = None
        self.methodDict[method] = [function,
                                   info,
                                   icon]
        self.__methodKeys.append(method)

        function = self.showShifts
        method = "Show shifts"
        text = "Calculates shifts by FFT alignment and displays them"
        info = text
        icon = None
        self.methodDict[method] = [function,
                                   info,
                                   icon]
        self.__methodKeys.append(method)

        function = self.calculateAndShowShifts
        method = "Calculate shifts"
        text = "Displays shifts previously calculated by FFT alignment"
        info = text
        icon = None
        self.methodDict[method] = [function,
                                   info,
                                   icon]
        self.__methodKeys.append(method)
        
        self.shiftDict = {}
        self.shiftMethod = self.fftShift
        self.alignmentMethod = self.calculateShiftsFFT
        
    #Methods to be implemented by the plugin
    def getMethods(self, plottype=None):
        """
        A list with the NAMES  associated to the callable methods
        that are applicable to the specified plot.

        Plot type can be "SCAN", "MCA", None, ...        
        """
#        if self.__randomization:
#            return self.__methodKeys[0:1] +  self.__methodKeys[2:]
#        else:
#            return self.__methodKeys[1:]
        return self.__methodKeys

    def getMethodToolTip(self, name):
        """
        Returns the help associated to the particular method name or None.
        """
        return self.methodDict[name][1]

    def getMethodPixmap(self, name):
        """
        Returns the pixmap associated to the particular method name or None.
        """
        return None

    def applyMethod(self, name):
        """
        The plugin is asked to apply the method associated to name.
        """
        return self.methodDict[name][0]()

    def dummy(self):
        print 'Dummy'

    def calculateAndApplyShifts(self):
        self.calculateShifts()
        self.applyShifts()
        self.shiftDict = {}

    def calculateAndShowShifts(self):
        self.calculateShifts()
        self.showShifts()

    def calculateShifts(self):
        self.shiftDict = self.alignmentMethod()
        return self.shiftDict

    def calculateShiftsMax(self):
        retDict = {}
        
        curves = self.getAllCurves()
        nCurves = len(curves)
        
        if nCurves < 2:
            raise ValueError("At least 2 curves needed")
            return
        
        # Check if plotwindow is zoomed in
        xmin, xmax = self.getGraphXLimits()
        # Determine largest overlap between curves
        xmin0, xmax0 = self.getXLimits(x for (x,y,leg,info) in curves)
        if xmin0 > xmin:
            xmin = xmin0
        if xmax0 < xmax:
            xmax = xmax0

        # Get active curve            
        activeCurve = self.getActiveCurve()
        if activeCurve is None:
            # If active curve is not set, continue with first curve
            activeCurve = curves[0]
        else:
            activeLegend = activeCurve[2]
            idx = list.index([curve[2] for curve in curves],
                             activeLegend)
            activeCurve = curves[idx]
        
        x0, y0 = activeCurve[0], activeCurve[1]
        idx = numpy.nonzero((xmin <= x0) & (x0 <= xmax))[0]
        x0 = numpy.take(x0, idx)
        y0 = numpy.take(y0, idx)
        
        # Determine the index of maximum in active curve
        shift0 = numpy.argmax(y0)
        for x,y,legend,info in curves:
            shifty = numpy.argmax(y)
            shift = x0[shift0] - x[shifty]
            key = info['selectionlegend']
            retDict[key] = shift
        return retDict
        
    def calculateShiftsFFT(self, threshold=.95):
        retDict = {}
        
        curves = self.interpolate()
        nCurves = len(curves)
        if nCurves < 2:
            raise ValueError("At least 2 curves needed")
            return
        
        # Check if scan window is zoomed in
        xmin, xmax = self.getGraphXLimits()
        # Determine largest overlap between curves
        xmin0, xmax0 = self.getXLimits(x for (x,y,leg,info) in curves)
        if xmin0 > xmin:
            xmin = xmin0
        if xmax0 < xmax:
            xmax = xmax0

        # Get active curve            
        activeCurve = self.getActiveCurve()
        if activeCurve is None:
            # If active curve is not set, continue with first curve
            activeCurve = curves[0]
        else:
            activeLegend = activeCurve[2]
            idx = list.index([curve[2] for curve in curves],
                             activeLegend)
            activeCurve = curves[idx]

        x0, y0 = activeCurve[0], activeCurve[1]
        idx = numpy.nonzero((xmin <= x0) & (x0 <= xmax))[0]
        x0 = numpy.take(x0, idx)
        y0 = numpy.take(y0, idx)
        
        fft0 = numpy.fft.fft(y0)
        fftList = []
        for x,y,legend,info in curves:
            idx = numpy.nonzero((x >= xmin) & (x <= xmax))[0]
            x = numpy.take(x, idx)
            y = numpy.take(y, idx)
            ffty = numpy.fft.fft(y)
            fftList += [ffty]
            shiftTmp = numpy.fft.ifft(fft0 * ffty.conjugate()).real
            shiftPhase = numpy.zeros(shiftTmp.shape, dtype=shiftTmp.dtype)
            m = shiftTmp.size//2
            shiftPhase[m:] = shiftTmp[:-m]
            shiftPhase[:m] = shiftTmp[-m:]
            
            # Thresholding
            thr = threshold * shiftPhase.max()
            idx = numpy.nonzero(shiftPhase > thr)[0]
            shiftTmp = (shiftPhase[idx] * idx/shiftPhase[idx].sum()).sum()
            shift = (shiftTmp - m) * (x[1] - x[0])
            
            # Remove counter name from legend
            #key = ' '.join(legend.split(' ')[:-1])
            key = info['selectionlegend']
            retDict[key] = shift
        return retDict

    def showShifts(self):
        if len(self.shiftDict) == 0:
            msg = qt.QMessageBox(None)
            msg.setWindowTitle('Alignment Error')
            msg.setText('No shift data present.\nDo you want to calculate shifts now?')
            msg.setStandardButtons(qt.QMessageBox.Ok | qt.QMessageBox.Close)
            if msg.exec_():
                self.calculateShifts()
            else:
                return
        widget = AlignmentWidget(None, self.shiftDict, self)
        ret = widget.exec_()
        if ret == 1:
            # Result code Accepted
            self.shiftDict = widget.getDict()
            self.setShiftMethod(widget.getShiftMethodName())
            self.applyShifts()
        elif ret == 2:    
            # Result code Store
            self.shiftDict = widget.getDict()
            self.setShiftMethod(widget.getShiftMethodName())
        else:
            # Dialog is canceled
            self.shiftDict = {}
#        widget.destroy() # Widget should be destroyed after finishing method
        return

    def setShiftMethod(self, methodName):
        print 'setShiftMethod --',methodName
        methodName = str(methodName)
        if methodName == 'FFT':
            self.shiftMethod = self.fftShift
        elif methodName == 'X':
            self.shiftMethod = self.xShift
        else:
            # Unknown method name, use fftShift as default
            self.shiftMethod = self.fftShift
        print 'setShiftMethod --',self.shiftMethod

    def setAlignmentMethod(self, methodName):
        print 'setAlignmentMethod --',methodName
        methodName = str(methodName)
        if methodName == 'FFT':
            self.alignmentMethod = self.calculateShiftsFFT
        elif methodName == 'MAX':
            self.alignmentMethod = self.calculateShiftsMax
        else:
            # Unknown method name, use fftShift as default
            self.alignmentMethod = self.calculateShiftsFFT
        print 'setAlignmentMethod --',self.alignmentMethod

    def getAllCurves(self):
        '''
        Ensures that the x-range of the curves
        is strictly monotonically increasing.
        Conserves curves legend and info dictionary.
        '''
        curves = Plugin1DBase.Plugin1DBase.getAllCurves(self)

        processedCurves = []
        for curve in curves:
            x, y, legend, info = curve[0:4]
            xproc = x[:]
            yproc = y[:]
            # Sort
            idx = numpy.argsort(xproc, kind='mergesort')
            xproc = numpy.take(xproc, idx)
            yproc = numpy.take(yproc, idx)
            # Ravel, Increasing
            xproc = xproc.ravel()
            idx = numpy.nonzero((xproc[1:] > xproc[:-1]))[0]
            xproc = numpy.take(xproc, idx)
            yproc = numpy.take(yproc, idx)
            processedCurves += [(xproc, yproc, legend, info)]
        return processedCurves

    def interpolate(self, factor=1.):
        '''
        Input
        -----
        factor : float
            factor used to oversample existing data, use
            with caution.
        
        Interpolates all existing curves to an equidistant
        x-range using the either the active or the first
        curve do determine the number of data points.
        
        Returns
        -------
        interpCurves : 
        '''
        curves = self.getAllCurves()
        if len(curves) < 1:
            raise ValueError("At least 1 curve needed")
            if DEBUG:
                print '-- interpolate: no curves present'
            return

        activeCurve = self.getActiveCurve()
        if not activeCurve:
            activeCurve = curves[0]
        else:
            activeLegend = activeCurve[2]
            idx = list.index([curve[2] for curve in curves],
                             activeLegend)
            activeCurve = curves[idx]
        activeX, activeY, activeLegend, activeInfo = activeCurve[0:4]
        
        # Determine average spaceing between Datapoints
        step = numpy.average(numpy.diff(activeX))
        xmin, xmax = self.getXLimits([x for (x,y,leg,info) in curves],
                                     overlap=False)
        num  = factor * numpy.ceil((xmax-xmin)/step)
        
        # Create equidistant x-range, exclude first and last point
        xeq = numpy.linspace(xmin, xmax, num, endpoint=False)[:-1]
        
        # Interpolate on sections of xeq
        interpCurves = []
        for (x,y,legend,info) in curves:
            idx = numpy.nonzero((x.min()<xeq) & (xeq<x.max()))[0]
            xi = numpy.take(xeq, idx)
            yi = SpecfitFuns.interpol([x], y, xi.reshape(-1,1), y.min())
            yi.shape = -1
            interpCurves += [(xi, yi, legend, info)]
        return interpCurves

    def getXLimits(self, values, overlap=True):
        if overlap:
            xmin0, xmax0 = -numpy.inf, numpy.inf
        else:
            xmin0, xmax0 = numpy.inf, -numpy.inf
        for x in values:
            xmin = x.min()
            xmax = x.max()
            if overlap:
                if xmin > xmin0:
                    xmin0 = xmin
                if xmax < xmax0:
                    xmax0 = xmax
            else:
                if xmin < xmin0:
                    xmin0 = xmin
                if xmax > xmax0:
                    xmax0 = xmax
        if DEBUG:
            print('-- getXLimits : overlap =',overlap,
                  ', xmin = ',xmin0,', xmax =',xmax0)
        return xmin0, xmax0

    def fftShift(self, shift, x, y):
        print 'fftShift -- shift =', shift
        freq = numpy.fft.fftfreq(len(x), d=(x[1]-x[0]))
        yfft = numpy.fft.fft(y)
#        y = numpy.fft.ifft(
#             numpy.exp(-2.0*numpy.pi*numpy.sqrt(numpy.complex(-1))*\
#                numpy.fft.fftfreq(len(x), d=x[1]-x[0])*shift)*numpy.fft.fft(y))
        yShifted = numpy.fft.ifft(
             numpy.exp(-2.0*numpy.pi*numpy.sqrt(numpy.complex(-1))*\
                numpy.fft.fftfreq(len(x), d=x[1]-x[0])*shift)*numpy.fft.fft(y))
#        yShifted = numpy.fft.ifft(
#                    numpy.exp(-2.0*numpy.pi*numpy.sqrt(numpy.complex(-1))*\
#                     freq*shift) * yfft)
        return x, yShifted.real

    def xShift(self, shift, x, y):
        return x+shift, y

    def applyShifts(self):
        if len(self.shiftDict) == 0:
            msg = qt.QMessageBox(None)
            msg.setWindowTitle('Alignment Error')
            msg.setText('No shift data present.')
            msg.setStandardButtons(qt.QMessageBox.Ok)
            return
#        if self.shiftMethod == 'FFT':
        if self.shiftMethod == self.fftShift:
            curves = self.interpolate()
        else:
            curves = self.getAllCurves()
        for idx, (x,y,legend,info) in enumerate(curves):
            replace=False
            selectionlegend = info.get('selectionlegend',legend)
            shift = self.shiftDict.get(selectionlegend,None)
            if shift is None:
                if DEBUG:
                    print('Curve \'%s\' not found in shiftDict\n%s'%(selectionlegend,str(self.shiftDict)))
                continue
            if shift == float('NaN'):
                if DEBUG:
                    print('Curve \'%s\' has NaN shift'%selectionlegend)
                continue
            xShifted, yShifted = self.shiftMethod(shift, x, y)
            if idx == 0:
                replace = True
            self.addCurve(xShifted, yShifted, (legend + ' SHIFT'), info, replace, True)
#            self.addCurve(x, y, legend=(legend + ' SHIFT'), info=info, replace=replace, replot=replot)
        
MENU_TEXT = "TAKE ME Alignment Plugin"
def getPlugin1DInstance(plotWindow, **kw):
    ob = AlignmentScanPlugin(plotWindow)
    return ob

if __name__ == "__main__":
    from PyMca import PyMcaQt as qt
#    app = qt.QApplication([])
    from PyMca.Plot1DQwt import Plot1DQwt as Plot1D
#    i = numpy.arange(100.)
#    y1 = 10.0 + 500.0 * numpy.exp(-0.01*(i-50)**2)
#    y2 = 10.0 + 500.0 * numpy.exp(-((i-55)/5.)**2)
#    x = numpy.arange(250, 750, 2, dtype=float)
#    y1 = 1.0 + 50.0 * numpy.exp(-0.001*(x-500)**2) + 2.*numpy.random.random(250.)
#    y2 = 1.0 + 20.5 * numpy.exp(-0.005*(x-600)**2) + 2.*numpy.random.random(250.)
#    plot = Plot1D()
##    plot.addCurve(i, y1, "y1")
#    plot.addCurve(x, y1, "y1")
##    plot.addCurve(i, y2, "y2")
#    plot.addCurve(x, y2, "y2")
#    plugin = getPlugin1DInstance(plot)
#    for method in plugin.getMethods():
#        print(method, ":", plugin.getMethodToolTip(method))
#    plugin.applyMethod(plugin.getMethods()[0])
#    curves = plugin.getAllCurves()
#    #for curve in curves:
#    #    print(curve[2])
#    print("LIMITS = ", plugin.getGraphYLimits())
    app = qt.QApplication([])
    plot = AlignmentWidget(None, {'foo': 3.14, 'bar': 2.74})
    plot.show()
    app.exec_()
#    print plugin.shiftDict
