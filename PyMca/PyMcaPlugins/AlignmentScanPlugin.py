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
from PyMca import specfilewrapper as SFW
from PyMca import SpecfitFunctions as SF
from PyMca import SNIPModule as snip
from PyMca.Gefit import LeastSquaresFit as LSF
from PyMca.SpecfitFuns import gauss

try:
    from PyMca import Plugin1DBase
except ImportError:
    print("WARNING:AlignmentScanPlugin import from somewhere else")
    from . import Plugin1DBase

from PyMca import SpecfitFuns

DEBUG = 0
class AlignmentWidget(qt.QDialog):
    def __init__(self, parent, ddict, plugin):
        qt.QDialog.__init__(self, parent)
        self.setWindowTitle('Alignment Window')
        
        nCols = 2
        nRows = len(ddict)
        self.plugin = plugin
        
        # Buttons
        buttonSave = qt.QPushButton('Save')
        buttonSave.setToolTip('Save shifts to file')
        buttonLoad = qt.QPushButton('Load')
        buttonLoad.setToolTip('Load shifts from file')
        buttonStore = qt.QPushButton('Store')
        buttonStore.setToolTip('Store shifts in memory.\n')
        buttonApply = qt.QPushButton('Apply')
        buttonApply.setToolTip('Apply shift to curves present'
                              +' in the plot window')
        buttonCancel = qt.QPushButton('Cancel')
        buttonCalc = qt.QPushButton('Calculate')

        # Table
        self.shiftTab = qt.QTableWidget(nRows, nCols)
        self.shiftTab.verticalHeader().hide()
        self.shiftTab.horizontalHeader().setStretchLastSection(True)
        self.shiftTab.setHorizontalHeaderLabels(['Legend','Shift'])
        
        # Shift Method selector
        self.shiftMethodComboBox = qt.QComboBox()
        self.shiftMethodComboBox.addItems(
            ['Shift x-range',
            'Inverse FFT shift'])
        shiftMethodToolTip =\
            ('Select the method that shifts the spectra\n\n'
            +'Shift x-range:\n'
            +'     Directly applies the shift to the data\'s\n'
            +'     x-range\n'
            +'Inverse FFT shift:\n'
            +'     Shifts the spectra by multiplying a\n'
            +'     phase factor to their Fourier transform. The result is\n'
            +'     transformed back to real space. Recommended for data with\n'
            +'     resp. regions with constant background.')
        self.shiftMethodComboBox.setToolTip(shiftMethodToolTip)
             
        # Alignment Method selector
        self.alignmentMethodComboBox = qt.QComboBox()
        self.alignmentMethodComboBox.addItems(
            ['FFT',
             'MAX',
             'FIT',
             'FIT DRV'])
        alignmentMethodToolTip =\
            ('Select the method used to calculate the shift is calculated.\n\n'
            +'FFT:\n'
            +'     Calculates the correlation between two curves using its\n'
            +'     Fourier transform. The shift is proportional to the distance of\n'
            +'     the correlation function\'s maxima.\n'
            +'MAX:\n'
            +'     Determines the shift as the distance between the maxima of\n'
            +'     two peaks\n'
            +'FIT:\n'
            +'     Guesses the most prominent feature in a spectrum and tries\n'
            +'     to fit it with a Gaussian peak. Before the fit is perform, the\n'
            +'     background is substracted. The shift is given by the difference\n'
            +'     of the center of mass between two peaks.\n'
            +'FIT DRV:\n'
            +'     Like FIT, but the fit is performed on the derivate of the\n'
            +'     spectrum. Recommended procedure for XAFS data.')
        self.alignmentMethodComboBox.setToolTip(alignmentMethodToolTip)
        
        # Fill table with data
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
        mainLayout.setContentsMargins(1,1,1,1)
        self.setLayout(mainLayout)
        
        # Connects
        self.shiftTab.cellChanged.connect(self.validateInput)
        buttonApply.clicked.connect(self.accept)
        buttonCancel.clicked.connect(self.reject)
        buttonStore.clicked.connect(self.store)
        buttonSave.clicked.connect(self.saveDict)
        buttonLoad.clicked.connect(self.loadDict)

        # ..to Plugin instance
        buttonCalc.clicked[()].connect(self.triggerCalculateShift)
        self.alignmentMethodComboBox.activated['QString'].\
                            connect(self.triggerCalculateShift)

    def triggerCalculateShift(self, methodName=None):
        if methodName != None:
            self.plugin.setAlignmentMethod(methodName)
        ddict = self.plugin.calculateShifts()
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
        if len(filename) == 0:
            return
        inDict = ConfigDict.ConfigDict()
        try:
            inDict.read(filename)
        except IOError:
            msg = qt.QMessageBox()
            msg.setTitle('FFTAlignment Load Error')
            msg.setText('Unable to read shifts form file \'%s\''%filename)
            msg.exec_()
            return
        if 'Shifts' not in inDict.keys():
            # Only if the shift file consists exclusively of ShiftList
            orderedLegends = [legend for legend in self.plugin.getOrder()]
            try:
                shiftList = inDict['Ordered Shifts']['ShiftList']
            except KeyError:
                msg = qt.QMessageBox()
                msg.setWindowTitle('FFTAlignment Load Error')
                msg.setText('No shift information found in file \'%s\''%filename)
                msg.exec_()
            ddict = dict(zip(orderedLegends, shiftList))
        else:
            ddict = (inDict['Shifts'])
        self.setDict(ddict)

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
            return False
        if len(filename) == 0:
            return False
        if not str(filename).endswith('.shift'):
            filename += '.shift'
        if DEBUG:
            print('saveOptions -- Filename: "%s"' % filename)
        outDict = ConfigDict.ConfigDict()
        ddict = self.getDict()
        outDict['Shifts'] = ddict
        outDict['ShiftList'] = [ddict[legend] for legend in self.plugin.getOrder()]
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
        return ddict
    
    def setDict(self, ddict):
        # Conserve order in which scans are added to plot window
        #dkeys = sorted(ddict.keys())
        tmp = self.plugin.getOrder()
        dkeys = [item for item in tmp if item in ddict.keys()]
        dvals = ['%f'%ddict[k] for k in dkeys]
        
        self.shiftTab.clear()
        self.shiftTab.setHorizontalHeaderLabels(['Legend','Shift'])
        self.shiftTab.setColumnCount(2)
        self.shiftTab.setRowCount(len(dkeys))
        if len(ddict) == 0:
            return

        for j, dlist in enumerate([dkeys, dvals]):
            for i in range(len(dlist)):
                if j == 0:
                    elem = qt.QTableWidgetItem(dlist[i])
                    elem.setFlags(qt.Qt.ItemIsSelectable)
                    #elem.setFlags(qt.Qt.ItemIsEnabled)
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

    def validateInput(self, row, col):
        if col == 0:
            return
        elif col == 1:
            item  = self.shiftTab.item(row, 1)
            try:
                floatValue = float(item.text())
                item.setText('%.6g'%floatValue)
            except:
                floatValue = float('NaN')
                item.setText(str(floatValue))
        
class AlignmentScanPlugin(Plugin1DBase.Plugin1DBase):
    def __init__(self, plotWindow, **kw):
        Plugin1DBase.Plugin1DBase.__init__(self, plotWindow, **kw)
        self.__randomization = True
        self.__methodKeys = []
        self.methodDict = {}

        function = self.calculateAndApplyShifts
        method = "Perform FFT Alignment"
        text = "Performs FFT based alignment and inverse FFT based shift"
        info = text
        icon = None
        self.methodDict[method] = [function,
                                   info,
                                   icon]
        self.__methodKeys.append(method)
    
        function = self.showShifts
        method = "Show Alignment Window"
        text  = "Displays the calculated shifts and"
        text += "allows to fine tune the plugin"
        info = text
        icon = None
        self.methodDict[method] = [function,
                                   info,
                                   icon]
        self.__methodKeys.append(method)
        
        self.alignmentMethod = self.calculateShiftsFFT
        self.shiftMethod     = self.fftShift
        self.shiftDict       = {}
        
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


    def calculateAndApplyShifts(self):
        # Assure that FFT alignment & shift methods are set
        self.alignmentMethod = self.calculateShiftsFFT
        self.shiftMethod     = self.fftShift
        self.calculateShifts()
        self.applyShifts()
        # Reset shift Dictionary
        self.shiftDict = {}

    def calculateShifts(self):
        '''
        Generic alignment method, executes the method
        that is set under self.alignmentMethod.
        
        Choices are:
        - calculateShiftsFit
        - calculateShiftsFFT
        - calculateShiftsMax
        '''
        self.shiftDict = self.alignmentMethod()
        return self.shiftDict

    def getOrder(self):
        ret = [legend for (x,y,legend,info) in self._plotWindow.getAllCurves()]
        return ret

    # BEGIN Alignment Methods
    def calculateShiftsFitDerivative(self):
        return self.calculateShiftsFit(derivative=True)
    
    def calculateShiftsFit(self, derivative=False, thr=30):
        retDict = {}
        
        curves = self.getAllCurves()
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
        if DEBUG:
            print('calculateShiftsFit -- xmin = %.3f, xmax = %.3f'%(xmin, xmax))

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
        
        if derivative:
            # Take first derivative
            y0 = numpy.diff(y0)/numpy.diff(x0)
            x0 = .5 * (x0[1:] + x0[:-1])
        
        peak0 = self.findPeaks(x0, y0, .80, derivative)
        if peak0:
            xp0, yp0, fwhm0, fitrange0 = peak0
        else:
            raise ValueError("No peak identified in '%s'"%activeCurve[2])
        fitp0, chisq0, sigma0 = LSF(gauss,
                                    numpy.asarray([yp0, xp0, fwhm0]), 
                                    xdata=x0[fitrange0], 
                                    ydata=y0[fitrange0])
        if DEBUG:
            if derivative:
                print('calculateShiftsFit -- Results (Leg, PeakPos, Shift):')
            else:
                print('calculateShiftsFitDerivative -- Results (Leg, PeakPos, Shift):')
        for x,y,legend,info in curves:
            idx = numpy.nonzero((xmin <= x) & (x <= xmax))[0]
            x = numpy.take(x, idx)
            y = numpy.take(y, idx)
            
            if derivative:
                # Take first derivative
                y = numpy.diff(y)/numpy.diff(x)
                x = .5 * (x[1:] + x[:-1])
            
            peak = self.findPeaks(x, y, .80, derivative)
            if peak:
                xp, yp, fwhm, fitrange = peak
            else:
                raise ValueError("No peak identified in '%s'"%activeCurve[2])
            try:
                fitp, chisq, sigma = LSF(gauss,
                                         numpy.asarray([yp, xp, fwhm]), 
                                         xdata=x[fitrange], 
                                         ydata=y[fitrange])
                # Shift is difference in peak's x position
                shift = fitp0[1] - fitp[1]
            except numpy.linalg.linalg.LinAlgError:
                msg = qt.QMessageBox(None)
                msg.setWindowTitle('Alignment Error')
                msg.setText('Singular matrix encountered during least squares fit.')
                msg.setStandardButtons(qt.QMessageBox.Ok)
                msg.exec_()
                shift = float('NaN')
            key = legend
            retDict[key] = shift
            if DEBUG:
                  print( '\t%s\t%.3f\t%.3f'%(legend, fitp[1], shift))
        return retDict

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
        if DEBUG:
            print('calculateShiftsMax -- Results:')
            print('\targmax(y) shift')
        for x,y,legend,info in curves:
            idx = numpy.nonzero((xmin <= x) & (x <= xmax))[0]
            x = numpy.take(x, idx)
            y = numpy.take(y, idx)
            
            shifty = numpy.argmax(y)
            shift = x0[shift0] - x[shifty]
            key = legend
            retDict[key] = shift
            if DEBUG:
                print('\t%d %.3f'%(x[shifty],shift))
        return retDict

    def calculateShiftsFFT(self, portion=.95):
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
        if DEBUG:
            print('calculateShiftsFFT -- xmin = %.3f, xmax = %.3f'%(xmin, xmax))

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
        y0 = self.normalize(y0)
        y0 = numpy.take(y0, idx)

        fft0 = numpy.fft.fft(y0)
        if DEBUG:
            print('calculateShiftsFFT -- results (Legend len(idx) shift):')
        for x,y,legend,info in curves:
            idx = numpy.nonzero((x >= xmin) & (x <= xmax))[0]
            x = numpy.take(x, idx)
            y = numpy.take(y, idx)
            ffty = numpy.fft.fft(y)
            shiftTmp = numpy.fft.ifft(fft0 * ffty.conjugate()).real
            shiftPhase = numpy.zeros(shiftTmp.shape, dtype=shiftTmp.dtype)
            m = shiftTmp.size//2
            shiftPhase[m:] = shiftTmp[:-m]
            shiftPhase[:m] = shiftTmp[-m:]
            # Normalize shiftPhase to standardize thresholding
            shiftPhase = self.normalize(shiftPhase)
            
            # Thresholding
            xShiftMax = shiftPhase.argmax()
            left, right = xShiftMax, xShiftMax
            threshold = portion * shiftPhase.max()
            while (shiftPhase[left] > threshold)&\
                  (shiftPhase[right] > threshold):
                left  -= 1
                right += 1
            idx = numpy.arange(left, right+1, 1, dtype=int)
            # The shift is determined by center-of-mass around shiftMax
            shiftTmp = (shiftPhase[idx] * idx/shiftPhase[idx].sum()).sum()
            shift = (shiftTmp - m) * (x[1] - x[0])

            key = legend
            retDict[key] = shift
            if DEBUG:
                print('\t%s\t%d\t%f'%(legend,len(idx),shift))
        return retDict
    # END Alignment Methods

    def applyShifts(self):
        '''
        Generic shift method. The method shifts curves
        according to the shift stored in self.shiftDict
        and executes the method stored in self.shiftMethod.
        
        Curves are sorted with respect to their legend,
        the values of self.shiftDict are sorted with
        respect to their key.
        '''
        if len(self.shiftDict) == 0:
            msg = qt.QMessageBox(None)
            msg.setWindowTitle('Alignment Error')
            msg.setText('No shift data present.')
            msg.setStandardButtons(qt.QMessageBox.Ok)
            msg.exec_()
            return False
            
        if self.shiftMethod == self.fftShift:
            curves = self.interpolate()
        else:
            curves = self.getAllCurves()
        
        # Sort shiftDict, curves
        shiftDictKeys = sorted(self.shiftDict.keys()) # Contains legends of shifts
        #curves = sorted(curves, key=lambda curve: curve[2])
        sortedLegends = sorted([legend for (x,y,legend,info) in curves])
        #sortedCurves  = sorted(curves, key=lambda curve: curve[2]) # Contains curves from plot window
        
        if len(shiftDictKeys) != len(curves):
            msg = qt.QMessageBox(None)
            msg.setWindowTitle('Alignment Error')
            msg.setText(
                '''Number of shifts does not match the number of curves.
                Do you want to continue anyway?''')
            msg.setStandardButtons(qt.QMessageBox.Ok)
            msg.setStandardButtons(qt.QMessageBox.Ok | qt.QMessageBox.Cancel)
            msg.setDefaultButton(qt.QMessageBox.Ok)
            
            if msg.exec_() != qt.QMessageBox.Ok:
                return False
        
        # Ensure correct assign between shift values and
        # curves that are currently present in plot window
        assignmentDict = dict(zip(sortedLegends, shiftDictKeys))

        if DEBUG:
            print('applyShifts -- Shifting ...')
        for idx, (x,y,legend,info) in enumerate(curves):
            selectionlegend = info.get('selectionlegend',legend)
            shiftLegend = assignmentDict[legend]
            shift = self.shiftDict[shiftLegend]
            if shift is None:
                if DEBUG:
                    print('\tCurve \'%s\' not found in shiftDict\n%s'\
                          %(legend,str(self.shiftDict)))
                continue
            if shift == float('NaN'):
                if DEBUG:
                    print('\tCurve \'%s\' has NaN shift'%selectionlegend)
                continue
            
            # Limit shift to zoomed in area
            xmin, xmax = self.getGraphXLimits()
            mask = numpy.nonzero((xmin<=x) & (x<=xmax))[0]
            # Execute method stored in self.shiftMethod
            xShifted, yShifted = self.shiftMethod(shift, x[mask], y[mask])
            
            if DEBUG:
                print('\t leg = %s, key = %s,shift = %f'%(legend,key,shift))
            
            if idx == 0:
                replace, replot = True, False
            elif idx == (len(curves)-1):
                replace, replot = False, True
            else:
                replace, replot = False, False
            # Check if scan number is adopted by new curve    
            #print "applyShifts -- %s -> '%s'"%(legend,info.get('Key'))
            self.addCurve(xShifted, yShifted, 
                          (selectionlegend + ' SHIFT'),
                          info,
                          replace, replot)
        return True


    # BEGIN Shift Methods
    def fftShift(self, shift, x, y):
        yShifted = numpy.fft.ifft(
             numpy.exp(-2.0*numpy.pi*numpy.sqrt(numpy.complex(-1))*\
                numpy.fft.fftfreq(len(x), d=x[1]-x[0])*shift)*numpy.fft.fft(y))
        return x, yShifted.real

    def xShift(self, shift, x, y):
        return x+shift, y
    # END Shift Methods

    def showShifts(self):
        '''
        Creates an instance of Alignment Widget that
        allows to
        
        - Calculate, display  & save/store shifts
        - Load existing shift data
        - Select different alignment and shift methods
        '''
        if 0:
        #if len(self.shiftDict) == 0:
            msg = qt.QMessageBox(None)
            msg.setWindowTitle('Alignment Error')
            msg.setText('No shift data present.\nDo you want to calculate shifts now?')
            msg.setStandardButtons(qt.QMessageBox.Ok | qt.QMessageBox.Close)
            msg.setDefaultButton(qt.QMessageBox.Ok)
            ret = msg.exec_()
            if ret == qt.QMessageBox.Ok:
                self.calculateShifts()
            else:
                return
        # Empty shift table in the beginning
        widget = AlignmentWidget(None, self.shiftDict, self)
        ret = widget.exec_()
        if ret == 1:
            # Result code Apply
            self.shiftDict = widget.getDict()
            self.setShiftMethod(widget.getShiftMethodName())
            self.applyShifts()
            self.shiftDict = {}
        elif ret == 2:    
            # Result code Store
            self.shiftDict = widget.getDict()
            self.setShiftMethod(widget.getShiftMethodName())
        else:
            # Dialog is canceled
            self.shiftDict = {}
        #widget.destroy() # Widget should be destroyed after finishing method
        return

    # BEGIN Helper Methods
    def setShiftMethod(self, methodName):
        '''
        Method receives methodName from AlignmentWidget
        instance and assigns the according shift method. 
        '''
        if DEBUG:
            print('setShiftMethod -- %s'%methodName)
        methodName = str(methodName)
        if methodName == 'Inverse FFT shift':
            self.shiftMethod = self.fftShift
        elif methodName == 'Shift x-range':
            self.shiftMethod = self.xShift
        else:
            # Unknown method name, use fftShift as default
            self.shiftMethod = self.fftShift

    def setAlignmentMethod(self, methodName):
        '''
        Method receives methodName from AlignmentWidget
        instance and assigns the according alignment method. 
        '''
        if DEBUG:
            print('setAlignmentMethod -- %s'%methodName)
        methodName = str(methodName)
        if methodName == 'FFT':
            self.alignmentMethod = self.calculateShiftsFFT
        elif methodName == 'MAX':
            self.alignmentMethod = self.calculateShiftsMax
        elif methodName == 'FIT':
            self.alignmentMethod = self.calculateShiftsFit
        elif methodName == 'FIT DRV':
            self.alignmentMethod = self.calculateShiftsFitDerivative
        else:
            # Unknown method name, use fftShift as default
            self.alignmentMethod = self.calculateShiftsFFT

    def getAllCurves(self, just_legend=False):
        '''
        Ensures that the x-range of the curves
        is strictly monotonically increasing.
        Conserves curves legend and info dictionary.
        '''
        curves = Plugin1DBase.Plugin1DBase.getAllCurves(self)
        if just_legend:
            return curves
    
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
            print "getAllCurves -- %s -> '%s'"%(legend,info.get('Key'))
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
        Use this method instead of self.getAllCurves() when
        performin FFT related tasks.
        
        Returns
        -------
        interpCurves : ndarray
            Array containing the interpolated curves shown
            in the plot window. 
            Format: [(x0, y0, legend0, info0), ...]
        '''
        curves = self.getAllCurves()
        if len(curves) < 1:
            raise ValueError("At least 1 curve needed")
            if DEBUG:
                print('interpolate -- no curves present')
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
        '''
        Input
        -----
        overlap : bool
            True  -> returns minimal and maximal x-values
                     that are that are still lie within the 
                     x-ranges of all curves in plot window
            False -> returns minimal and maximal x-values of
                     all curves in plot window
                     
        Returns
        -------
        xmin0, xmax0 : float
        '''
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
            print('getXLimits -- overlap = %s, xmin = %.3f, xmax =%.3f'\
                  %(overlap,xmin0,xmax0))
        return xmin0, xmax0

    def normalize(self, y):
        ymax, ymin = y.max(), y.min()
        return (y-ymin)/(ymax-ymin)

    def findPeaks(self, x, y, thr, derivative):
        '''
        Input
        -----
        x,y : ndarrays
            Arrays contain curve intformation
        thr : float
            Threshold in percent of normalized maximum
        derivative : bool
            The derivative of a curve is being fitted

        Finds most prominent feature contained in y
        and tries to estimate starting parameters for a
        Gaussian least squares fit (LSF). Recommends values
        used to fit the Gaussian.
        
        Return
        ------
        xpeak, ypeak, fwhm : float
            Estimated values for x-position, amplitude
            and width of the Gaussian
        fwhmIdx : ndarray
            Indices determine the range on which the LSF
            is performed
        '''
        # Use SNIP algorithm for background substraction &
        # seek method for peak detection
        sffuns = SF.SpecfitFunctions()
        if derivative:
            # Avoid BG substraction & normalization if
            # fitting the derivate of a curve
            ybg = y
            ynorm = y/(abs(y.max())+abs(y.min()))
        else:
            ybg = y-snip.getSnip1DBackground(y, len(y)//thr) # USER INPUT!!!
            # Normalize background substracted data to
            # standardize the yscaling of seek method
            #ynorm = (ybg - ybg.min())/(ybg.max()-ybg.min())
            ynorm = self.normalize(ybg)

        # Replace by max()?
        try:
            # Calculate array woth all peak indices
            peakIdx = numpy.asarray(sffuns.seek(ybg, yscaling=1000.), dtype=int)
            # Extract highest peak
            sortIdx = y[peakIdx].argsort()[-1]
        except IndexError:
            if DEBUG:
                print('No peaks found..')
            return None
        except SystemError:
            if DEBUG:
                print('Peak search failed. Continue with y maximum')
            peakIdx = [ybg.argmax()]
            sortIdx = 0
        xpeak = float(x[peakIdx][sortIdx])
        ypeak = float(y[peakIdx][sortIdx])
        ypeak_norm = float(ynorm[peakIdx][sortIdx])
        ypeak_bg   = float(ybg[peakIdx][sortIdx])
        
        # Estimate FWHM
        fwhmIdx = numpy.nonzero(ynorm >= thr*ypeak_norm)[0]
        #fwhmIdx = numpy.nonzero(ybg >= thr*ypeak_bg)[0]
        # Underestimates FWHM
        x0, x1 = x[fwhmIdx].min(), x[fwhmIdx].max()
        fwhm = x1 - x0
        
        return xpeak, ypeak, fwhm, fwhmIdx
    # END Helper Methods
    
MENU_TEXT = "Alignment Plugin"
def getPlugin1DInstance(plotWindow, **kw):
    ob = AlignmentScanPlugin(plotWindow)
    return ob

if __name__ == "__main__":
    from PyMca import PyMcaQt as qt
    app = qt.QApplication([])
    from PyMca.Plot1DQwt import Plot1DQwt as Plot1D

    x = numpy.arange(250, 750, 2, dtype=float)
    y1 = 1.0 + 50.0 * numpy.exp(-0.001*(x-500)**2) + 2.*numpy.random.random(250.)
    y2 = 1.0 + 20.5 * numpy.exp(-0.005*(x-600)**2) + 2.*numpy.random.random(250.)

    plot = Plot1D()
    plot.addCurve(x, y1, "y1", {'selectionlegend': 'y1'})
    plot.addCurve(x, y2, "y2", {'selectionlegend': 'y2'})

    plugin = getPlugin1DInstance(plot)
    for method in plugin.getMethods():
        print(method, ":", plugin.getMethodToolTip(method))
    plugin.applyMethod(plugin.getMethods()[0])
