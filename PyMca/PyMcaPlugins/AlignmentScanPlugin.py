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

try:
    from PyMca import Plugin1DBase
except ImportError:
    print("WARNING:AlignmentScanPlugin import from somewhere else")
    from . import Plugin1DBase

from PyMca import SpecfitFuns

#class TableItem(qt.QTableWidgetItem):
#    def __init__(self, parent, content):
        

class AlignmentWidget(qt.QDialog):
    def __init__(self, parent, ddict):
        qt.QDialog.__init__(self, parent)
        self.setWindowTitle('FFT Shifts')
        
        nCols = 2
        nRows = len(ddict)
        
        # Buttons
        buttonApply  = qt.QPushButton('Apply')
        buttonCancel = qt.QPushButton('Cancel')

        # Table
        shiftTab = qt.QTableWidget(nRows, nCols)
        shiftTab.verticalHeader().hide()
        shiftTab.setHorizontalHeaderLabels(['Legend','Shift'])
        shiftTab.horizontalHeader().setStretchLastSection(True)
        
        # Fill table with data
        dkeys = sorted(ddict.keys())
        dvals = ['%.4f'%ddict[k] for k in dkeys]
        for j, dlist in enumerate([dkeys, dvals]):
            for i in range(len(dlist)):
                if j == 0:
                    elem = qt.QTableWidgetItem(dlist[i])
                elif j == 1:
                    elem = qt.QTableWidgetItem(str(dlist[i]))
                    elem.setTextAlignment(qt.Qt.AlignRight)
                    elem.setTextAlignment(qt.Qt.AlignRight + qt.Qt.AlignVCenter)
                else:
                    elem = qt.QTableWidgetItem('')
                elem.setFlags(qt.Qt.ItemIsEnabled)
                shiftTab.setItem(i,j, elem)
        shiftTab.resizeColumnToContents(0)    
        shiftTab.resizeRowsToContents()
        
        buttonLayout = qt.QHBoxLayout()
        buttonLayout.addWidget(buttonApply)
        buttonLayout.addWidget(qt.HorizontalSpacer())
        buttonLayout.addWidget(buttonCancel)
        
        mainLayout = qt.QVBoxLayout()
        mainLayout.addWidget(shiftTab)
        mainLayout.addLayout(buttonLayout)
        self.setLayout(mainLayout)
        
        # Connects
        buttonApply.clicked.connect(self.accept)
        buttonCancel.clicked.connect(self.close)

DEBUG = False
class AlignmentScanPlugin(Plugin1DBase.Plugin1DBase):
    def __init__(self, plotWindow, **kw):
        Plugin1DBase.Plugin1DBase.__init__(self, plotWindow, **kw)
        self.__randomization = True
        self.__methodKeys = []
        self.methodDict = {}
        text = "Perform FFT based alignment"
        info = text
        icon = None
        function = self.fftAlignment
        method = "Perform FFT Alignment"
        self.methodDict[method] = [function,
                                   info,
                                   icon]
        self.__methodKeys.append(method)
        text = "Apply calculated shifts by FFT based alignment to curves"
        info = text
        icon = None
        function = self.applyAlignment
        method = "Apply FFT Alignment"
        self.methodDict[method] = [function,
                                   info,
                                   icon]
        self.__methodKeys.append(method)
        text = "Displays shifts calculated by FFT based alignment to curves"
        info = text
        icon = None
        function = self.showAlignment
        method = "Calculate shifts"
        self.methodDict[method] = [function,
                                   info,
                                   icon]
        self.__methodKeys.append(method)
        
        self.shiftDict = {}
        
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

    def resetAlignment(self):
        self.shiftDict = {}

    def showAlignment(self):
        if len(self.shiftDict) == 0:
            self.fftAlignment(plot=False)
        ret = AlignmentWidget(None, self.shiftDict)
        if ret.exec_():
            self.applyAlignment()

    def applyAlignment(self):
        if len(self.shiftDict) == 0:
            msg = qt.QMessageBox()
            msg.setWindowTitle('Alignment Error')
            msg.setText('No shifts calculated yet!')
            msg.exec_()
            return
        curves = self.getAllCurves()
        active = self.getActiveCurve()
        if not active:
            msg = qt.QMessageBox()
            msg.setWindowTitle('Alignment Error')
            msg.setText('No active curve selected!')
            msg.exec_()
            return
        else:
            (xa,ya,legendActive,infoActive) = active
        count = 0
        for (x, y, legend, info) in curves:
            if not count:
                replace = True
            else:
                replace = False
            selLegend = info['selectionlegend']
            if selLegend == infoActive.get('selectionlegend'):
                self.addCurve(x,
                          y,
                          legend=legend,
                          info=None,
                          replot=True,
                          replace=replace)
            else:
                shift = self.shiftDict.get(selLegend,None)
                try:
                    y = numpy.fft.ifft(numpy.exp(-2.0*numpy.pi*numpy.sqrt(numpy.complex(-1))*\
                                        numpy.fft.fftfreq(len(x), d=x[1]-x[0])*shift)*numpy.fft.fft(y))
                except ValueError:
                    print 'Shapes:'
                    print '\td = ',(x[1]-x[0])
                    print '\tlen(x) =',len(x)
                    print '\tlen(y) =',len(y)
                    return
                y = y.real
                if not count:
                    replace = True
                self.addCurve(x,
                              y,
                              legend=legend,
                              info=None,
                              replot=True,
                              replace=replace)
            count += 1
        return

    def fftAlignment(self, plot=True):
        self.shiftDict = {}
        
        curves = self.getAllCurves()
        nCurves = len(curves)
        if nCurves < 2:
            raise ValueError("At least 2 curves needed")
            return

        processedCurves = []
        ### CURVE PREPROCESSING ###
        ### Sort
        ### Ravel
        ### Increasing
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
        
        # Determine larges overlap between curves
        xmin, xmax = self.getGraphXLimits()
        for x,y,legend,info in processedCurves:
            try:
                if xmin<x.min():
                    xmin = x.min()
                if x.max()<xmax:
                    xmax = x.max()
            except ValueError as e:
                print e
                print len(x)

        # get active curve            
#        activeCurve = self.getActiveCurve()
#        if activeCurve is None:
##            activeCurve = curves[0]
#            activeCurve = processedCurves[0]
        activeCurve = self.getActiveCurve()
        if activeCurve is None:
#            activeCurve = curves[0]
            activeCurve = processedCurves[0]
        else:
            (xa,ya,activeLegend,activeInfo) = activeCurve
            legends = [legend for (x,y,legend,info) in processedCurves]
            idx = legends.index(activeLegend)
            activeCurve = processedCurves[idx]
            
        # apply between graph limits
        x0 = activeCurve[0][:]
        y0 = activeCurve[1][:]

        #sort the values
#        idx = numpy.argsort(x0, kind='mergesort')
#        x0 = numpy.take(x0, idx)
#        y0 = numpy.take(y0, idx)

        #remove duplicates
#        x0 = x0.ravel()
#        idx = numpy.nonzero((x0[1:] > x0[:-1]))[0]
#        x0 = numpy.take(x0, idx)
#        y0 = numpy.take(y0, idx)
        


        print 'xmin =',xmin,', xmax =',xmax
#        idx = numpy.nonzero((x0 >= xmin) & (x0 <= xmax))[0]
        idx = numpy.nonzero((xmin < x0) & (x0 < xmax))[0]
        x0 = numpy.take(x0, idx)
        y0 = numpy.take(y0, idx)

        #make sure values are regularly spaced
        xi = numpy.linspace(x0[0], x0[-1], len(idx)).reshape(-1, 1)
        print 'Equidistant'
        print '\tlen(x) =',len(xi)
        yi = SpecfitFuns.interpol([x0], y0, xi, y0.min())
        x0 = xi
        y0 = yi
#        print 'After Interp1'
#        print '\tlen(x) =',len(x0)
#        print '\tlen(y) =',len(y0)
        y0.shape = -1
        fft0 = numpy.fft.fft(y0)
        y0.shape = -1, 1
        x0.shape = -1, 1
        nChannels = x0.shape[0]

        # built a couple of temporary array of spectra for handy access
        ret = {}
        tmpArray = numpy.zeros((nChannels, nCurves), numpy.float)
        fftList = []
        shiftList = []
        curveList = []
        i = 0
        for idx in range(nCurves):
#            x, y, legend, info = curves[idx][0:4]
            x, y, legend, info = processedCurves[idx][0:4]
            print 'From Curves'
            print '\tlen(x) =',len(x)
            print '\tlen(y) =',len(y)
            #sort the values
            x = x[:]
#            idx = numpy.argsort(x, kind='mergesort')
#            x = numpy.take(x, idx)
#            y = numpy.take(y, idx)

            #take the portion of x between limits
            idx = numpy.nonzero((x>=xmin) & (x<=xmax))[0]
            if not len(idx):
                # no overlap
                if DEBUG:
                    print 'fftAlignment -- no overlap between curves'
                continue
#            x = numpy.take(x, idx)
#            y = numpy.take(y, idx)

            #remove duplicates
#            x = x.ravel()
#            idx = numpy.nonzero((x[1:] > x[:-1]))[0]
#            x = numpy.take(x, idx)
#            y = numpy.take(y, idx)
#            x.shape = -1, 1
#            print 'Before Interp'
#            print '\tlen(x) =',len(x)
#            print '\tlen(y) =',len(y)
            if numpy.allclose(x, x0):
                # no need for interpolation
                pass
            else:
                # we have to interpolate
                x.shape = -1
                y.shape = -1
                xi = x0[:]
                y = SpecfitFuns.interpol([x], y, xi, y0.min())
                x = xi
            y.shape = -1
            tmpArray[:, i] = y
            i += 1
            print 'After Interp'
            print '\tlen(x) =',len(x)
            print '\tlen(y) =',len(y)

            # now calculate the shift
            ffty = numpy.fft.fft(y)
            fftList.append(ffty)
            if 0:
                self.addCurve(x,
                      y,
                      legend="NEW Y%d" % i,
                      info=None,
                      replot=True,
                      replace=False)
            elif numpy.allclose(fft0, ffty):
                shiftList.append(0.0)
                shift = 0.0
            else:
                shift = numpy.fft.ifft(fft0 * ffty.conjugate()).real
                shift2 = numpy.zeros(shift.shape, dtype=shift.dtype)
                m = shift2.size//2
                shift2[m:] = shift[:-m]
                shift2[:m] = shift[-m:]
                if 0:
                    self.addCurve(numpy.arange(len(shift2)),
                          shift2,
                          legend="SHIFT",
                          info=None,
                          replot=True,
                          replace=False)
                threshold = 0.50*shift2.max()
                #threshold = shift2.mean()
                idx = numpy.nonzero(shift2 > threshold)[0]
                #print("max indices = %d" % (m - idx))
                shift = (shift2[idx] * idx/shift2[idx].sum()).sum()
                #print("shift = ", shift - m, "in x units = ", (shift - m) * (x[1]-x[0]))

                # shift the curve
                shift = (shift - m) * (x[1]-x[0])
                x.shape = -1
                try:
                    y = numpy.fft.ifft(numpy.exp(-2.0*numpy.pi*numpy.sqrt(numpy.complex(-1))*\
                                numpy.fft.fftfreq(len(x), d=x[1]-x[0])*shift)*numpy.fft.fft(y))
                except ValueError:
                    print 'Shapes:'
                    print '\td = ',(x[1]-x[0])
                    print '\tlen(x) =',len(x)
                    print '\tlen(y) =',len(y)
                    return
                y = y.real
                y.shape = -1
            curveList.append([x, y, legend + " SHIFT", info, False, False])
            self.shiftDict[' '.join(legend.split(' ')[:-1])] = shift
        tmpArray = None
        curveList[-1][-2] = True
        curveList[-1][-1] = False
        if plot:
            x, y, legend, info, replot, replace = curveList[0]
            self.addCurve(x, y, legend=legend, replot=True, replace=True)
            for i in range(1, len(curveList)):
                x, y, legend, info, replot, replace = curveList[i]
                self.addCurve(x,
                              y,
                              legend=legend,
                              info=info,
                              replot=replot,
                              replace=False)
        return

        
        # now get the final spectrum
        y = medianSpectra.sum(axis=1) / nCurves
        x0.shape = -1
        y.shape = x0.shape
        legend = "%d Median from %s to %s" % (width,
                                              curves[0][2],
                                              curves[-1][2])
        self.addCurve(x0,
                      y,
                      legend=legend,
                      info=None,
                      replot=True,
                      replace=True)

MENU_TEXT = "Alignment Plugin"
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
