import numpy
from PyMca import PyMcaQt as qt
from PyMca import ScanWindow
from PyMca import specfilewrapper as sf
from PyMca import SimpleFitModule as SFM
from PyMca import SpecfitFunctions
from PyMca import Elements
from PyMca import ConfigDict

try:
    from PyMca import Plugin1DBase
except ImportError:
    print("WARNING:SumRulesPlugin import from somewhere else")
    from . import Plugin1DBase

DEBUG = 1

class Mathematics(object):
    def __init__(self):
        self.simpleFit = SFM.SimpleFit()
        self.simpleFit.importFunctions(SpecfitFunctions)
        
    def ricker(self, points, a):
        """
        SciPy implementation of the ricker wavelet
        
        From https://github.com/scipy/scipy/blob/v0.13.0/scipy/signal/wavelets.py
        """
        A = 2 / (np.sqrt(3 * a) * (np.pi**0.25))
        wsq = a**2
        vec = np.arange(0, points) - (points - 1.0) / 2
        tsq = vec**2
        mod = (1 - tsq / wsq)
        gauss = np.exp(-tsq / (2 * wsq))
        total = A * mod * gauss
        return total
    
    def continousWaveletTransform(self, data, widths):
        """
        SciPy implementation of cwt
        
        From https://github.com/scipy/scipy/blob/v0.13.0/scipy/signal/wavelets.py
        """
        wavelet = self.ricker
        nCols = len(data)
        nRows = len(widths)
        out = numpy.zeros([nRows, nCols], dtype=numpy.float)
        for idx, width in enumerate(widths):
            waveletData = wavelet(min(10 * width, len(data)), width)
            out[idx, :] = numpy.convolve(data, waveletData, mode='same')
        return out

    def estimateStuff(self, x, y, edgeList=[], preEdgeIdx=None, postEdgeIdx=None):
        ddict = {}
        if len(edgeList) == 0:
            # Track peaks
            pass
        else:
            # Edge positions given by the interface
            pass
        if preEdgeIdx:
            # Assume constant background in the pre edge region
            preEdgeLevel = y[0:preEdgeIdx]
        else:
            # Just average over the first 20 data points
            preEdgeLevel = y[0:20]
        if postEdgeIdx:
            # Assume constant background in the pre edge region
            preEdgeLevel = y[0:postEdgeIdx]
        else:
            # Just average over the first 20 data points
            preEdgeLevel = y[0:20]
            
            
        ddict['PreEdgeLevel']  = numpy.average(preEdgeLevel)
        ddict['PostEdgeLevel'] = numpy.average(postEdgeIdx)
        return ddict

    def cumtrapz(self, y, x=None, dx=1.0):
        y = y[:]
        if x is None:
            x = numpy.arange(len(y), dtype=y.dtype) * dx
        else:
            x = x[:]
        
        if not numpy.all(numpy.diff(x) > 0.):
            # assure monotonically increasing x
            idx = numpy.argsort(x)
            x = numpy.take(x, idx)
            y = numpy.take(y, idx)
            # Avoid dublicates
            x.ravel()
            idx = numpy.nonzero(numpy.diff(x) > 0)[0]
            x = numpy.take(x, idx)
            y = numpy.take(y, idx)
        
        return numpy.cumsum(.5 * numpy.diff(x) * (y[1:] + y[:-1]))
        
    def rndDataQuad(self):
        x = 5 + numpy.random.rand(500) + numpy.arange(500, dtype=float)
        y = 50*np.exp(-0.005*(x-250)**2) + 5 + 0.00005*x**2
        return x, y
        
    def rndDataLin(self):
        x = 5 + numpy.random.rand(500) + numpy.arange(500, dtype=float)
        y = 50*np.exp(-0.005*(x-250)**2) + 5 + 0.03*x
        return x, y        
        
    def detrend(self, x, y, order='linear'):
        if order not in ['linear', 'quadratic', 'cubic']:
            raise ValueError('Order must be linear, quadratic or cubic')
        if order == 'linear':
            ord = 1
        elif order == 'quadratic':
            ord = 2
        elif order == 'cubic':
            ord = 3
        coeff = numpy.polyfit(x,y,ord)
        poly  = numpy.zeros(x.shape)
        for a in coeff:
            poly *= x
            poly += a
        return y-poly
    
    def run(self):
        from matplotlib import pyplot as plt
        xLin, yLin   = self.rndDataLin()
        xQuad, yQuad = self.rndDataQuad()
        yLinCorr  = self.detrend(xLin, yLin, 'linear')
        yQuadCorr = self.detrend(xQuad, yQuad, 'quadratic')
        plt.plot(xLin, yLin, xLin, yLinCorr)
        plt.plot(xQuad, yQuad, xQuad, yQuadCorr)
        plt.show()

class MarkerSpinBox(qt.QDoubleSpinBox):
    
    valueChangedSignal = qt.pyqtSignal(object)
    
    def __init__(self, window, graph, label='', parent=None):
        qt.QDoubleSpinBox.__init__(self, parent)
        
        # Attributes
        self.window = window
        self.graph = graph
        self.markerID = self.graph.insertX1Marker(0., label=label)
        
        # Initialize
        self.setMinimum(0.)
        self.setMaximum(10000.)
        self.setValue(0.)
        
        # Connects
        self.connect(self.graph,
             qt.SIGNAL("QtBlissGraphSignal"),
             self._markerMoved)
        self.valueChanged['double'].connect(self._valueChanged)
        self.valueChanged['QString'].connect(self._valueChanged)

    def _setMarkerFollowMouse(self, windowTitle):
        windowTitle = str(windowTitle)
        if self.window == windowTitle:
            self.graph.setmarkerfollowmouse(self.markerID, True)
        else:
            self.graph.setmarkerfollowmouse(self.markerID, False)

    def _markerMoved(self, ddict):
        if 'marker' not in ddict:
            return
        else:
            if ddict['marker'] != self.markerID:
                return  
        if DEBUG:
            print "_markerMoved -- ddict:\n\t",ddict
        if ddict['event'] == 'markerMoving':
            self.setValue(ddict['x'])
            
    def _valueChanged(self, val):
        try:
            val = float(val)
        except ValueError:
            print '_valueChanged -- Sorry, it ain\'t gonna float:',val
            return
        self.graph.setMarkerXPos(self.markerID, val)
        self.graph.replot()
        
class LineEditDisplay(qt.QLineEdit):
    def __init__(self, combobox, parent=None):
        qt.QLineEdit.__init__(self, parent)
    
    

class SumRulesWindow(qt.QWidget):

    tabChangedSignal = qt.pyqtSignal('QString')
    tabList = ['element','xas', 'xmcd'] # TODO: expand to ['element','xas', 'xmcd']
    xasMarkerList  = ['Pre Min','Pre Max','Post Min','Post Max']
    xmcdMarkerList = ['p','q','r']
    
    # 3d
    transistionMetals = ['Ti', 'V', 'Cr', 'Mn',\
                         'Fe', 'Co', 'Ni', 'Cu']
    # 4f
    rareEarths = ['Ce', 'Pr', 'Nd', 'Pm', 'Sm',\
                  'Eu', 'Gd', 'Tb', 'Dy', 'Ho',\
                  'Er', 'Tm', 'Yb']
    elementsDict = {
            ''  : [],
            '3d': transistionMetals,
            '4f': rareEarths
    }
    # Electron configurations
    electronConfs = ['3d','4f']

    def __init__(self, parent=None):
        qt.QWidget.__init__(self, parent)
        self.setWindowTitle('Sum Rules')
        self.plotWindow = ScanWindow.ScanWindow(self)
        self.plotWindow.scanWindowInfoWidget.hide()
        self.plotWindow.toolBar.hide()
        self.plotWindow.graph.enablemarkermode()
        
        # Marker Handling
        # spinboxDict connects marker movement to spinbox
        # keys() -> id(MarkerSpinBox)
        # values() -> MarkerSpinBox
        self.spinboxDict = {}
        self.valuesDict = dict(
                    [(item, {}) for item in self.tabList])
        
        # Tab Widget
        self.tabWidget = qt.QTabWidget()
        for window in self.tabList:
            if window == 'element':
                # electron configuration combo box
                self.elementEConfCB = qt.QComboBox()
                self.elementEConfCB.addItems(['']+self.electronConfs)
                self.elementEConfCB.currentIndexChanged['QString'].connect(self.setElectronConf)
                elementEConfLayout = qt.QHBoxLayout()
                elementEConfLayout.setContentsMargins(0,0,0,0)
                elementEConfLayout.addWidget(qt.QLabel('Electron configuration'))
                elementEConfLayout.addWidget(qt.HorizontalSpacer())
                elementEConfLayout.addWidget(self.elementEConfCB)
                elementEConfWidget = qt.QWidget()
                elementEConfWidget.setLayout(elementEConfLayout)
                # element selection combo box
                self.elementCB = qt.QComboBox()
                self.elementCB.addItems([''])
                self.elementCB.currentIndexChanged['QString'].connect(self.getElementInfo)
                elementLayout = qt.QHBoxLayout()
                elementLayout.setContentsMargins(0,0,0,0)
#                elementLayout.addWidget(qt.QLabel('Electron configuration'))
#                elementLayout.addWidget(qt.HorizontalSpacer())
#                elementLayout.addWidget(self.elementEConfCB)
                elementLayout.addWidget(qt.QLabel('Element'))
                elementLayout.addWidget(qt.HorizontalSpacer())
                elementLayout.addWidget(self.elementCB)
                elementWidget = qt.QWidget()
                elementWidget.setLayout(elementLayout)
                # X-ray absorption edge selection combo box
                # TODO: Add LineEdits to display edge energies, c.f. LineEditDisplay class!
                self.edge1CB = qt.QComboBox()
                self.edge2CB = qt.QComboBox()
                self.edge1CB.addItems([''])
                self.edge2CB.addItems([''])
                #self.edgeCB.currentIndexChanged['QString'].connect()
                edge1Layout = qt.QHBoxLayout()
                edge1Layout.setContentsMargins(0,0,0,0)
                edge1Layout.addWidget(qt.QLabel('Edge 1'))
                edge1Layout.addWidget(qt.HorizontalSpacer())
                edge1Layout.addWidget(self.edge1CB)
                edge2Layout = qt.QHBoxLayout()
                edge2Layout.setContentsMargins(0,0,0,0)
                edge2Layout.addWidget(qt.QLabel('Edge 2'))
                edge2Layout.addWidget(qt.HorizontalSpacer())
                edge2Layout.addWidget(self.edge2CB)
                edge1Widget = qt.QWidget()
                edge1Widget.setLayout(edge1Layout)
                edge2Widget = qt.QWidget()
                edge2Widget.setLayout(edge2Layout)
                # electron occupation number
                self.electronOccupation = qt.QLineEdit('e.g. 3.14')
                electronOccupationValidator = qt.QDoubleValidator()
                electronOccupationValidator.setBottom(0.)
                electronOccupationValidator.setTop(14.)
                self.electronOccupation.setValidator(electronOccupationValidator)
                electronOccupationLayout = qt.QHBoxLayout()
                electronOccupationLayout.setContentsMargins(0,0,0,0)
                electronOccupationLayout.addWidget(qt.QLabel('Electron Occupation Number'))
                electronOccupationLayout.addWidget(qt.HorizontalSpacer())
                electronOccupationLayout.addWidget(self.electronOccupation)
                electronOccupationWidget = qt.QWidget()
                electronOccupationWidget.setLayout(electronOccupationLayout)
                # tab layouting
                elementTabLayout = qt.QVBoxLayout()
                elementTabLayout.addWidget(elementEConfWidget)
                elementTabLayout.addWidget(elementWidget)
                elementTabLayout.addWidget(electronOccupationWidget)
                elementTabLayout.addWidget(qt.QLabel('X-ray absorption edges'))
                elementTabLayout.addWidget(edge1Widget)
                elementTabLayout.addWidget(edge2Widget)
                elementTabLayout.addWidget(qt.VerticalSpacer())
                elementTabWidget = qt.QWidget()
                elementTabWidget.setLayout(elementTabLayout)
                # Add to self.valuesDict
                self.valuesDict['element']\
                        ['element'] = self.elementCB
                self.valuesDict['element']\
                        ['electron configuration'] = self.elementEConfCB
                self.valuesDict['element']\
                        ['electron occupation'] = self.electronOccupation
                self.valuesDict['element']['edge1'] = self.edge1CB
                self.valuesDict['element']['edge2'] = self.edge2CB
                self.valuesDict['element']['info'] = {}
                self.tabWidget.addTab(
                            elementTabWidget,
                            window.upper())
            elif window == 'xas':
                # Creat Pre/Post edge group box
                prePostLayout = qt.QVBoxLayout()
                prePostLayout.setContentsMargins(0,0,0,0)
                for markerLabel in self.xasMarkerList:
                    # TODO: Fix intial xpos
                    markerWidget, spinbox = self.addMarker(window=window,
                                                  label=markerLabel,
                                                  xpos=630.)
                    self.valuesDict[window][markerLabel] = spinbox
                    markerWidget.setContentsMargins(0,-8,0,-8)
                    prePostLayout.addWidget(markerWidget)
                prePostGB = qt.QGroupBox('Pre/Post edge')
                prePostGB.setLayout(prePostLayout)
                # Creat Edge group box
                # TODO: Determine number of edges
                numberOfEdges = 2
                addDelLayout = qt.QHBoxLayout()
                addDelLayout.setContentsMargins(0,0,0,0)
                buttonAdd = qt.QPushButton('Add')
                buttonDel = qt.QPushButton('Del')
                buttonAdd.clicked.connect(self.addEdgeMarker)
                buttonDel.clicked.connect(self.delEdgeMarker)
                addDelLayout.addWidget(qt.HorizontalSpacer())
                addDelLayout.addWidget(buttonAdd)
                addDelLayout.addWidget(buttonDel)
                addDelWidget = qt.QWidget()
                addDelWidget.setLayout(addDelLayout)
                edgeLayout = qt.QVBoxLayout()
                edgeLayout.setContentsMargins(0,0,0,0)
                edgeLayout.addWidget(addDelWidget)
                for i in range(numberOfEdges):
                    markerWidget, spinbox = self.addMarker(window=window,
                                                  label=markerLabel,
                                                  xpos=700.)
                    self.valuesDict[window][markerLabel] = spinbox
                    markerWidget.setContentsMargins(0,-8,0,-8)
                    edgeLayout.addWidget(markerWidget)
                edgeGB = qt.QGroupBox('Edge positions')
                edgeGB.setLayout(edgeLayout)
                # Insert into tab
                xasTabLayout = qt.QVBoxLayout()
                xasTabLayout.setContentsMargins(0,0,0,0)
                xasTabLayout.addWidget(prePostGB)
                xasTabLayout.addWidget(edgeGB)    
                xasTabLayout.addWidget(qt.VerticalSpacer())
                xasWidget = qt.QWidget()
                xasWidget.setLayout(xasTabLayout)
                self.tabWidget.addTab(
                            xasWidget,
                            window.upper())
            #markerWidget, spinbox = self.addMarker(window=window,label=window,xpos=700.)
            elif window == 'xmcd':
                pqLayout = qt.QVBoxLayout()
                pqLayout.setContentsMargins(0,0,0,0)
                for markerLabel in self.xmcdMarkerList:
                    # TODO: Fix intial xpos
                    markerWidget, spinbox = self.addMarker(window=window,
                                                  label=markerLabel,
                                                  xpos=800.)
                    self.valuesDict[window][markerLabel] = spinbox
                    markerWidget.setContentsMargins(0,-8,0,-8)
                    pqLayout.addWidget(markerWidget)
                pqGB = qt.QGroupBox('XMCD integrals')
                pqGB.setLayout(pqLayout)
                xmcdTabLayout = qt.QVBoxLayout()
                xmcdTabLayout.addWidget(pqGB)
                xmcdTabLayout.addWidget(qt.VerticalSpacer())
                xmcdWidget = qt.QWidget()
                xmcdWidget.setLayout(xmcdTabLayout)
                self.tabWidget.addTab(
                            xmcdWidget,
                            window.upper())
            #self.tabWidget.addTab(markerWidget, window.upper())
        self.tabWidget.currentChanged['int'].connect(
                            self._handleTabChangedSignal)
        
        
        # Add/Remove marker Buttons
        #buttonAddMarker = qt.QPushButton('Add')
        #buttonDelMarker = qt.QPushButton('Del')
        buttonPrint = qt.QPushButton('Print ElementInfo')
        buttonPrint.clicked.connect(self.printElementInfo)
        self.plotWindow.graphBottomLayout.addWidget(qt.HorizontalSpacer())
        #self.plotWindow.graphBottomLayout.addWidget(buttonAddMarker)
        #self.plotWindow.graphBottomLayout.addWidget(buttonDelMarker)
        self.plotWindow.graphBottomLayout.addWidget(buttonPrint)
        
        # Layout
        mainLayout = qt.QVBoxLayout()
        mainLayout.addWidget(self.plotWindow)
        mainLayout.addWidget(self.tabWidget)
        mainLayout.setContentsMargins(1,1,1,1)
        
        # Data handling
        self.xmcdData = None
        self.xasData  = None
        
        self.setLayout(mainLayout)
        tmpDict = {
                'xas' : {
                    'Pre Min': 630,
                    'Pre Max': 650,
                    'Post Min': 670,
                    'Post Max': 690
                },
                'xmcd': {
                    'p': 750,
                    'q': 770
                },
                'element': {
                    'electron configuration': '3d',
                    'element': 'Fe'
                }
        }
        self.setValuesDict(tmpDict)
        print self.getValuesDict()

    def setElectronConf(self, eConf):
        # updates the element combo box
        eConf = str(eConf)
        if len(eConf) == 0:
            self.electronOccupation.setDisabled(True)
        else:
            self.electronOccupation.setDisabled(False)
        self.elementCB.clear()
        elementsList = self.elementsDict[eConf]
        self.elementCB.addItems(['']+elementsList)

    def getElementInfo(self, symbol):
        #eConfTmp = self.valuesDict['element']\
        #               ['Electron Configuration']
        #self.valuesDict['element']['element']
        #self.valuesDict['element']\
        #        ['Electron Configuration'] = eConfTmp
        ddict = {}
        symbol = str(symbol)
        try:
            ddict = Elements.Element[symbol]
        except KeyError:
            msg  = ('setElement -- %s not found in '%symbol)
            msg += 'Elements.Element dictionary'
            print(msg)
        self.valuesDict['element']['info'] = ddict
    
    def printElementInfo(self):
        if len(self.valuesDict['element']['info']) == 0:
            print 'self.valuesDict[\'element\'][\'info\'] is empty'
            return
        for k,v in self.valuesDict['element']['info'].items():
            print k,':',v

    def getCurrentTab(self):
        idx = self.tabWidget.currentIndex()
        return self.tabList[idx]

    def getValuesDict(self):
        ddict = {}
        for tab, tabDict in self.valuesDict.items():
            if tab not in ddict.keys():
                ddict[tab] = {}
            for key, obj in tabDict.items():
                value = None
                if isinstance(obj, MarkerSpinBox):
                    value = obj.value()
                if isinstance(obj, qt.QComboBox):
                    tmp = obj.currentText()
                    value = str(tmp)
                ddict[tab][key] = value
        return ddict

    def setValuesDict(self, ddict):
        markerList  = (self.xasMarkerList + self.xmcdMarkerList)
        elementList = (self.transistionMetals 
                       + self.rareEarths 
                       + self.electronConfs)
        for tab, tabDict in ddict.items():
            if tab not in self.valuesDict.keys():
                raise KeyError('setValuesDict -- Tab not found')
            for key, value in tabDict.items():
                if not isinstance(key, str):
                    raise KeyError('setValuesDict -- key is not str instance')
                obj = self.valuesDict[tab][key]
                if isinstance(obj, MarkerSpinBox):
                    try:
                        tmp = float(value)
                    except ValueError:
                        xmin, xmax = self.plotWindow.graph.getX1AxisLimits()
                        tmp = xmin + (xmax-xmin)/10.
                        msg  = 'setValuesDict -- Float conversion failed'
                        msg += ' while setting marker positions. Value:', value
                        print(msg)
                    obj.setValue(tmp)
                elif isinstance(obj, qt.QComboBox):
                    idx = obj.findText(qt.QString(value))
                    obj.setCurrentIndex(idx)
                else:
                    raise KeyError('setValuesDict -- \'%s\' not found'%key)
                
    def addEdgeMarker(self):
        print 'addEdgeMarker clicked'
        
    def delEdgeMarker(self):
        print 'delEdgeMarker clicked'

    def setData(self, x, y, identifier, xlabel='ene_st', ylabel='zratio'):
        if identifier not in ['xmcd', 'xas']:
            raise ValueError('Identifier must either be \'xmcd\' or \'xas\'!')
        # Sort energy range
        sortedIdx = x.argsort()
        xSorted = x.take(sortedIdx)[:]
        ySorted = y.take(sortedIdx)[:]
        # Ensure strictly monotonically increasing energy range
        dx = numpy.diff(x)
        if not numpy.all(dx > 0.):
            mask = numpy.nonzero(dx)
            xSorted = numpy.take(xSorted, mask)
            ySorted = numpy.take(ySorted, mask)
        # Add spectrum to plotWindow using the 
        self.plotWindow.newCurve(
                x=xSorted, 
                y=ySorted,
                legend=identifier,
                xlabel=xlabel, 
                ylabel=ylabel, 
                info={}, 
                replot=False, 
                replace=False)
        specLegend = self.plotWindow.dataObjectsList[-1]
        # Calculate the cumulative intergral
#        mathObj = Mathematics()
#        yInt = mathObj.cumtrapz(y=ySorted, x=xSorted)
#        xInt = .5 * (xSorted[1:] + xSorted[:-1])
#        # Add integral to plotWindow
#        self.plotWindow.newCurve(
#                x=xInt, 
#                y=yInt,
#                legend=identifier+' Integral',
#                xlabel=xlabel, 
#                ylabel=ylabel, 
#                info={}, 
#                replot=False, 
#                replace=False)
#        intLegend = self.plotWindow.dataObjectsList[-1]
        if identifier == 'xmcd':
            self.xmcdData = self.plotWindow.dataObjectsDict[specLegend]
            self.plotWindow.graph.mapToY2(specLegend)
            #self.plotWindow.graph.mapToY2(intLegend)
        elif identifier == 'xas':
            self.xasData  = self.plotWindow.dataObjectsDict[specLegend]
        
    def addMarker(self, window, label='X MARKER', xpos=None):
        # Add spinbox controlling the marker
        graph = self.plotWindow.graph
        spinbox = MarkerSpinBox(window, graph, label)

        # Connects
        self.tabChangedSignal.connect(spinbox._setMarkerFollowMouse)
        
        # Widget & Layout
        spinboxWidget = qt.QWidget()
        spinboxLayout = qt.QHBoxLayout()
        spinboxLayout.addWidget(qt.QLabel(label))
        spinboxLayout.addWidget(qt.HorizontalSpacer())
        spinboxLayout.addWidget(spinbox)
        spinboxWidget.setLayout(spinboxLayout)
        
        #self.spinboxDict[id(spinbox)] = spinbox
        print 'spinboxID:',id(spinbox)
#        graph.replot()
        
        return spinboxWidget, spinbox

    def _handleTabChangedSignal(self, idx):
        if idx >= len(self.tabList):
            print 'Tab changed -- idx:',idx,'..Abort'
            return
        tab = self.tabList[idx]
        print 'Tab changed -- idx:',idx,'tab:',tab
        self.tabChangedSignal.emit(tab)
        

class SumRulesPlugin(Plugin1DBase.Plugin1DBase):
    def __init__(self, plotWindow, **kw):
        Plugin1DBase.Plugin1DBase.__init__(self, plotWindow, **kw)
        self.__randomization = True
        self.__methodKeys = []
        self.methodDict = {}

        function = self.dummy
        method = "Dummy function"
        text = "Does nothing"
        info = text
        icon = None
        self.methodDict[method] = [function,
                                   info,
                                   icon]
        self.__methodKeys.append(method)
        
        self._widget = None
        
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

    def showWidget(self):
        if self._widget is None:
            self._widget = SumRulesWindow()
            self._widget.show()
        else:
            self._widget.show()
            self._widget.raise_()

MENU_TEXT = "Sum Rules Plugin"
def getPlugin1DInstance(plotWindow, **kw):
    ob = AlignmentScanPlugin(plotWindow)
    return ob

def getData(fn='/home/truter/lab/datasets/sum_rules/Fe_L23/xld_analysis.spec'):
    analysis = sf.specfile.Specfile(fn)
    xmcdArr = []
    xasArr  = []
    avgA, avgB = [], []
    spec = analysis[0]
    x = spec[0][:]
    avgA = spec[1][:]
    avgB = spec[2][:]
    xmcdArr = spec[3][:]
    xasArr  = spec[4][:]
    return x, avgA, avgB, xmcdArr, xasArr

if __name__ == '__main__':
   
    app = qt.QApplication([])
    win = SumRulesWindow()
    x, avgA, avgB, xmcd, xas = getData()
    #win.plotWindow.newCurve(x,xmcd, legend='xmcd', xlabel='ene_st', ylabel='zratio', info={}, replot=False, replace=False)
    win.setData(x,xmcd, identifier='xmcd')
    #win.plotWindow.newCurve(x,xas, legend='xas', xlabel='ene_st', ylabel='zratio', info={}, replot=False, replace=False)
    win.setData(x,xas, identifier='xas')
    win.show()
    app.exec_()
