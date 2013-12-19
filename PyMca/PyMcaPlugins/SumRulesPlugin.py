import numpy
from PyMca import PyMcaQt as qt
from PyMca import ScanWindow
from PyMca import specfilewrapper as sf
from PyMca import SimpleFitModule as SFM
from PyMca import SpecfitFunctions
from PyMca import Elements
from PyMca import ConfigDict
from PyMca.SpecfitFuns import upstep, downstep
from PyMca.Gefit import LeastSquaresFit as LSF

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
        A = 2 / (numpy.sqrt(3 * a) * (numpy.pi**0.25))
        wsq = a**2
        vec = numpy.arange(0, points) - (points - 1.0) / 2
        tsq = vec**2
        mod = (1 - tsq / wsq)
        gauss = numpy.exp(-tsq / (2 * wsq))
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
        y = 50*numpy.exp(-0.005*(x-250)**2) + 5 + 0.00005*x**2
        return x, y
        
    def rndDataLin(self):
        x = 5 + numpy.random.rand(500) + numpy.arange(500, dtype=float)
        y = 50*numpy.exp(-0.005*(x-250)**2) + 5 + 0.03*x
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
            self.graph.setmarkercolor(self.markerID, 'blue')
            self.graph.setmarkerfollowmouse(self.markerID, True)
            self.graph.replot()
        else:
            self.graph.setmarkercolor(self.markerID, 'black')
            self.graph.setmarkerfollowmouse(self.markerID, False)
            self.graph.replot()

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
    def __init__(self, combobox, ddict={}, parent=None):
        qt.QLineEdit.__init__(self, parent)
        self.setReadOnly(True)
        self.setAlignment(qt.Qt.AlignRight)
        self.ddict = ddict
        self.setMaximumWidth(120)
        self.combobox = combobox
        self.combobox.currentIndexChanged['QString'].connect(self.setText)
        #self.combobox.destroyed.connect(self.destroy)
        
    def updateDict(self, ddict):
        self.ddict = ddict
        
    def checkComboBox(self):
        tmp = self.combobox.currentText()
        self.setText(tmp)
        
    def setText(self, inp):
        inp = str(inp)
        if inp == '':
            text = ''
        else:
            tmp = self.ddict.get(inp,None)
            if tmp is not None:
                try:
                    text = '%.2f meV'%(1000. * float(tmp))
                except ValueError:
                    text = 'NaN' 
            else:
                text = '---'
        qt.QLineEdit.setText(self, text)
    

class SumRulesWindow(qt.QWidget):

    tabChangedSignal = qt.pyqtSignal('QString')
    tabList = ['element','background', 'integration']
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
                sampleGB         = qt.QGroupBox('Sample definition')
                sampleLayout     = qt.QVBoxLayout()
                sampleGB.setLayout(sampleLayout)
                
                absorptionGB     = qt.QGroupBox('X-ray absorption edges')
                absorptionLayout = qt.QVBoxLayout()
                absorptionGB.setLayout(absorptionLayout)
                
                # BEGIN sampleGB
                # electron configuration combo box
                self.elementEConfCB = qt.QComboBox()
                self.elementEConfCB.setMinimumWidth(100)
                self.elementEConfCB.addItems(['']+self.electronConfs)
                self.elementEConfCB.currentIndexChanged['QString'].connect(self.setElectronConf)
                elementEConfLayout = qt.QHBoxLayout()
                elementEConfLayout.setContentsMargins(0,0,0,0)
                elementEConfLayout.addWidget(qt.QLabel('Electron configuration'))
                elementEConfLayout.addWidget(qt.HorizontalSpacer())
                elementEConfLayout.addWidget(self.elementEConfCB)
                elementEConfWidget = qt.QWidget()
                elementEConfWidget.setLayout(elementEConfLayout)
                sampleLayout.addWidget(elementEConfWidget)
                # element selection combo box
                self.elementCB = qt.QComboBox()
                self.elementCB.setMinimumWidth(100)
                self.elementCB.addItems([''])
                self.elementCB.currentIndexChanged['QString'].connect(self.getElementInfo)
                elementLayout = qt.QHBoxLayout()
                elementLayout.setContentsMargins(0,0,0,0)
                elementLayout.addWidget(qt.QLabel('Element'))
                elementLayout.addWidget(qt.HorizontalSpacer())
                elementLayout.addWidget(self.elementCB)
                elementWidget = qt.QWidget()
                elementWidget.setLayout(elementLayout)
                sampleLayout.addWidget(elementWidget)
                # electron occupation number
                self.electronOccupation = qt.QLineEdit('e.g. 3.14')
                self.electronOccupation.setMaximumWidth(120)
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
                sampleLayout.addWidget(electronOccupationWidget)
                # END sampleGB
                
                # BEGIN absorptionGB: X-ray absorption edge 
                # selection combo box by transition (L3M1, etc.)
                self.edge1CB = qt.QComboBox()
                self.edge1CB.setMinimumWidth(100)
                self.edge1CB.addItems([''])
                self.edge1Line = LineEditDisplay(self.edge1CB)
                edge1Layout = qt.QHBoxLayout()
                edge1Layout.setContentsMargins(0,0,0,0)
                edge1Layout.addWidget(qt.QLabel('Edge 1'))
                edge1Layout.addWidget(qt.HorizontalSpacer())
                edge1Layout.addWidget(self.edge1CB)
                edge1Layout.addWidget(self.edge1Line)
                edge1Widget = qt.QWidget()
                edge1Widget.setLayout(edge1Layout)
                absorptionLayout.addWidget(edge1Widget)
                
                self.edge2CB = qt.QComboBox()
                self.edge2CB.setMinimumWidth(100)
                self.edge2CB.addItems([''])
                self.edge2Line = LineEditDisplay(self.edge2CB)
                edge2Layout = qt.QHBoxLayout()
                edge2Layout.setContentsMargins(0,0,0,0)
                edge2Layout.addWidget(qt.QLabel('Edge 2'))
                edge2Layout.addWidget(qt.HorizontalSpacer())
                edge2Layout.addWidget(self.edge2CB)
                edge2Layout.addWidget(self.edge2Line)
                edge2Widget = qt.QWidget()
                edge2Widget.setLayout(edge2Layout)
                absorptionLayout.addWidget(edge2Widget)
                # END absorptionGB
                
                # BEGIN tab layouting
                elementTabLayout = qt.QVBoxLayout()
                #elementTabLayout.addWidget(elementEConfWidget)
                #elementTabLayout.addWidget(elementWidget)
                #elementTabLayout.addWidget(electronOccupationWidget)
                #elementTabLayout.addWidget(qt.QLabel('X-ray absorption edges'))
                #elementTabLayout.addWidget(edge1Widget)
                #elementTabLayout.addWidget(edge2Widget)
                elementTabLayout.addWidget(sampleGB)
                elementTabLayout.addWidget(absorptionGB)
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
                self.valuesDict['element']['edge1Transistion'] = self.edge1CB
                self.valuesDict['element']['edge2Transistion'] = self.edge2CB
                self.valuesDict['element']['edge1Energy'] = self.edge1Line
                self.valuesDict['element']['edge2Energy'] = self.edge2Line
                self.valuesDict['element']['info'] = {}
                self.tabWidget.addTab(
                            elementTabWidget,
                            window.upper())
                # END tab layouting
            elif window == 'background':
                # BEGIN Pre/Post edge group box
                prePostLayout = qt.QGridLayout()
                prePostLayout.setContentsMargins(0,0,0,0)
                for idx, markerLabel in enumerate(self.xasMarkerList):
                    # TODO: Fix intial xpos
                    markerWidget, spinbox = self.addMarker(window=window,
                                                  label=markerLabel,
                                                  xpos=0.)
                    self.valuesDict[window][markerLabel] = spinbox
                    markerWidget.setContentsMargins(0,-8,0,-8)
                    if idx == 0: posx, posy = 0,0
                    if idx == 1: posx, posy = 1,0
                    if idx == 2: posx, posy = 0,1
                    if idx == 3: posx, posy = 1,1                    
                    prePostLayout.addWidget(markerWidget, posx, posy)
                prePostGB = qt.QGroupBox('Pre/Post edge')
                prePostGB.setLayout(prePostLayout)
                # END Pre/Post edge group box
                
                # BEGIN Edge group box
                numberOfEdges = 2
                #addDelLayout = qt.QHBoxLayout()
                #addDelLayout.setContentsMargins(0,0,0,0)
                #buttonAdd = qt.QPushButton('Add')
                #buttonDel = qt.QPushButton('Del')
                #buttonAdd.clicked.connect(self.addEdgeMarker)
                #buttonDel.clicked.connect(self.delEdgeMarker)
                #addDelLayout.addWidget(qt.HorizontalSpacer())
                #addDelLayout.addWidget(buttonAdd)
                #addDelLayout.addWidget(buttonDel)
                #addDelWidget = qt.QWidget()
                #addDelWidget.setLayout(addDelLayout)
                edgeLayout = qt.QVBoxLayout()
                edgeLayout.setContentsMargins(0,0,0,0)
                #edgeLayout.addWidget(addDelWidget)
                for idx in range(numberOfEdges):
                    markerLabel = 'Edge %d'%(idx+1)
                    markerWidget, spinbox = self.addMarker(window=window,
                                                  label=markerLabel,
                                                  xpos=0.)
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
                # END Edge group box
            #markerWidget, spinbox = self.addMarker(window=window,label=window,xpos=700.)
            elif window == 'integration':
                pqLayout = qt.QVBoxLayout()
                pqLayout.setContentsMargins(0,0,0,0)
                for markerLabel in self.xmcdMarkerList:
                    # TODO: Fix intial xpos
                    markerWidget, spinbox = self.addMarker(window=window,
                                                  label=markerLabel,
                                                  xpos=0.)
                    self.valuesDict[window][markerLabel] = spinbox
                    markerWidget.setContentsMargins(0,-8,0,-8)
                    pqLayout.addWidget(markerWidget)
                pqGB = qt.QGroupBox('XAS/XMCD integrals')
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
        buttonPrint = qt.QPushButton('Estimate')
        #buttonPrint.clicked.connect(self.estimatePrePostEdgePositions)
        buttonPrint.clicked.connect(self.estimateBG)
        self.plotWindow.graphBottomLayout.addWidget(qt.HorizontalSpacer())
        #self.plotWindow.graphBottomLayout.addWidget(buttonAddMarker)
        #self.plotWindow.graphBottomLayout.addWidget(buttonDelMarker)
        self.plotWindow.graphBottomLayout.addWidget(buttonPrint)
        
        # Layout
        mainLayout = qt.QVBoxLayout()
        mainLayout.addWidget(self.plotWindow)
        mainLayout.addWidget(self.tabWidget)
        mainLayout.setContentsMargins(1,1,1,1)
        
        # Data handling:
        # Each is Tuple (x,y)
        # type(x),type(y) == ndarray
        self.xmcdData    = None
        self.xasData     = None
        self.xasDataCorr = None
        self.xasDataBG   = None
        
        self.setLayout(mainLayout)
        tmpDict = {
                'background' : {
                    'Pre Min': 658.02,
                    'Pre Max': 703.75,
                    'Post Min': 730.5,
                    'Post Max': 808.7,
                    'Edge 1': 721.44,
                    'Edge 2': 708.7,
                },
                'integration': {
                    'p': 900,
                    'q': 900
                },
                'element': {
                    'electron configuration': '3d',
                    'element': 'Fe',
                    'edge1Transistion': 'L3M4',
                    'edge2Transistion': 'L2M4'
                }
        }
        self.setValuesDict(tmpDict)

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
        if len(symbol) == 0:
            self.valuesDict['element']['info'] = {}
            return
        try:
            ddict = Elements.Element[symbol]
        except KeyError:
            msg  = ('setElement -- \'%s\' not found in '%symbol)
            msg += 'Elements.Element dictionary'
            print(msg)
        # Update valuesDict
        self.valuesDict['element']['info'] = ddict
        # Update the EdgeCBs
        # Lookup all keys ending in 'xrays'
        keys = [item for item in ddict.keys() if item.endswith('xrays')]
        keys.sort()
        # keys is list of list, flatten it..
        transitions = sum([ddict[key] for key in keys],[])
        tmpDict = dict( [(transition, ddict[transition]['energy']) for transition in transitions])
        for cb, ed in [(self.edge1CB, self.edge1Line),
                       (self.edge2CB, self.edge2Line)]:
            curr = cb.currentText()
            cb.clear()
            ed.clear()
            ed.updateDict(tmpDict)
            cb.addItems(['']+transitions)
            # Try to set to old entry
            idx = cb.findText(qt.QString(curr))
            if idx < 0: idx = 0
            cb.setCurrentIndex(idx)

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
                elif isinstance(obj, qt.QComboBox):
                    tmp = obj.currentText()
                    value = str(tmp)
                elif isinstance(obj, LineEditDisplay) or\
                     isinstance(obj, qt.QLineEdit):
                    tmp = obj.text()
                    value = str(tmp)
                elif isinstance(obj, dict):
                    value = obj
                ddict[tab][key] = value
        return ddict

    def setValuesDict(self, ddict):
        markerList  = (self.xasMarkerList + self.xmcdMarkerList)
        elementList = (self.transistionMetals 
                       + self.rareEarths 
                       + self.electronConfs)
        # Check as early as possible if element symbol is present
        try:
            symbol = ddict['element']['element']
            self.getElementInfo(symbol)
        except KeyError:
            pass
        for tab, tabDict in ddict.items():
            if tab not in self.valuesDict.keys():
                raise KeyError('setValuesDict -- Tab not found')
            succession = []
            for key, value in tabDict.items():
                if not isinstance(key, str):
                    raise KeyError('setValuesDict -- Key is not str instance')
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
                elif isinstance(obj, LineEditDisplay):
                    # Must be before isinstance(obj, qt.QLineEdit)
                    # since LineEditDisplay inherits from QLineEdit
                    obj.checkComboBox()
                elif isinstance(obj, qt.QLineEdit):
                    if value:
                        obj.setText(value) 
                    else:
                        obj.setText('')
                elif isinstance(obj, dict):
                    obj = value
                else:
                    raise KeyError('setValuesDict -- \'%s\' not found'%key)
                
    def addEdgeMarker(self):
        print 'addEdgeMarker clicked'
        
    def delEdgeMarker(self):
        print 'delEdgeMarker clicked'

    def setRawData(self, x, y, identifier):
        if identifier not in ['xmcd', 'xas']:
            msg  = 'Identifier must either be \'xmcd\' or \'xas\''
            raise ValueError(msg)
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
        if identifier == 'xmcd':
            self.xmcdData = (xSorted, ySorted)
            #self.plotWindow.graph.mapToY2(intLegend)
        elif identifier == 'xas':
            self.xasData  = (xSorted, ySorted)
        # Trigger replot when data is added
        currIdx = self.tabWidget.currentIndex()
        self._handleTabChangedSignal(currIdx)

    def estimatePrePostEdgePositions(self):
    #def estimatePrePostEdgePositions(self, x, y, edgeList=[]):
        if self.xasData is None:
            return

        ddict = self.getValuesDict()
        edgeList = [ddict['element']['edge1Energy'],
                    ddict['element']['edge2Energy']]
        edgeList = [float(tmp.replace('meV','')) for tmp in edgeList
                    if len(tmp)>0 and tmp!='---']
        x, y = self.xasData
        xLimMin, xLimMax = self.plotWindow.getGraphXLimits()
            
        xMin = x[0]
        xMax = x[-1]
        xLen = xMax - xMin
        xMiddle = .5 *(xMax + xMin)
        # Average step length (Watch out for uneven data!)
        xStep = (xMax + xMin) / float(len(x))
        # Look for the index closest to the physical middle
        mask = numpy.nonzero(x <= xMiddle)[0]
        idxMid = mask[-1]

        factor = 20./100.
        if len(edgeList) == 0.:
            preMax = xMiddle - factor*xLen
            preMin = xMiddle + factor*xLen
            edge1  = xMiddle
            edge2  = xMiddle
        elif len(edgeList) == 1:
            edge = edgeList[0]
            preMax = edge - factor*xLen
            preMin = edge + factor*xLen
            edge1  = edge
            edge2  = edge
        else:
            edge1 = min(edgeList)
            edge2 = max(edgeList)
            preMax = edge1 - factor*xLen
            preMin = edge2 + factor*xLen

        ddict['background']['Pre Min']  = max(xMin,xLimMin+xStep)
        ddict['background']['Pre Max']  = preMax
        ddict['background']['Post Min'] = preMin
        ddict['background']['Post Max'] = min(xMax,xLimMax-xStep)
        ddict['background']['Edge 1'] = edge1
        ddict['background']['Edge 2'] = edge2
        
        self.setValuesDict(ddict)

    def plotOnDemand(self, window, xlabel='ene_st', ylabel='zratio'):
        # TODO: Errorhandling if self.xas, self.xmcd are none
        # Remove all curves
        legends = self.plotWindow.getAllCurves(just_legend=True)
        self.plotWindow.removeCurves(legends, replot=False)
        if (self.xmcdData is None) or (self.xasData is None):
            # Nothing to do
            return
        xyList  = []
        mapToY2 = False
        window = window.lower()
        if window == 'element':
            xmcdX, xmcdY = self.xmcdData
            xasX,  xasY  = self.xasData
            xyList = [(xmcdX, xmcdY, 'xmcd'), 
                      (xasX, xasY, 'xas')]
            # At least one of the curve is going
            # to get plotted on secondary y axis
            mapToY2 = True 
        elif window == 'background':
            xasX, xasY= self.xasData
            xyList = [(xasX, xasY, 'xas')]
            if self.xasDataBG is not None:
                xasBGX, xasBGY = self.xasDataBG
                xyList += [(xasBGX, xasBGY, 'xas Background')]
        elif window == 'integration':
            if self.xasDataCorr is None:
                print 'plotOnDemand -- xasDataCorr is None. Calculate corrected XAS data first!'
                return
            mathObj = Mathematics()
            xmcdX, xmcdY = self.xmcdData
            xasX,  xasY  = self.xasData
            xmcdIntY = mathObj.cumtrapz(y=xmcdY, x=xmcdX)
            xmcdIntX = .5 * (xmcdX[1:] + xmcdX[:-1])
            
            xasIntY  = mathObj.cumtrapz(y=xasY,  x=xasX)
            xasIntX  = .5 * (xasX[1:] + xasX[:-1])
            ylabel += ' integral'
            xyList = [(xmcdIntX, xmcdIntY, 'xmcd'),
                      (xasIntX,  xasIntY,  'xas')]
        for x,y,legend in xyList:
            self.plotWindow.newCurve(
                    x=x, 
                    y=y,
                    legend=legend,
                    xlabel=xlabel, 
                    ylabel=ylabel, 
                    info={}, 
                    replot=False, 
                    replace=False)
            if mapToY2:
                specLegend = self.plotWindow.dataObjectsList[-1]
                self.plotWindow.graph.mapToY2(specLegend)
                mapToY2 = False
        self.plotWindow.graph.replot()
        
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

    def estimateBG(self):
        if self.xasData is None:
            return
        if self.tabWidget.currentIndex() != 1:
            return
            
        x, y = self.xasData
        #self.estimatePrePostEdgePositions()
        ddict = self.getValuesDict()
        x01 = ddict['background']['Edge 1']
        x02 = ddict['background']['Edge 2']
        preMin = ddict['background']['Pre Min']
        preMax = ddict['background']['Pre Max']
        postMin = ddict['background']['Post Min']
        postMax = ddict['background']['Post Max']
        
        idxPre  = numpy.nonzero((preMin <= x) & (x <= preMax))[0]
        idxPost = numpy.nonzero((postMin <= x) & (x <= postMax))[0]
        
        xPreMax  = x[idxPre.max()]
        xPostMin = x[idxPost.min()]
        gap = abs(xPreMax - xPostMin)
        
        avgPre  = numpy.average(y[idxPre])
        avgPost = numpy.average(y[idxPost])
        bottom  = min(avgPre,avgPost)
        top     = max(avgPre,avgPost)
        if avgPost >= avgPre:
            sign = 1.
            erf  = upstep
        else:
            sign = -1.
            erf  = downstep
        diff = abs(avgPost - avgPre)
        
        ymin = y.min()
        ymax = y.max()
        par1 = (2./3., x01, gap/4.)
        par2 = (1./3., x02, gap/6.)
        
        print (xPreMax)
        print (xPostMin)
        print (gap)
        
        print (avgPre)
        print (avgPost)
        print (bottom)
        print (top)
        print (sign)
        print (diff)
        
        model = bottom + sign * diff * (erf(par1, x) + erf(par2, x))
        preModel  = numpy.asarray(len(x)*[avgPre])
        postModel = numpy.asarray(len(x)*[avgPost])
        
        self.plotWindow.addCurve(x,model,'BG model',{})
        self.plotWindow.addCurve(x,preModel,'Pre BG model',{})
        self.plotWindow.addCurve(x,postModel,'Post BG model',{})
        

    def _handleTabChangedSignal(self, idx):
        if idx >= len(self.tabList):
            print 'Tab changed -- idx:',idx,'..Abort'
            return
        tab = self.tabList[idx]
        print 'Tab changed -- idx:',idx,'tab:',tab
        self.plotOnDemand(window=tab)
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
    win.setRawData(x,xmcd, identifier='xmcd')
    #win.plotWindow.newCurve(x,xas, legend='xas', xlabel='ene_st', ylabel='zratio', info={}, replot=False, replace=False)
    win.setRawData(x,xas, identifier='xas')
    win.show()
    app.exec_()
