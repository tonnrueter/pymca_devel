import numpy
from PyMca import PyMcaQt as qt
from PyMca import ScanWindow
from PyMca import specfilewrapper as sf

try:
    from PyMca import Plugin1DBase
except ImportError:
    print("WARNING:SumRulesPlugin import from somewhere else")
    from . import Plugin1DBase

DEBUG = 1
class DoubleSpinBox(qt.QDoubleSpinBox):
    
    valueChangedSignal = qt.pyqtSignal(object)
    
    def __init__(self, marker=None, parent=None):
        qt.QDoubleSpinBox.__init__(self, parent)
        self.markerID = marker
        # Connects
        self.valueChanged['double'].connect(self.valueChangedEmitter)
        self.valueChanged['QString'].connect(self.valueChangedEmitter)
        
    def valueChangedEmitter(self, value):
        ddict = {}
        if isinstance(value, qt.QString):
            value = float(value)
        ddict['emitter'] = id(self)
        ddict['value'] = value
        self.valueChangedSignal.emit(ddict)
        
    def _markerMoved(self, ddict):
        if 'marker' not in ddict:
            return
        else:
            markerID = ddict['marker']
        if DEBUG:
            print(ddict)
        if ddict['event'] == 'markerMoving':
            self.setValue(ddict['x'])

class SumRulesWindow(qt.QWidget):
    
#    QtBlissGraphSignal = qt.SIGNAL("QtBlissGraphSignal")
    
    def __init__(self, parent=None):
        qt.QWidget.__init__(self, parent)
        self.plotWindow = ScanWindow.ScanWindow(self)
        
        # Marker Handling
        self.spinboxDict = {}
        # Marker can be clicked
        self.plotWindow.graph.enablemarkermode()
        
        # Layout
        mainLayout = qt.QVBoxLayout()
        mainLayout.addWidget(self.plotWindow)
        
        # Data handling
        xasList, xmcdList = [], []
        
        self.setLayout(mainLayout)

    def addMarkerBox(self, label='X MARKER', xpos=None):
        xmin, xmax = self.plotWindow.getGraphXLimits()
        if xpos is None:
            xpos = .5*(xmin + xmax)
        marker = self.plotWindow.graph.insertX1Marker(xpos, label=label)
        # Marker can be moved when clicked
        self.plotWindow.graph.setmarkerfollowmouse(marker, True)
        
        spinbox = DoubleSpinBox(marker)
        spinbox.setMaximum(xmax)
        spinbox.setMinimum(xmin)
        spinbox.setValue(xpos)
        spinbox.valueChangedSignal.connect(self._moveMarker)
        self.connect(self.plotWindow.graph,
                     qt.SIGNAL("QtBlissGraphSignal"),
                     spinbox._markerMoved)
        
        # Widget & Layout
        spinboxWidget = qt.QWidget()
        spinboxLayout = qt.QHBoxLayout()
        spinboxLayout.addWidget(qt.QLabel(label))
        spinboxLayout.addWidget(qt.HorizontalSpacer())
        spinboxLayout.addWidget(spinbox)
        spinboxWidget.setLayout(spinboxLayout)
        
        self.spinboxDict[id(spinbox)] = spinbox
        print 'spinboxID:',id(spinbox)
        print 'markerID:',marker
        return spinboxWidget

    def _moveMarker(self, ddict):
        emitter = ddict['emitter'] # = id(spinbox)
        xpos    = ddict['value']
        spinbox = self.spinboxDict[emitter]
        print 'type(xpos) =',type(xpos)
        print 'xpos =',xpos
        self.plotWindow.graph.setMarkerXPos(spinbox.markerID, xpos)

    def cumtrapz(y, x=None, dx=1.0):
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
            idx = numpy.nonzero(np.diff(x) > 0)[0]
            x = numpy.take(x, idx)
            y = numpy.take(y, idx)
        
        return numpy.cumsum(.5 * numpy.diff(x) * (y[1:] + y[:-1]))
        

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

def getData(fn='/home/truter/lab/datasets/xld_130902.spec'):
    analysis = sf.specfile.Specfile(fn)
    xmcdArr = []
    xasArr  = []
    avgA, avgB = [], []
    spec = analysis[3]
    x = spec[0][:]
    avgA = spec[1][:]
    avgB = spec[2][:]
    xmcdArr = spec[3][:]
    xasArr  = spec[4][:]
    return x, avgA, avgB, xmcdArr, xasArr

if __name__ == '__main__':
#    from 
    
    # Dummt data
    
    app = qt.QApplication([])
    win = SumRulesWindow()
    x, avgA, avgB, xmcd, xas = getData()
    win.plotWindow.newCurve(x,xmcd, legend='xmcd', xlabel='ene_st', ylabel='zratio', info={}, replot=False, replace=False)
    win.plotWindow.newCurve(x,xas, legend='xas', xlabel='ene_st', ylabel='zratio', info={}, replot=False, replace=False)
#    win.addMarker(spinbox=win.spinboxDict.values()[0])
    win.layout().addWidget(win.addMarkerBox(label='Foobar'))
    win.show()
    app.exec_()
