import numpy
try:
    from PyMca import Plugin1DBase
except ImportError:
    from . import Plugin1DBase
    
try:
    from PyMca import SortPlotsWindow
except ImportError:
    print("SortPlotsWindow importing from somewhere else")
    import SortPlotsWindow
    
DEBUG = 0
class SortPlots(Plugin1DBase.Plugin1DBase):
    def __init__(self,  plotWindow,  **kw):
        Plugin1DBase.Plugin1DBase.__init__(self,  plotWindow,  **kw)
        self.methodDict = {}
        text = 'Sort plots for motor value.'
        function = self.showSortPlotsWindow
        icon = None
        info = text
        self.methodDict["Sort plots"] =[function, info, icon]
        self.widget = None
#        print( str(type(foo)) )
#        print( foo.__func__.func_name ) # Confirms
    
    def getMethods(self, plottype=None):
        names = list(self.methodDict.keys())
        names.sort()
        return names

    def getMethodToolTip(self, name):
        return self.methodDict[name][1]

    def getMethodPixmap(self, name):
        return self.methodDict[name][2]

    def applyMethod(self, name):
        self.methodDict[name][0]()
        return

    def addActions(self):
        pass

    def showSortPlotsWindow(self):
        self._setLists()
        if self.widget is None:
            self._createWidget()
        else:
            self._setLists()
            self.widget.updatePlots(self.legendsList,  self.motorValuesList)
        self.widget.show()
        self.widget.raise_()

    def _setLists(self):
        curves = self.getAllCurves()
        nCurves = len(curves)
        if DEBUG:
            print ("Received %d curve(s).." % nCurves)
        self.legendsList = [leg for (xvals, yvals,  leg,  info) in curves] 
        infoList = [info for (xvals, yvals,  leg,  info) in curves] 
        self.motorValuesList = self._convertInfoDictionary( infoList )

    def _convertInfoDictionary(self,  infosList):
        ret = []
        for info in infosList :
            motorNames = info.get('MotorNames',  None)
            if motorNames is not None:
                if type(motorNames) == str:
                    namesList = motorNames.split()
                elif type(motorNames) == list:
                    namesList = motorNames
                else:
                    namesList = []
            else:
                namesList = []
            motorValues = info.get('MotorValues',  None)
            if motorNames is not None:
                if type(motorValues) == str:
                    valuesList = motorValues.split()
                elif type(motorValues) == list:
                    valuesList = motorValues
                else:
                    valuesList = []
            else:
                valuesList = []
            if len(namesList) == len(valuesList):
                ret.append( dict( zip( namesList,  valuesList ) ) )
            else:
                print("Number of motors and values does not match!")
        return ret

    def _createWidget(self):
        parent = None
        self.widget = SortPlotsWindow.SortPlotsWidget(parent,  
                                                                                 list(self.legendsList),  
                                                                                 list(self.motorValuesList))
        self.widget.buttonUpdate.clicked.connect(self.showSortPlotsWindow)

MENU_TEXT = "Sort Plots"
def getPlugin1DInstance(plotWindow,  **kw):
    ob = SortPlots(plotWindow)
    return ob
    
if __name__ == "__main__":
    from PyMca import Plot1D
    x = numpy.arange(100.)
    y = x * x
    plot = Plot1D.Plot1D()
    plot.addCurve(x, y, "Curve1", {'MotorNames': "foo bar",  'MotorValues': "3.14 2.97"})
    plot.addCurve(x+100, -x*x, "Curve2", {'MotorNames': "baz",  'MotorValues': "6.28"})
    plugin = getPlugin1DInstance(plot)
    # Create app..
