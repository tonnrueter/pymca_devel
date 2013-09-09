import numpy

try:
    from PyMca import Plugin1DBase
except ImportError:
    from . import Plugin1DBase

try:
    from PyMca import MotorInfoWindow
except ImportError:
    print("MotorInfoPlugin importing from somewhere else")
    import MotorInfoWindow
    
DEBUG = 0
class MotorInfo(Plugin1DBase.Plugin1DBase):
    def __init__(self,  plotWindow,  **kw):
        Plugin1DBase.Plugin1DBase.__init__(self,  plotWindow,  **kw)
        self.methodDict = {}
        text = 'Show values of different motors.'
        function = self.showMotorInfo
        icon = None
        info = text
        self.methodDict["Show Motor Info"] =[function, info, icon]
        self.widget = None
    
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

    def showMotorInfo(self):
        curves = self.getAllCurves()
        nCurves = len(curves)
        if DEBUG:
            print ("Received %d curve(s).." % nCurves)
        if nCurves:
            self.legendsList = [leg for (xvals, yvals,  leg,  info) in curves] 
            infoList = [info for (xvals, yvals,  leg,  info) in curves] 
            self.motorValuesList = self._convertInfoDictionary( infoList )
        else:
            raise ValueError("No Curve(s) present.")
        self._createWidget()
        ret = self.widget.show() # non-modal QDialog

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
        self.widget = MotorInfoWindow.MotorInfoWidget(parent,  self.legendsList,  self.motorValuesList)

MENU_TEXT = "Motor Info"
def getPlugin1DInstance(plotWindow,  **kw):
    ob = MotorInfo(plotWindow)
    return ob

if __name__ == "__main__":
    # Basic test setup
    from PyMca import Plot1D
    x = numpy.arange(100.)
    y = x * x
    plot = Plot1D.Plot1D()
    plot.addCurve(x, y, "Curve1", {'MotorNames': "foo bar",  'MotorValues': "3.14 2.97"})
    plot.addCurve(x+100, -x*x, "Curve2", {'MotorNames': "baz",  'MotorValues': "6.28"})
    plugin = getPlugin1DInstance(plot)
