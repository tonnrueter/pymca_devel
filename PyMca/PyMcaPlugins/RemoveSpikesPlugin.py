try:
    from PyMca import Plugin1DBase
except ImportError:
    from . import Plugin1DBase

try:
    from PyMca.PyMcaSciPy.signal import medfilt1d
except ImportError:
    print("RemoveSpikesPlugin importing directly")
    from PyMca.PyMcaSciPy.signal import medfilt1d

import numpy
    
DEBUG = True
class RemoveSpikes(Plugin1DBase.Plugin1DBase):
    def __init__(self,  plotWindow,  **kw):
        Plugin1DBase.Plugin1DBase.__init__(self,  plotWindow,  **kw)
        self.methodDict = {}
        text = 'Remove spikes from active curve'
        function = self.removeSpikesActive
        icon = None
        info = text
        self.methodDict["Remove spikes"] =[function, info, icon]
        text = 'Remove spikes from all curve'
        function = self.removeSpikesAll
        icon = None
        info = text
        self.methodDict["Remove spikes"] =[function, info, icon]
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

    def setThreshold(self):
        # TODO: Spawn QDialog, ask for threshold % window width
        pass

    def removeSpikesAll(self):
        self.medianThresholdFilter(False)

    def removeSpikesActive(self):
        self.medianThresholdFilter(True)

    def medianThresholdFilter(self, activeOnly, threshold=None, length=9):
        if activeOnly:
            active = self._plotWindow.getActiveCurve()
            if not active:
                return
        else:
            spectra = self._plotWindow.getAllCurves()
        for (i,spec) in enumerate(spectra):
            x, y, legend, info = spec
            filtered = medfilt1d(y, length)
            diff = abs(filtered-y)
            if not threshold:
                threshold = .5 * diff.max() # OR: diff.mean()
            ynew = numpy.where(diff>threshold, filtered, y)
            legend = info.get('selectionlegend','') + ' SR'
            if i==0 and (len(spectra)!=1):
                self._plotWindow.addCurve(x,ynew,legend,info, replace=True)
            else:
                self._plotWindow.addCurve(x,ynew,legend,info)
        

MENU_TEXT = "Remove spikes"
def getPlugin1DInstance(plotWindow,  **kw):
    ob = RemoveSpikes(plotWindow)
    return ob
    
if __name__ == "__main__":
    from PyMca import ScanWindow
    from PyMca import PyMcaQt as qt
    import numpy
    app = qt.QApplication([])
    
    sw = ScanWindow.ScanWindow()

    x = numpy.arange(1000.)
    y0 =  10 * x + 10000. * numpy.exp(-0.5*(x-500)*(x-500)/400) + 1500 * numpy.random.random(1000.)
    y1 =  10 * x + 10000. * numpy.exp(-0.5*(x-600)*(x-600)/400) + 1500 * numpy.random.random(1000.)
    y2 =  10 * x + 10000. * numpy.exp(-0.5*(x-400)*(x-400)/400) + 1500 * numpy.random.random(1000.)
    y2[320:322] = 50000.
    info0 = {'FileHeader':['#F /data/id08/inhouse/somerandomname'],'xlabel': 'foo', 'ylabel': 'arb', 'MotorNames': 'oxPS Motor11 Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3',  'MotorValues': '1 8.69271399699 21.9836418539 0.198068826612 0.484475455792 0.350252217264 0.663925270933 0.813033264421 0.221149410218 0.593188258866 0.678010392881 0.267389247833 0.677890617858'}
    info1 = {'MotorNames': 'PhaseD oxPS Motor16 Motor15 Motor14 Motor13 Motor12 Motor11 Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3', 'MotorValues': '0.470746882688 0.695816070299 0.825780811755 0.25876374531 0.739264467436 0.090842892619 2 0.213445659833 0.823400550314 0.020278096857 0.568744021322 0.85378115537 0.696730386891 0.269196313956 0.793293334395 0.769216567757 0.959092709527 0.0109264683697 0.538264972553'}
    info2 = {'MotorNames': 'PhaseD oxPS Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3',  'MotorValues': '2 0.44400576644 0.613870067852 0.901968648111 0.319768771085 0.571432278628 0.278675836163 0.154436774878 0.416231999332 0.294201017231 0.813913587748 0.577572903105 0.869045182568'}
    
    sw.addCurve(x, y0, legend="Curve0", info=info0, replot=False, replace=False)
    sw.addCurve(x, y1, legend="Curve1", info=info1, replot=False, replace=False)
    sw.addCurve(x, y2, legend="Curve2", info=info2, replot=False, replace=False)

    plugin = getPlugin1DInstance(sw)
    plugin.applyMethod(plugin.getMethods()[0])
    
    sw.show()
    
    app.exec_()
    
