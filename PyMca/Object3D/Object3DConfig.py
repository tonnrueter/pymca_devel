import Object3DQt as qt
import Object3DMovement
import Object3DProperties
import ClippingPlaneConfiguration
import Object3DColormap
from VerticalSpacer import VerticalSpacer

class Object3DConfig(qt.QWidget):
    def __init__(self, parent = None):
        qt.QWidget.__init__(self, parent)
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setMargin(4)
        self.mainLayout.setSpacing(4)

        #drawing
        self.mainTab = qt.QTabWidget(self)
        self.tabDrawing = qt.QWidget(self.mainTab)
        self.tabDrawing.mainLayout = qt.QVBoxLayout(self.tabDrawing) 
        self.tabDrawing.mainLayout.setMargin(0)
        self.tabDrawing.mainLayout.setSpacing(0)
        self.movementsWidget  = Object3DMovement.Object3DMovement(self.tabDrawing)
        self.scaleWidget      = Object3DProperties.Object3DScale(self.tabDrawing)
        self.propertiesWidget = Object3DProperties.Object3DProperties(self.tabDrawing)
        self.tabDrawing.mainLayout.addWidget(self.movementsWidget)
        self.tabDrawing.mainLayout.addWidget(self.scaleWidget)
        self.tabDrawing.mainLayout.addWidget(self.propertiesWidget)
        self.mainTab.addTab(self.tabDrawing, "DRAWING")

        #clipping
        self.tabClipping = qt.QWidget(self.mainTab)
        self.tabClipping.mainLayout = qt.QVBoxLayout(self.tabClipping)
        self.tabClipping.mainLayout.setMargin(0)
        self.tabClipping.mainLayout.setSpacing(0)
        self.clippingPlaneWidget = ClippingPlaneConfiguration.ClippingPlaneWidget(self.tabClipping)
        self.colormapWidget = Object3DColormap.Object3DColormap(self.tabClipping)
        self.tabClipping.mainLayout.addWidget(self.clippingPlaneWidget)
        self.tabClipping.mainLayout.addWidget(self.colormapWidget)
        self.tabClipping.mainLayout.addWidget(VerticalSpacer())
        self.mainTab.addTab(self.tabClipping, "CLIP and COLOR")

        self.mainLayout.addWidget(self.mainTab)

        self.connect(self.movementsWidget,
                     qt.SIGNAL('Object3DMovementSignal'),
                     self._movementsSlot)
        self.connect(self.scaleWidget,
                     qt.SIGNAL('Object3DScaleSignal'),
                     self._scaleSlot)
        self.connect(self.propertiesWidget,
                     qt.SIGNAL('Object3DPropertiesSignal'),
                     self._propertiesSlot)
        self.connect(self.clippingPlaneWidget,
                     qt.SIGNAL('ClippingPlaneWidgetSignal'),
                     self._clippingPlaneSlot)

        self.connect(self.colormapWidget,
                     qt.SIGNAL('Object3DColormapSignal'),
                     self._colormapSlot)

    def _movementsSlot(self, ddict0):
        ddict = {}
        ddict['common'] = ddict0
        ddict['common'].update(self.propertiesWidget.getParameters())
        self._signal(ddict)

    def _scaleSlot(self, ddict0):
        ddict = {}
        event = ddict0['event']
        ddict['common'] = ddict0
        ddict['common'].update(self.scaleWidget.getParameters())

        #get the movement positions
        movementsDict = self.movementsWidget.getParameters()
        ddict['common'].update(movementsDict)
        if event == "xScaleUpdated":
            if ddict['common']['scale'][0] != 0.0:
                ddict['common']['translation'][0] /= ddict0['magnification']
        if event == "yScaleUpdated":
            if ddict['common']['scale'][1] != 0.0:
                ddict['common']['translation'][1] /= ddict0['magnification']
        if event == "zScaleUpdated":
            if ddict['common']['scale'][2] != 0.0:
                ddict['common']['translation'][2] /= ddict0['magnification']
        self.movementsWidget.setParameters(ddict['common'])
        self._signal(ddict)

    def _propertiesSlot(self, ddict):
        ddict['common'].update(self.movementsWidget.getParameters()) 
        ddict['common'].update(self.scaleWidget.getParameters()) 
        self._signal(ddict)

    def _clippingPlaneSlot(self, ddict0):
        ddict = {}
        ddict['common'] = ddict0
        self._signal(ddict)

    def _colormapSlot(self, ddict0):
        ddict = {}
        ddict['common'] = ddict0
        self._signal(ddict)

    def _signal(self, ddict):
        self.emit(qt.SIGNAL('Object3DConfigSignal'), ddict)

    def getConfiguration(self):
        ddict = self.propertiesWidget.getParameters()
        ddict['common'].update(self.movementsWidget.getParameters())
        ddict['common'].update(self.scaleWidget.getParameters())
        ddict['common'].update(self.clippingPlaneWidget.getParameters())
        ddict['common'].update(self.colormapWidget.getParameters())
        return ddict

    def setConfiguration(self, ddict):
        self.movementsWidget.setParameters(ddict['common'])
        self.scaleWidget.setParameters(ddict['common'])
        self.propertiesWidget.setParameters(ddict)
        self.clippingPlaneWidget.setParameters(ddict['common'])
        self.colormapWidget.setParameters(ddict['common'])
 
if __name__ == "__main__":
    import sys
    app = qt.QApplication(sys.argv)
    def myslot(ddict):
        print "Signal received"
        print "dict = ", ddict

    w = Object3DConfig()
    qt.QObject.connect(w,
                       qt.SIGNAL('Object3DConfigSignal'),
                       myslot)
    w.show()    
    app.exec_()
