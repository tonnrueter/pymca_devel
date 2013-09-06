import sys
import Object3DQt as qt
from HorizontalSpacer import HorizontalSpacer
from VerticalSpacer import VerticalSpacer
import numpy
DEBUG = 0

class SceneCoordinates(qt.QWidget):
    def __init__(self, parent = None):
        qt.QWidget.__init__(self, parent)
        self.mainLayout = qt.QGridLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(0)
        self.limitsWidget = SceneLimitsWidget(self)
        self.observerWidget = ObserverPositionWidget(self)
        self.axesVectorsWidget = SceneAxesVectorsWidget(self)
        #self.viewOrientationWidget = ViewOrientationWidget(self)
        self.mainLayout.addWidget(self.limitsWidget, 0, 0)
        self.mainLayout.addWidget(self.observerWidget, 0, 1)
        self.mainLayout.addWidget(self.axesVectorsWidget, 1, 0)
        self.mainLayout.addWidget(HorizontalSpacer(self), 1, 2, 1, 2)
        #self.mainLayout.addWidget(self.viewOrientationWidget, 1, 0)
        #self.mainLayout.addWidget(HorizontalSpacer(self), 0, 1, 1, 3)
        #self.mainLayout.addWidget(VerticalSpacer(self), 1, 0, 1, 3)
        self.mainLayout.setColumnStretch(0, 0)
        self.mainLayout.setColumnStretch(1, 0)
        self.mainLayout.setColumnStretch(2, 1)
        #self.mainLayout.setRowStretch(0, 0)
        #self.mainLayout.setRowStretch(1, 1)

        #connections
        self.connect(self.limitsWidget,
                     qt.SIGNAL('SceneLimitsSignal'),
                     self._emitSignal)
        self.connect(self.observerWidget,
                     qt.SIGNAL('ObserverPositionSignal'),
                     self._emitSignal)
        self.connect(self.axesVectorsWidget,
                     qt.SIGNAL('SceneAxesVectorsSignal'),
                     self._emitSignal)

    def getParameters(self):
        ddict = self.limitsWidget.getParameters()
        ddict.update(self.observerWidget.getParameters())
        ddict.update(self.axesVectorsWidget.getParameters())
        return ddict

    def setParameters(self, ddict):
        self.limitsWidget.setParameters(ddict)
        self.observerWidget.setParameters(ddict)
        self.axesVectorsWidget.setParameters(ddict)

    def _emitSignal(self, ddict=None):
        if ddict is None:
            event = 'SceneCoordinatesChanged'
        else:
            event = ddict['event']
        ddict = self.getParameters()
        ddict['event'] = event
        self.emit(qt.SIGNAL('SceneCoordinatesSignal'), ddict)

class SceneAxesVectorsWidget(qt.QGroupBox):
    def __init__(self, parent = None):
        qt.QGroupBox.__init__(self, parent)
        self.setTitle('Axes Vectors')
        self.mainLayout = qt.QGridLayout(self)
        self.mainLayout.setMargin(4)
        self.mainLayout.setSpacing(4)
        self.build()
        ddict = {'u':[1.0, 0.0, 0.0],
                 'v':[0.0, 1.0, 0.0],
                 'w':[0.0, 0.0, 1.0],
                 'uvw':False}
        self.setParameters(ddict)

        
    def build(self):
        self._entryDict = {'u':[],'v':[],'w':[]}
        alignment = qt.Qt.AlignCenter
        i = 0
        c = 0
        for t in ["i", "j", "k"]:
            c += 1
            l = qt.QLabel(self)
            l.setAlignment(alignment)
            l.setText("%s" % t)
            self.mainLayout.addWidget(l, i, c)
            
        i = 0
        for t in ["u", "v", "w"]:
            i += 1
            l = qt.QLabel(self)
            l.setText("%s" % t)
            self.mainLayout.addWidget(l, i, 0)
            for c in range(3):
                l = qt.QLineEdit(self)
                l._v = qt.QDoubleValidator(l)
                l.setValidator(l._v)
                self._entryDict[t].append(l)
                self.mainLayout.addWidget(l, i, c + 1)

        self.buildActions()

    def buildActions(self):
        self._useCheckBox = qt.QCheckBox(self)
        self._useCheckBox.setText('Use Transformation')
        self._useCheckBox.setChecked(True)
        self.mainLayout.addWidget(self._useCheckBox, 4, 0,
                                                     1, 2,
                                                     qt.Qt.AlignCenter)

        self._updateButton = qt.QPushButton(self)
        self._updateButton.setText('Update')
        self.mainLayout.addWidget(self._updateButton, 4, 2,
                                                     1, 2,
                                                     qt.Qt.AlignCenter)

        
        self.connect(self._useCheckBox,
                     qt.SIGNAL('clicked()'),
                     self._emitSignal)

        self.connect(self._updateButton,
                     qt.SIGNAL('clicked()'),
                     self._emitSignal)

    def _emitSignal(self):
        ddict = self.getParameters()
        ddict['event'] = 'SceneAxesVectorsChanged'               
        self.emit(qt.SIGNAL('SceneAxesVectorsSignal'), ddict)

    def getParameters(self):
        ddict = {}
        for t in ["u", "v", "w"]:
            ddict[t] = []
            entryList = self._entryDict[t]
            for i in range(3):
                ddict[t].append(float(entryList[i].text()))
        if self._useCheckBox.isChecked():
            ddict['uvw'] = True
        else:
            ddict['uvw'] = False
        return ddict

    def setParameters(self, ddict):
        for t in ["u", "v", "w"]:
            if t in ddict:
                if len(ddict[t]) == 3:
                    entryList = self._entryDict[t]
                    for i in range(3):
                       entryList[i].setText("%f" % ddict[t][i])
        if 'uvw' in ddict:
            if ddict['uvw']:
                self._useCheckBox.setChecked(True)
            else:
                self._useCheckBox.setChecked(False)
        return


class SceneLimitsWidget(qt.QGroupBox):
    def __init__(self, parent = None):
        qt.QGroupBox.__init__(self, parent)
        self.setTitle('Visual Volume')
        self.mainLayout = qt.QGridLayout(self)
        self.mainLayout.setMargin(4)
        self.mainLayout.setSpacing(4)
        self.build()
        ddict = {'xmin': -10.0,
                 'xmax':  10.0,
                 'ymin': -10.0,
                 'ymax':  10.0,
                 'zmin': -1.0,
                 'zmax':  1.0,
                 'autoscale':True}
        self.setParameters(ddict)

    def build(self):
        self._entryList = []
        self._centerList = []
        self._deltaList  = []
        i = 0
        for t in ['X', 'Y', 'Z']:
            l = qt.QLabel(self)
            l.setText("%s Min:" % t)
            self.mainLayout.addWidget(l, i, 0)

            l = qt.QLineEdit(self)
            l._v = qt.QDoubleValidator(l)
            l.setValidator(l._v)
            self._entryList.append(l)
            self.mainLayout.addWidget(l, i, 1)
            self.connect(l,
                         qt.SIGNAL('editingFinished()'),
                         self._emitSignal)
            
            l = qt.QLabel(self)
            l.setText("%s Max:" % t)
            self.mainLayout.addWidget(l, i, 2)

            l = qt.QLineEdit(self)
            l._v = qt.QDoubleValidator(l)
            l.setValidator(l._v)
            self._entryList.append(l)
            self.mainLayout.addWidget(l, i, 3)
            self.connect(l,
                         qt.SIGNAL('editingFinished()'),
                         self._emitSignal)

            l = qt.QLabel(self)
            l.setText("%s Center:" % t)
            self.mainLayout.addWidget(l, i, 4)

            l = qt.QLineEdit(self)
            l._v = qt.QDoubleValidator(l)
            l.setValidator(l._v)
            self._centerList.append(l)
            self.mainLayout.addWidget(l, i, 5)

            if i == 0:  
                self.connect(l,
                         qt.SIGNAL('editingFinished()'),
                         self.__xCenterChanged)
            elif i == 1:
                self.connect(l,
                         qt.SIGNAL('editingFinished()'),
                         self.__yCenterChanged)
            else:
                self.connect(l,
                         qt.SIGNAL('editingFinished()'),
                         self.__zCenterChanged)

            l = qt.QLabel(self)
            l.setText("%s Delta:" % t)
            self.mainLayout.addWidget(l, i, 6)

            l = qt.QLineEdit(self)
            l._v = qt.QDoubleValidator(l)
            l.setValidator(l._v)
            self._deltaList.append(l)
            if i == 0:  
                self.connect(l,
                         qt.SIGNAL('editingFinished()'),
                         self.__xDeltaChanged)
            elif i == 1:
                self.connect(l,
                         qt.SIGNAL('editingFinished()'),
                         self.__yDeltaChanged)
            else:
                self.connect(l,
                         qt.SIGNAL('editingFinished()'),
                         self.__zDeltaChanged)

            self.mainLayout.addWidget(l, i, 7)
            i += 1
            
        self._autoCheckBox = qt.QCheckBox(self)
        self._autoCheckBox.setText('Automatic')
        self._autoCheckBox.setChecked(True)
        self.mainLayout.addWidget(self._autoCheckBox, i, 0,
                                                      1, 4,
                                                      qt.Qt.AlignCenter)
        self.connect(self._autoCheckBox,
                     qt.SIGNAL('clicked()'),
                     self._autoCheckBoxClicked)
        return
        self._applyButton = qt.QPushButton(self)
        self._applyButton.setText('Apply')
        self._applyButton.setAutoDefault(False)
        self.mainLayout.addWidget(self._applyButton, i, 3)
        
    def getParameters(self):
        ddict = {'xmin': 0.0,
                 'xmax': 0.0,
                 'ymin': 0.0,
                 'ymax': 0.0,
                 'zmin': 0.0,
                 'zmax': 0.0,
                 'autoscale':True}
        labelList = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']
        i = 0
        for key in labelList:
            try:
                ddict[key] = float(self._entryList[i].text())
            except:                
                pass
            i += 1
        if self._autoCheckBox.isChecked():
            ddict['autoscale'] = True
        else:
            ddict['autoscale'] = False
        for i in range(3):
            self.__updateCenterAndDelta(i*2)
        return ddict

    def setParameters(self, ddict):
        labelList = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']
        i = 0
        updateList = []
        for key in labelList:
            if ddict.has_key(key):
                self._entryList[i].setText("%g" % ddict[key])
                if int(i/2) not in updateList:
                    updateList.append(int(i/2))
            i += 1
        for i in updateList:
            self.__updateCenterAndDelta(2*i)
        flag = ddict.get('autoscale', True)
        if flag:
            self._autoCheckBox.setChecked(True)
        else:
            self._autoCheckBox.setChecked(False)
        self._updateEntries()

    def _autoCheckBoxClicked(self):
        self._updateEntries()
        self._emitSignal()

    def _updateEntries(self):
        if self._autoCheckBox.isChecked():
            for e in self._entryList:
                e.setDisabled(True)
            for e in self._centerList:
                e.setDisabled(True)
            for e in self._deltaList:
                e.setDisabled(True)
        else:
            for e in self._entryList:
                e.setDisabled(False)
            for e in self._centerList:
                e.setDisabled(False)
            for e in self._deltaList:
                e.setDisabled(False)
            

    def _emitSignal(self):
        ddict = self.getParameters()
        ddict['event'] = 'SceneLimitsChanged'               
        self.emit(qt.SIGNAL('SceneLimitsSignal'), ddict)

    def __xCenterChanged(self):
        self.__updateFromCenterAndDelta(0)

    def __yCenterChanged(self):
        self.__updateFromCenterAndDelta(1)

    def __zCenterChanged(self):
        self.__updateFromCenterAndDelta(2)

    def __xDeltaChanged(self):
        self.__updateFromCenterAndDelta(0)

    def __yDeltaChanged(self):
        self.__updateFromCenterAndDelta(1)

    def __zDeltaChanged(self):
        self.__updateFromCenterAndDelta(2)

    def __updateFromCenterAndDelta(self, idx):
        center = float(self._centerList[idx].text())
        delta = float(self._deltaList[idx].text())
        self._entryList[2*idx].setText("%g" % (center-0.5*delta))
        self._entryList[2*idx+1].setText("%g" % (center+0.5*delta))
        self._emitSignal()

    def __updateCenterAndDelta(self, idx):
        i = int(idx/2)
        try:
            vmin = float(self._entryList[2*i].text())
            vmax = float(self._entryList[2*i+1].text())
            self._centerList[i].setText("%g" % (0.5*(vmin+vmax)))
            self._deltaList[i].setText("%g" % (vmax-vmin))
        except:
            #This can happen if there is no text.
            pass

class ViewOrientationWidget(qt.QGroupBox):
    def __init__(self, parent = None):
        qt.QGroupBox.__init__(self, parent)
        self.setTitle('View Orientation')
        self.mainLayout = qt.QGridLayout(self)
        self.mainLayout.setMargin(4)
        self.mainLayout.setSpacing(4)
        self.build()
        viewMatrix = numpy.zeros((4,4), numpy.float32)
        for i in [0, 1, 2, 3]:
            viewMatrix[i, i] = 1.0
        ddict = {'theta': 0.0,
                 'phi':0.0,
                 'view':viewMatrix} 
        self.setParameters(ddict)

    def build(self):
        self._entryList = []
        for i in range(4):
            for j in range(4):
                l = qt.QLineEdit(self)
                l.setReadOnly(True)
                l._v = qt.QDoubleValidator(l)
                l.setValidator(l._v)
                self._entryList.append(l)
                self.mainLayout.addWidget(l, i, j)

        # Theta
        i = 4
        l = qt.QLabel(self)
        l.setText("Theta:")
        self.mainLayout.addWidget(l, i, 0)
        l = qt.QLineEdit(self)
        l.setReadOnly(True)
        l._v = qt.QDoubleValidator(l)
        l.setValidator(l._v)
        self._theta = l
        self.mainLayout.addWidget(l, i, 1)

        #Phi
        l = qt.QLabel(self)
        l.setText("Phi:")
        self.mainLayout.addWidget(l, i, 2)
        l = qt.QLineEdit(self)
        l.setReadOnly(True)
        l._v = qt.QDoubleValidator(l)
        l.setValidator(l._v)
        self._phi = l
        self.mainLayout.addWidget(l, i, 3)        

    def setParameters(self, ddict):
        if ddict.has_key('view'):
            k = 0
            for i in range(4):
                for j in range(4):
                    self._entryList[k].setText("%g" % float(ddict['view'][i, j]))
                    k = k + 1

        if ddict.has_key('theta'):
            self._theta.setText("%.3g" % float(ddict['theta']))

        if ddict.has_key('phi'):
            self._phi.setText("%.3g" % float(ddict['phi']))

    def getParameters(self):
        m = numpy.zeros((4,4), numpy.float32)
        k = 0
        for i in range(4):
            for j in range(4):
                m[i,j] = float(self._entryList[k].text())
                k = k + 1
        theta = float(self._theta.text())
        phi = float(self._phi.text())
        ddict = {}
        ddict['view'] = m
        ddict['theta'] = theta
        ddict['phi'] = phi
  

class ObserverPositionWidget(qt.QGroupBox):
    def __init__(self, parent = None):
        qt.QGroupBox.__init__(self, parent)
        self.setTitle('Observer Position')
        self.mainLayout = qt.QGridLayout(self)
        self.mainLayout.setMargin(4)
        self.mainLayout.setSpacing(4)
        self.build()
        ddict = {'observer': [0.0, 0.0, 0.0],
                 'autolocate':True}
        self.setParameters(ddict)

    def build(self):
        self._entryList = []
        i = 0
        for t in ['X:', 'Y:', 'Z:']:
            l = qt.QLabel(self)
            l.setText(t)
            self.mainLayout.addWidget(l, i, 0)

            l = qt.QLineEdit(self)
            l._v = qt.QDoubleValidator(l)
            l.setValidator(l._v)
            self._entryList.append(l)
            self.mainLayout.addWidget(l, i, 1)
            self.connect(l,
                         qt.SIGNAL('editingFinished()'),
                         self._emitSignal)
            i += 1
        self._autoCheckBox = qt.QCheckBox(self)
        self._autoCheckBox.setText('Automatic')
        self._autoCheckBox.setChecked(True)
        self.mainLayout.addWidget(self._autoCheckBox, i, 0,
                                                      1, 2,
                                                      qt.Qt.AlignCenter)
        self.connect(self._autoCheckBox,
                     qt.SIGNAL('clicked()'),
                     self._emitSignal)


    def getParameters(self):
        position = [0.0, 0.0, 0.0]
        for i in range(3):
            try:
                position[i] = float(self._entryList[i].text())
            except:     
                pass
        ddict = {'observer': position}
        if self._autoCheckBox.isChecked():
            ddict['autolocate'] = True
        else:
            ddict['autolocate'] = False
        return ddict

    def setParameters(self, ddict):
        if not ddict.has_key('observer'):
            return
        position = ddict['observer']
        for i in range(3):
            self._entryList[i].setText("%f" % position[i])
        self._autoCheckBox.setChecked(ddict['autolocate'])

    def _emitSignal(self):
        ddict = self.getParameters()
        ddict['event'] = 'ObserverPositionChanged'
        self.emit(qt.SIGNAL('ObserverPositionSignal'), ddict)

if __name__ == "__main__":
    app = qt.QApplication([])
    w = SceneCoordinates()
    w.show()
    def slot(ddict):
        for key in ddict:
            print("Key = %s key Content = %s" % (key, ddict[key]))
    qt.QObject.connect(w, qt.SIGNAL('SceneCoordinatesSignal'),slot)
    app.exec_()
