import Object3DQt as qt

ORIGIN     = 0
BBOXMIN    = 1
BBOXCENTER = 2
BBOXMAX    = 3

class Object3DAnchorWidget(qt.QGroupBox):
    def __init__(self, parent = None, anchor = None):
        qt.QGroupBox.__init__(self, parent)
        self.setTitle('Anchor')
        text  = 'Specify what origin to take for 3D object\n'
        text += 'translations and rotations.\n'
        text += 'It can be the origin of coordinates\n'
        text += 'or relative to the object bounding box'
        self.setToolTip(text)
        self.build()
        if anchor is None:
            anchor = [BBOXCENTER, BBOXCENTER, BBOXCENTER]
        self.setAnchor(anchor)
        

    def build(self):
        self.l = qt.QGridLayout(self)
        self.l.setMargin(4)
        self.l.setSpacing(4)
        i = 0
        j = 0
        for header in ['Axis', 'Origin', 'Min', 'Center', 'Max']:
            label = qt.QLabel(header)
            self.l.addWidget(label, i, j)
            if j != 0:
                self.l.setAlignment(label, qt.Qt.AlignHCenter)
            j += 1
        i = 1
        self.buttonGroupList = []
        for axis in ['X', 'Y', 'Z']:
            label = qt.QLabel(axis)
            self.l.addWidget(label, i, 0)
            j = 1
            buttonGroup = qt.QButtonGroup(self)
            for position in ['ORIGIN', 'MIN', 'CENTER', 'MAX']:
                rButton = qt.QRadioButton(self)
                self.l.addWidget(rButton, i, j)
                self.l.setAlignment(rButton, qt.Qt.AlignHCenter)
                buttonGroup.addButton(rButton)
                buttonGroup.setId(rButton, j - 1)
                j += 1
            self.buttonGroupList.append(buttonGroup)
            self.connect(buttonGroup,
                         qt.SIGNAL('buttonPressed(QAbstractButton *)'),
                         self._slot)
            i += 1

        """
        for header in ['Translation']:
            label = qt.QLabel(header)
            self.l.addWidget(label, 0, j+1)
            self.l.setAlignment(label, qt.Qt.AlignHCenter)

        i = 1
        self.spinBoxList = []
        for axis in ['X', 'Y', 'Z']:
            label = qt.QLabel(axis)
            self.l.addWidget(label, i, j)
            slider = qt.QDoubleSpinBox(self)
            slider.setRange(-1000., 1000.)
            self.l.addWidget(slider, i, j + 1)
            self.l.setAlignment(slider, qt.Qt.AlignVCenter)
            self.spinBoxList.append(slider)
            #self.connect(slider,
            #             qt.SIGNAL('valueChanged(double)'),
            #             self._slot)
            i += 1
        """
        self.l.setRowStretch(0, 0)
        self.l.setRowStretch(1, 1)
        self.l.setRowStretch(2, 1)
        self.l.setRowStretch(3, 1)


    def _slot(self, button):
        button.setChecked(True)
        self._signal()

    def _signal(self, event = None):
        if event is None:
            event = 'AnchorUpdated'
        anchor = self.getAnchor()
        ddict = {}
        ddict['event']  = 'AnchorUpdated'
        ddict['anchor'] = anchor 
        self.emit(qt.SIGNAL('Object3DAnchorSignal'), ddict)

    def setAnchor(self, anchor):
        for i in range(3):
            self.buttonGroupList[i].button(anchor[i]).setChecked(True)
    
    def getAnchor(self):
        anchorList = [0, 0, 0]
        for i in range(3):
            n = self.buttonGroupList[i].checkedId()
            if n >= 0:
                anchorList[i] = n
            else:
                print "WARNING: getAnchor -> Unselected button"
        return anchorList

class Object3DTranslationWidget(qt.QGroupBox):
    def __init__(self, parent = None, translation = None, labels = None):
        qt.QGroupBox.__init__(self, parent)
        text  = 'Specify translations of 3D objects\n'
        text += 'in the OpenGL scene.\n'
        self.setToolTip(text)
        self.setTitle('Translation')
        self.build(labels)
        if translation is None:
            translation = [0.0, 0.0, 0.0]
        self.setTranslation(translation)
        

    def build(self, labels):
        self.l = qt.QGridLayout(self)
        self.l.setMargin(4)
        self.l.setSpacing(4)
        if labels is None:labels = ['Axis','Amount']
        i = 0
        if 1:
            j = 0
            for header in labels:
                label = qt.QLabel(header)
                self.l.addWidget(label, i, j)
                #self.l.setAlignment(label, qt.Qt.AlignHCenter)
                j += 1
            i += 1
        self.spinBoxList = []
        for axis in ['X', 'Y', 'Z']:
            label = qt.QLabel(axis)
            self.l.addWidget(label, i, 0)
            slider = qt.QDoubleSpinBox(self)
            slider.setRange(-10000., 10000.)
            self.l.addWidget(slider, i, 1)
            self.l.setAlignment(slider, qt.Qt.AlignHCenter)
            self.spinBoxList.append(slider)
            self.connect(slider,
                         qt.SIGNAL('valueChanged(double)'),
                         self._slot)
            i += 1
        self.l.setRowStretch(0, 0)
        self.l.setRowStretch(1, 1)
        self.l.setRowStretch(2, 1)
        self.l.setRowStretch(3, 1)
        self.l.setColumnStretch(0, 0)

    def _slot(self, value):
        self._signal()

    def _signal(self, event = None):
        if event is None:
            event = 'TranslationUpdated'
        translation = self.getTranslation()
        ddict = {}
        ddict['event']  = event
        ddict['translation'] = translation 
        self.emit(qt.SIGNAL('Object3DTranslationSignal'), ddict)

    def setTranslation(self, translation):
        for i in range(3):
            self.spinBoxList[i].setValue(translation[i])
    
    def getTranslation(self):
        translationList = [0, 0, 0]
        for i in range(3):
            translationList[i] = self.spinBoxList[i].value()
        return translationList


class Object3DRotationWidget(qt.QGroupBox):
    def __init__(self, parent = None, rotation = None):
        qt.QGroupBox.__init__(self, parent)
        self.setTitle('Rotation')
        text  = 'Specify rotations of 3D objects\n'
        text += 'in the OpenGL scene.\n'
        self.setToolTip(text)
        self.build()
        if rotation is None:
            rotation = [0.0, 0.0, 0.0]
        self.setRotation(rotation)
        

    def build(self):
        self.l = qt.QGridLayout(self)
        self.l.setMargin(0)
        self.l.setSpacing(0)
        i = 0
        j = 0
        for header in ['X Rot','Y Rot', 'Z Rot']:
            label = qt.QLabel(header)
            self.l.addWidget(label, i, j)
            self.l.setAlignment(label, qt.Qt.AlignHCenter)
            j += 1
        i = 1
        j = 0
        self.spinList = []
        self.lineEditList = []
        self.validatorList = []
        for axis in ['X', 'Y', 'Z']:
            slider = qt.QDial(self)
            slider.setRange(0., 36000.)
            slider.setWrapping(1)
            self.l.addWidget(slider, i, j)
            self.l.setAlignment(slider, qt.Qt.AlignHCenter)
            self.spinList.append(slider)
            self.connect(slider,
                         qt.SIGNAL('valueChanged(int)'),
                         self._slot)
            lineEdit = qt.QLineEdit(self)
            v = qt.QDoubleValidator(lineEdit)
            lineEdit.setValidator(v)
            self.validatorList.append(v)
            self.l.addWidget(lineEdit, i+1, j)
            self.lineEditList.append(lineEdit)
            self.connect(lineEdit,
                         qt.SIGNAL('editingFinished()'),
                         self._lineSlot)
            j += 1
        """
        self.l.setRowStretch(0, 0)
        self.l.setRowStretch(1, 1)
        self.l.setRowStretch(2, 1)
        self.l.setRowStretch(3, 1)
        self.l.setColumnStretch(0, 0)
        """
        
    def _slot(self, value):
        self._signal(emit = True)

    def _lineSlot(self):
        rotation = [0.0, 0.0, 0.0]
        for i in range(3):
            value = float(str(self.lineEditList[i].text()))
            rotation[i]=value % 360
        self.setRotation(rotation)
        self._signal(emit = True)


    def _signal(self, event = None, emit = None):
        if emit is None : emit = False
        if not emit:return
        if event is None:
            event = 'RotationUpdated'
        rotation = self.getRotation()
        for i in range(3):
            self.lineEditList[i].setText("%.2f" % rotation[i])
        ddict = {}
        ddict['event']  = event
        ddict['rotation'] = rotation 
        self.emit(qt.SIGNAL('Object3DRotationSignal'), ddict)

    def setRotation(self, rotation):
        for i in range(3):
            rotation[i] = rotation[i] % 360
            value = (int(rotation[i] * 100) + 18000) % 36000
            self.spinList[i].setValue(value)
            self.lineEditList[i].setText("%f" % rotation[i])

    
    def getRotation(self):
        rotationList = [0, 0, 0]
        for i in range(3):
            rotationList[i] = (self.spinList[i].value() - 18000)/100.0
            if rotationList[i] < 0:
                rotationList[i] = rotationList[i] + 360.
        return rotationList

class Object3DMovement(qt.QWidget):
    def __init__(self, parent = None,
                       anchor = None,
                       translation = None,
                       rotation = None,
                       connect = True):
        qt.QWidget.__init__(self, parent)
        self.l = qt.QHBoxLayout(self)
        self.l.setMargin(0)
        self.l.setSpacing(0)
        self.anchorWidget      = Object3DAnchorWidget(self, anchor)
        self.translationWidget = Object3DTranslationWidget(self, translation)
        self.rotationWidget    = Object3DRotationWidget(self, rotation)

        self.setAnchor = self.anchorWidget.setAnchor
        self.getAnchor = self.anchorWidget.getAnchor
        self.setTranslation = self.translationWidget.setTranslation
        self.getTranslation = self.translationWidget.getTranslation
        self.setRotation = self.rotationWidget.setRotation
        self.getRotation = self.rotationWidget.getRotation

        self.l.addWidget(self.anchorWidget)
        self.l.addWidget(self.translationWidget)
        self.l.addWidget(self.rotationWidget)

        if connect:
            self.connect(self.anchorWidget,
                     qt.SIGNAL('Object3DAnchorSignal'),
                     self._anchorSlot)
            self.connect(self.translationWidget,
                     qt.SIGNAL('Object3DTranslationSignal'),
                     self._translationSlot)
            self.connect(self.rotationWidget,
                     qt.SIGNAL('Object3DRotationSignal'),
                     self._rotationSlot)

    def _anchorSlot(self, ddict):
        self._emitSignal()

    def _translationSlot(self, ddict):
        self._emitSignal()

    def _rotationSlot(self, ddict):
        self._emitSignal()

    def _emitSignal(self):
        ddict = self.getParameters()
        ddict['event']  = 'Object3DMovementUpdated'
        self.emit(qt.SIGNAL('Object3DMovementSignal'), ddict)

    def getParameters(self):
        ddict= {}
        ddict['anchor'] = self.anchorWidget.getAnchor()
        ddict['translation'] = self.translationWidget.getTranslation()
        ddict['rotation']    = self.rotationWidget.getRotation()
        return ddict

    def setParameters(self, ddict):
        translation = ddict.get('translation', [0.0, 0.0, 0.0])
        rotation    = ddict.get('rotation',    [0.0, 0.0, 0.0])
        anchor      = ddict.get('anchor',      [0 , 0 , 0])
        self.anchorWidget.setAnchor(anchor)
        self.rotationWidget.setRotation(rotation)
        self.translationWidget.setTranslation(translation)

if __name__ == "__main__":
    import sys
    app = qt.QApplication(sys.argv)
    def myslot(ddict):
        print "Signal received"
        print "Anchor      = ", ddict['anchor']
        print "Translation = ", ddict['translation']
        print "Rotation    = ", ddict['rotation']

    w = Object3DMovement()
    qt.QObject.connect(w,
                       qt.SIGNAL('Object3DMovementSignal'),
                       myslot)
    w.show()    
    app.exec_()
