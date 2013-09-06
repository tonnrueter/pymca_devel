#/*##########################################################################
# Copyright (C) 2004-2012 European Synchrotron Radiation Facility
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
import sys
import numpy
from PyMca import PyMcaQt as qt
SCANWINDOW = False
QWT = False
# This strange looking import is to workaround an endless import
try:
    from PyMca import ScanWindow
    SCANWINDOW = True
except:
    try:
        from PyMca import Plot1DQwt
        QWT = True
    except ImportError:
        from PyMca import Plot1DMatplotlib
from PyMca import SpecfitFuns


class StripParametersWidget(qt.QWidget):
    def __init__(self, parent=None):
        qt.QWidget.__init__(self, parent)
        self.build()

    def build(self):
        self.mainLayout = qt.QGridLayout(self)
        self.mainLayout.setMargin(11)
        self.mainLayout.setSpacing(6)

        #strip algorithm
        self.stripComboLabel = qt.QLabel(self)
        self.stripComboLabel.setText("Non-analytical (or estimation) background algorithm")
        self.stripCombo = qt.QComboBox(self)
        self.stripCombo.addItem(str("Strip"))
        self.stripCombo.addItem(str("SNIP"))
        self.connect(self.stripCombo,
                     qt.SIGNAL("activated(int)"),
                     self._stripComboActivated)

        #SNIP width
        self.snipWidthLabel = qt.QLabel(self)
        self.snipWidthLabel.setText(str("SNIP Background Width"))

        self.snipWidthSpin = qt.QSpinBox(self)
        self.snipWidthSpin.setMaximum(300)
        self.snipWidthSpin.setMinimum(0)
        self.connect(self.snipWidthSpin,
                     qt.SIGNAL('valueChanged(int)'),
                     self._emitSignal)

        #Strip width
        self.stripWidthLabel = qt.QLabel(self)
        self.stripWidthLabel.setText(str("Strip Background Width"))

        self.stripWidthSpin = qt.QSpinBox(self)
        self.stripWidthSpin.setMaximum(100)
        self.stripWidthSpin.setMinimum(1)
        self.connect(self.stripWidthSpin,
                     qt.SIGNAL('valueChanged(int)'),
                     self._emitSignal)

        #Strip iterations
        self.stripIterLabel = qt.QLabel(self)
        self.stripIterLabel.setText(str("Strip Background Iterations"))
        self.stripIterValue = qt.QLineEdit(self)
        validator = qt.QIntValidator(self.stripIterValue)
        self.stripIterValue._v = validator
        self.connect(self.stripIterValue,
                     qt.SIGNAL('editingFinished()'),
                     self._emitSignal)

        #Strip smoothing
        self.stripFilterLabel = qt.QLabel(self)
        self.stripFilterLabel.setText(str("Strip Background Smoothing Width (Savitsky-Golay)"))

        self.stripFilterSpin = qt.QSpinBox(self)
        self.stripFilterSpin.setMinimum(1)
        self.stripFilterSpin.setMaximum(40)
        self.stripFilterSpin.setSingleStep(2)
        self.connect(self.stripFilterSpin,
                     qt.SIGNAL('valueChanged(int)'),
                     self._emitSignal)

        #anchors
        self.anchorsContainer = qt.QWidget(self)
        anchorsContainerLayout = qt.QHBoxLayout(self.anchorsContainer)
        anchorsContainerLayout.setMargin(0)
        anchorsContainerLayout.setSpacing(2)
        self.stripAnchorsFlagCheck = qt.QCheckBox(self.anchorsContainer)
        self.stripAnchorsFlagCheck.setText(str("Strip Background use Anchors"))
        self.connect(self.stripAnchorsFlagCheck,
                     qt.SIGNAL('stateChanged(int)'),
                     self._emitSignal)
        anchorsContainerLayout.addWidget(self.stripAnchorsFlagCheck)

        #self.iterSpin = qt.QSpinBox(self)
        #self.iterSpin.setMinimum(1)

        maxnchannel  = 16384*4
        self.stripAnchorsList = []
        for i in range(4):
            anchorSpin = qt.QSpinBox(self.anchorsContainer)
            anchorSpin.setMinimum(0)
            anchorSpin.setMaximum(maxnchannel)
            self.connect(anchorSpin,
                     qt.SIGNAL('valueChanged(int)'),
                     self._emitSignal)
            anchorsContainerLayout.addWidget(anchorSpin)
            self.stripAnchorsList.append(anchorSpin)

        self.mainLayout.setColumnStretch(0, 1)
        row  = 0
        self.mainLayout.addWidget(self.stripComboLabel,  row, 0)
        self.mainLayout.addWidget(self.stripCombo, row, 4)

        row += 1
        self.mainLayout.addWidget(self.snipWidthLabel,row, 0)
        self.mainLayout.addWidget(self.snipWidthSpin, row, 4)

        row += 1
        self.mainLayout.addWidget(self.stripWidthLabel, row, 0)
        self.mainLayout.addWidget(self.stripWidthSpin,  row, 4)

        row += 1
        self.mainLayout.addWidget(self.stripIterLabel, row, 0)
        self.mainLayout.addWidget(self.stripIterValue, row, 4)
        
        row += 1
        self.mainLayout.addWidget(self.stripFilterLabel, row, 0)
        self.mainLayout.addWidget(self.stripFilterSpin,  row, 4)

        row += 1
        self.mainLayout.addWidget(self.anchorsContainer, row, 0, 1, 5)

        self._stripComboActivated(0)

    def _stripComboActivated(self, iValue):
        if iValue == 1:
            self.setSNIP(True)
        else:
            self.setSNIP(False)

    def setSNIP(self, bValue):
        if bValue:
            self.snipWidthSpin.setEnabled(True)
            self.stripWidthSpin.setEnabled(False)
            #self.stripFilterSpin.setEnabled(False)
            self.stripIterValue.setEnabled(False)
            self.stripCombo.setCurrentIndex(1)
        else:
            self.snipWidthSpin.setEnabled(False)
            #self.stripFilterSpin.setEnabled(True)
            self.stripWidthSpin.setEnabled(True)
            self.stripIterValue.setEnabled(True)
            self.stripCombo.setCurrentIndex(0)

    def setParameters(self, ddict):
        if 'fit' in ddict:
            pars = ddict['fit']
        else:
            pars = ddict

        key = "stripalgorithm"
        if key in pars:
            stripAlgorithm = int(pars[key])
            self.setSNIP(stripAlgorithm)
            
        key = "snipwidth"            
        if key in pars:
            self.snipWidthSpin.setValue(int(pars[key]))

        key = "stripwidth"            
        if key in pars:
            self.stripWidthSpin.setValue(int(pars[key]))

        key = "stripiterations"
        if key in pars:
            self.stripIterValue.setText("%d" % int(pars[key]))

        key = "stripfilterwidth"
        if key in pars:
            self.stripFilterSpin.setValue(int(pars[key]))

        key = "stripanchorsflag"
        if key in pars:
            self.stripAnchorsFlagCheck.setChecked(int(pars[key]))

        key = "stripanchorslist"
        if key in pars:
            anchorslist = pars[key]
            if anchorslist in [None, 'None']:
                anchorslist = []
            for spin in self.stripAnchorsList:
                spin.setValue(0)

            i = 0
            for value in anchorslist:
                self.stripAnchorsList[i].setValue(int(value))
                i += 1

    def getParameters(self):
        pars = {}
        pars["stripalgorithm"] = int(self.stripCombo.currentIndex())
        pars["stripconstant"]= 1.0
        pars["snipwidth"] = self.snipWidthSpin.value()
        txt = str(self.stripIterValue.text())
        if len(txt):
            pars["stripiterations"]= int(txt)
        else:
            pars["stripiterations"] = 0
        pars["stripwidth"]= self.stripWidthSpin.value()
        pars["stripfilterwidth"] = self.stripFilterSpin.value()
        pars["stripanchorsflag"] = int(self.stripAnchorsFlagCheck.isChecked())
        pars["stripanchorslist"] = []
        for spin in self.stripAnchorsList:
            pars["stripanchorslist"].append(spin.value())
        return pars

    def _emitSignal(self, dummy=None):
        ddict= {}
        ddict['event']='ParametersChanged'
        ddict['parameters'] = self.getParameters()
        self.emit(qt.SIGNAL('StripParametersWidgetSignal'), ddict)

class StripBackgroundWidget(qt.QWidget):
    def __init__(self, parent=None):
        qt.QWidget.__init__(self, parent)
        self.setWindowTitle("Strip and SNIP Configuration Window")
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(2)
        self.parametersWidget = StripParametersWidget(self)
        if SCANWINDOW:
            self.graphWidget = ScanWindow.ScanWindow(self)
        elif QWT:
            self.graphWidget = Plot1DQwt.Plot1DQwt(self)
        else:
            self.graphWidget = Plot1DMatplotlib.Plot1DMatplotlib(self)
        try:
            self.graphWidget.fitButton.hide()
            self.graphWidget.scanWindowInfoWidget.hide()
        except:
            pass
        self.mainLayout.addWidget(self.parametersWidget)
        self.mainLayout.addWidget(self.graphWidget)
        self.getParameters = self.parametersWidget.getParameters
        self.setParameters = self.parametersWidget.setParameters
        self._x = None
        self._y = None
        self.connect(self.parametersWidget,
                     qt.SIGNAL('StripParametersWidgetSignal'),
                     self._slot)

    def setData(self, x, y):
        self._x = x
        self._y = y
        self.update()

    def _slot(self, ddict):
        self.update()

    def update(self):
        if self._y is None:
            return

        pars = self.getParameters()

        #smoothed data
        y = numpy.ravel(numpy.array(self._y)).astype(numpy.float)
        ysmooth = SpecfitFuns.SavitskyGolay(y, pars['stripfilterwidth'])        
        f=[0.25,0.5,0.25]
        ysmooth[1:-1] = numpy.convolve(ysmooth,f,mode=0)
        ysmooth[0] = 0.5 *(ysmooth[0] + ysmooth[1])
        ysmooth[-1] = 0.5 * (ysmooth[-1] + ysmooth[-2])

        #loop for anchors
        x = self._x
        niter = pars['stripiterations']
        anchorslist = []
        if pars['stripanchorsflag']:
            if pars['stripanchorslist'] is not None:
                ravelled = x
                for channel in pars['stripanchorslist']:
                    if channel <= ravelled[0]:continue
                    index = numpy.nonzero(ravelled >= channel)[0]
                    if len(index):
                        index = min(index)
                        if index > 0:
                            anchorslist.append(index)
        if niter > 1000:
            stripBackground = SpecfitFuns.subac(ysmooth,
                                            pars['stripconstant'],
                                            niter,
                                            pars['stripwidth'],
                                            anchorslist)
            #final smoothing
            stripBackground = SpecfitFuns.subac(stripBackground,
                                  pars['stripconstant'],
                                  500,1,
                                  anchorslist)
        elif niter > 0:
            stripBackground = SpecfitFuns.subac(ysmooth,
                                            pars['stripconstant'],
                                            niter,
                                            pars['stripwidth'],
                                            anchorslist)
        else:
            stripBackground = 0.0 * ysmooth + ysmooth.min()

        if len(anchorslist) == 0:
            anchorslist = [0, len(ysmooth)-1]
        anchorslist.sort()
        snipBackground = 0.0 * ysmooth
        lastAnchor = 0
        width = pars['snipwidth']
        for anchor in anchorslist:
            if (anchor > lastAnchor) and (anchor < len(ysmooth)):
                snipBackground[lastAnchor:anchor] =\
                            SpecfitFuns.snip1d(ysmooth[lastAnchor:anchor], width, 0)
                lastAnchor = anchor
        if lastAnchor < len(ysmooth):
            snipBackground[lastAnchor:] =\
                            SpecfitFuns.snip1d(ysmooth[lastAnchor:], width, 0)

        self.graphWidget.addCurve(x, y, \
                                  legend='Input Data',\
                                  replace=True,
                                  replot=False)
        self.graphWidget.addCurve(x, stripBackground,\
                                  legend='Strip Background',\
                                  replot=False)
        self.graphWidget.addCurve(x, snipBackground,\
                                  legend='SNIP Background',
                                  replot=True)
        
class StripBackgroundDialog(qt.QDialog):
    def __init__(self, parent=None):
        qt.QDialog.__init__(self, parent)
        self.setWindowTitle("Strip and SNIP Configuration Window")
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(2)
        self.parametersWidget = StripBackgroundWidget(self)
        self.setData = self.parametersWidget.setData
        self.getParameters = self.parametersWidget.getParameters
        self.setParameters = self.parametersWidget.setParameters
        self.mainLayout.addWidget(self.parametersWidget)
        hbox = qt.QWidget(self)
        hboxLayout = qt.QHBoxLayout(hbox)
        hboxLayout.setMargin(0)
        hboxLayout.setSpacing(2)
        self.okButton = qt.QPushButton(hbox)
        self.okButton.setText("OK")
        self.okButton.setAutoDefault(False)   
        self.dismissButton = qt.QPushButton(hbox)
        self.dismissButton.setText("Cancel")
        self.dismissButton.setAutoDefault(False)
        hboxLayout.addWidget(qt.HorizontalSpacer(hbox))
        hboxLayout.addWidget(self.okButton)
        hboxLayout.addWidget(self.dismissButton)
        self.mainLayout.addWidget(hbox)
        self.connect(self.dismissButton, qt.SIGNAL("clicked()"), self.reject)
        self.connect(self.okButton, qt.SIGNAL("clicked()"), self.accept)
        
    def sizeHint(self):
        return qt.QSize(int(1.5*qt.QDialog.sizeHint(self).width()),
                        qt.QDialog.sizeHint(self).height())
        
if __name__ == "__main__":
    a = qt.QApplication(sys.argv)
    qt.QObject.connect(a,qt.SIGNAL("lastWindowClosed()"),a,qt.SLOT("quit()"))
    w = StripBackgroundDialog()
    def mySlot(ddict):
        print(ddict)
    qt.QObject.connect(w.parametersWidget, qt.SIGNAL('StripParametersWidgetSignal'), mySlot)
    x = numpy.arange(1000.).astype(numpy.float32)
    y = 100 + x + 100 * numpy.exp(-0.5*(x-500) * (x-500)/ 30.)
    w.setData(x, y)
    w.exec_()
    #a.exec_()
