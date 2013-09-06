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
__author__ = "V.A. Sole - ESRF BLISS Group"

from PyMca import PyMcaQt as qt
QTVERSION = qt.qVersion()

if QTVERSION < '4.0.0':
    import Qwt5
else:
    from PyQt4 import Qwt5

DEBUG = 0
    
class DoubleSlider(qt.QWidget):
    def __init__(self, parent = None, scale = False):
        qt.QWidget.__init__(self, parent)
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setMargin(6)
        self.mainLayout.setSpacing(1)
        orientation = qt.Qt.Horizontal
        
        self.minSliderContainer = MySlider(self, orientation)
        self.minSlider = self.minSliderContainer.slider
        if scale:
            self.minSlider.setScalePosition(Qwt5.QwtSlider.BottomScale)
        self.minSlider.setRange(0.0, 100.0, 0.01)
        self.minSlider.setValue(0.0)
        self.maxSliderContainer = MySlider(self, orientation)
        self.maxSlider = self.maxSliderContainer.slider
        self.maxSlider.setRange(0.0, 100.0, 0.01)
        self.maxSlider.setValue(100.)
        self.mainLayout.addWidget(self.maxSliderContainer)
        self.mainLayout.addWidget(self.minSliderContainer)
        self.connect(self.minSlider,
                     qt.SIGNAL("valueChanged(double)"),
                     self._sliderChanged)
        self.connect(self.maxSlider,
                     qt.SIGNAL("valueChanged(double)"),
                     self._sliderChanged)

    def __getDict(self):
        ddict = {}
        ddict['event'] = "doubleSliderValueChanged"
        m   = self.minSlider.value()
        M   = self.maxSlider.value()
        if m > M:
            ddict['max'] = m
            ddict['min'] = M
        else:
            ddict['min'] = m
            ddict['max'] = M
        return ddict

    def _sliderChanged(self, value):
        if DEBUG:
            print("DoubleSlider._sliderChanged()")
        ddict = self.__getDict()
        if QTVERSION < '4.0.0':
            self.emit(qt.PYSIGNAL("doubleSliderValueChanged"), (ddict,))
        else:
            self.emit(qt.SIGNAL("doubleSliderValueChanged"), ddict)

    def setMinMax(self, m, M):
        self.minSlider.setValue(m)
        self.maxSlider.setValue(M)

    def getMinMax(self):
        m = self.minSlider.value()
        M = self.maxSlider.value()
        if m > M:
            return M, m
        else:
            return m, M
        

class MySlider(qt.QWidget):
    def __init__(self, parent = None, orientation=qt.Qt.Horizontal):
        qt.QWidget.__init__(self, parent)
        if orientation == qt.Qt.Horizontal:
            alignment = qt.Qt.AlignHCenter | qt.Qt.AlignTop
            layout = qt.QHBoxLayout(self)
        else:
            alignment = qt.Qt.AlignVCenter | qt.Qt.AlignLeft
            layout = qt.QVBoxLayout(self)
        layout.setMargin(0)
        layout.setSpacing(0)
        self.slider = Qwt5.QwtSlider(self, orientation)
        self.label  = qt.QLabel("0", self)
        self.label.setAlignment(alignment)
        self.label.setFixedWidth(self.label.fontMetrics().width('100.99'))

        layout.addWidget(self.slider)
        layout.addWidget(self.label)
        self.connect(self.slider,
                     qt.SIGNAL('valueChanged(double)'),
                     self.setNum)

    def setNum(self, value):
        self.label.setText('%s' % value)

def test():
    app = qt.QApplication([])
    qt.QObject.connect(app,
                       qt.SIGNAL("lastWindowClosed()"),
                       app,
                       qt.SLOT('quit()'))
    
    w = DoubleSlider()
    w.show()
    if QTVERSION < '4.0.0':
        app.exec_loop()
    else:
        app.exec_()

if __name__ == "__main__":
    test()
 
