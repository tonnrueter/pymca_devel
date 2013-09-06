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
from PyQt4 import QtCore, QtGui
if hasattr(QtCore, 'QStringList'):
    MyQVariant = QtCore.QVariant
else:
    def MyQVariant(x=None):
        return x
from PyMca import NumpyArrayTableModel
import sys

class HorizontalHeader(QtCore.QAbstractItemModel):
    def __init__(self, parent=None):
        QtGui.QHeaderView.__init__(self, parent)

    def columnCount(self, modelIndex):
        return self.parent().columnCount()

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole:
            return MyQVariant("%d" % section)
        return MyQVariant()

class VerticalHeader(QtCore.QAbstractItemModel):
    def __init__(self, parent=None):
        QtGui.QHeaderView.__init__(self, parent)

    def rowCount(self, modelIndex):
        return self.parent().rowCount()

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole:
            return MyQVariant("%d" % section)
        return MyQVariant()

class NumpyArrayTableView(QtGui.QTableView):
    def __init__(self, parent=None):
        QtGui.QTableView.__init__(self, parent)
        self._model = NumpyArrayTableModel.NumpyArrayTableModel(self)
        self.setModel(self._model)
        self._horizontalHeaderModel = HorizontalHeader(self._model)
        self._verticalHeaderModel = VerticalHeader(self._model)
        self.horizontalHeader().setModel(self._horizontalHeaderModel)
        self.verticalHeader().setModel(self._verticalHeaderModel)

    def setArrayData(self, data):
        t = "%s" % data.dtype
        if '|' in t:
            fmt = "%s"
        else:
            fmt = "%g"
        self._model.setFormat(fmt)
        self._model.setArrayData(data)
        #some linux distributions need this call
        self.setModel(self._model)
        if sys.platform not in ['win32']:
            self._horizontalHeaderModel = HorizontalHeader(self._model)
            self._verticalHeaderModel = VerticalHeader(self._model)
        self.horizontalHeader().setModel(self._horizontalHeaderModel)
        self.verticalHeader().setModel(self._verticalHeaderModel)
        
    def setCurrentArrayIndex(self, index):
        return self._model.setCurrentArrayIndex(index)

if __name__ == "__main__":
    import numpy
    a = QtGui.QApplication([])
    d = numpy.random.normal(0,1, (5, 1000,1000))
    for i in range(5):
        d[i, :, :] += i
    w = NumpyArrayTableView()
    w.setArrayData(d)
    w.show()
    a.exec_()
