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
import PyQt4.QtCore as QtCore
import PyQt4.QtGui as QtGui

class HDFInfoCustomEvent(QtCore.QEvent):
    def __init__(self, ddict):
        if ddict is None:
            ddict = {}
        self.dict = ddict
        QtCore.QEvent.__init__(self, QtCore.QEvent.User)
        
class CloseEventNotifyingWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self._notifyCloseEventToWidget = []

    def notifyCloseEventToWidget(self, widget):
        if widget not in self._notifyCloseEventToWidget:
            self._notifyCloseEventToWidget.append(widget)

    def closeEvent(self, event):
        if len(self._notifyCloseEventToWidget):
            for widget in self._notifyCloseEventToWidget:
                ddict={}
                ddict['event'] = 'closeEventSignal'
                ddict['id']    = id(self)
                newEvent = HDFInfoCustomEvent(ddict)
                QtGui.QApplication.postEvent(widget,
                                      newEvent)
            self._notifyCloseEventToWidget = []
        return QtGui.QWidget.closeEvent(self, event)
