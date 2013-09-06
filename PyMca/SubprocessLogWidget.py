#/*##########################################################################
# Copyright (C) 2004-2013 European Synchrotron Radiation Facility
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
import os
import subprocess
import time
from PyMca import PyMcaQt as qt

class SubprocessLogWidget(qt.QWidget):
    def __init__(self, parent=None, args=None):
        qt.QWidget.__init__(self, parent)
        self.setWindowTitle("Subprocess Log Widget")
        self.mainLayout = qt.QVBoxLayout(self)
        self._p = None
        self.__timer = qt.QTimer()
        self._args = args
        self.connect(self.__timer, qt.SIGNAL("timeout()"), self._timerSlot)
        self.logWidget = qt.QTextEdit(self)
        self.logWidget.setReadOnly(True)
        self.mainLayout.addWidget(self.logWidget)

    def setSubprocessArgs(self, args):
        self._args = args

    def start(self, args=None, timing=0.1):
        if args is None:
            if self._args is None:
                raise ValueError("Subprocess command not defined")
            else:
                self._args = args
        else:
            self._args = args
        self._startTimer(timing=timing)

    def stop(self):
        if self.isSubprocessRunning():
            #print("MESSAGE TO KILL PROCESS")
            #print("HOW TO KILL IT IN A GOOD WAY?")
            self._p.kill()

    def isSubprocessRunning(self):
        running = False
        if self._p is not None:
            if self._p.poll() is None:
                running = True
        return running

    def _startTimer(self, timing=0.1):
        if self._args is None:
            raise ValueError("Subprocess command not defined")
        self._p = subprocess.Popen(self._args,
                                   bufsize=0,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True)
        ddict = {}
        ddict["subprocess"] = self._p
        ddict["event"] = "ProcessStarted"
        self.emit(qt.SIGNAL("SubprocessLogWidgetSignal"), ddict)
        self.__timer.start(timing)

    def _timerSlot(self):
        ddict = {}
        ddict["subprocess"] = self._p
        if self._p.poll() is None:
            # process did not finish yet
            line = self._p.stdout.readline()
            if len(line) > 1:
                self.logWidget.append(line[:-1])
                qt.qApp.processEvents()
            ddict["event"] = "ProcessRunning"
        else:
            self.__timer.stop()
            returnCode = self._p.returncode
            ddict['event'] = "ProcessFinished"
            ddict['code'] = returnCode
            ddict["message"] = []
            if returnCode == 0:
                line = self._p.stdout.readline()
                while len(line) > 1:
                    self.logWidget.append(line[:-1])
                    line = self._p.stdout.readline()
            else:
                line = self._p.stderr.readline()
                while len(line) > 1:
                    ddict["message"].append(line)
                    self.logWidget.append(line[:-1])
                    line = self._p.stderr.readline()
            self._p = None
        self.emit(qt.SIGNAL("SubprocessLogWidgetSignal"), ddict)

    def clear(self):
        self.logWidget.clear()

    def append(self, text):
        self.logWidget.append(text)

    def closeEvent(self, event):
        if self._p is not None:
            try:
                self.stop()
            except:
                # this may happen if the process finished in the mean time
                pass
        qt.QWidget.closeEvent(self, event)
