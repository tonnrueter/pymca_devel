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
import os
from PyMca import PyMcaQt as qt
QTVERSION = qt.qVersion()
import PyMca.PyMca_Icons as icons
from PyMca import spswrap as sps
from PyMca import PyMcaDirs

DEBUG = 0

class QSourceSelector(qt.QWidget):
    def __init__(self, parent=None, filetypelist=None):
        qt.QWidget.__init__(self, parent)
        self.mainLayout= qt.QVBoxLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(0)
        if filetypelist is None:
            self.fileTypeList = ["Spec Files (*mca)",
                                "Spec Files (*dat)",
                                "Spec Files (*spec)",
                                "SPE Files (*SPE)",
                                "EDF Files (*edf)",
                                "EDF Files (*ccd)",
                                "CSV Files (*csv)",
                                "All Files (*)"]
        else:
            self.fileTypeList = filetypelist
        self.lastFileFilter = self.fileTypeList[0]

        # --- file combo/open/close
        self.lastInputDir = PyMcaDirs.inputDir
        self.fileWidget= qt.QWidget(self)
        fileWidgetLayout= qt.QHBoxLayout(self.fileWidget)
        fileWidgetLayout.setMargin(0)
        fileWidgetLayout.setSpacing(0)
        self.fileCombo  = qt.QComboBox(self.fileWidget)
        self.fileCombo.setEditable(0)
        self.mapCombo= {}
        openButton= qt.QToolButton(self.fileWidget)
        if QTVERSION < '4.0.0':
            self.openIcon= qt.QIconSet(qt.QPixmap(icons.fileopen))
            self.closeIcon= qt.QIconSet(qt.QPixmap(icons.fileclose))
            self.specIcon= qt.QIconSet(qt.QPixmap(icons.spec))
            self.reloadIcon = qt.QIconSet(qt.QPixmap(icons.reload))
            openButton.setIconSet(self.openIcon)
            openButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
            closeButton= qt.QToolButton(self.fileWidget)
            closeButton.setIconSet(self.closeIcon)
            specButton= qt.QToolButton(self.fileWidget)
            specButton.setIconSet(self.specIcon)
            refreshButton= qt.QToolButton(self.fileWidget)
            refreshButton.setIconSet(self.reloadIcon)
        else:
            self.openIcon   = qt.QIcon(qt.QPixmap(icons.fileopen))
            self.closeIcon  = qt.QIcon(qt.QPixmap(icons.fileclose))
            self.reloadIcon = qt.QIcon(qt.QPixmap(icons.reload))
            self.specIcon   = qt.QIcon(qt.QPixmap(icons.spec))

            openButton.setIcon(self.openIcon)
            openButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
            openButton.setToolTip("Open new file data source")

            closeButton= qt.QToolButton(self.fileWidget)
            closeButton.setIcon(self.closeIcon)
            closeButton.setToolTip("Close current data source")

            refreshButton= qt.QToolButton(self.fileWidget)
            refreshButton.setIcon(self.reloadIcon)
            refreshButton.setToolTip("Refresh data source")

            specButton= qt.QToolButton(self.fileWidget)
            specButton.setIcon(self.specIcon)
            specButton.setToolTip("Open new shared memory source")

        closeButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        specButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        refreshButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))

        self.connect(openButton, qt.SIGNAL("clicked()"), self._openFileSlot)
        self.connect(closeButton, qt.SIGNAL("clicked()"), self.closeFile)
        self.connect(refreshButton, qt.SIGNAL("clicked()"),
                         self._reload)
            
        self.connect(specButton, qt.SIGNAL("clicked()"), self.openSpec)
        self.connect(self.fileCombo, qt.SIGNAL("activated(const QString &)"),
                                                     self._fileSelection)

        fileWidgetLayout.addWidget(self.fileCombo)
        fileWidgetLayout.addWidget(openButton)            
        fileWidgetLayout.addWidget(closeButton)
        fileWidgetLayout.addWidget(specButton)
        if sys.platform == "win32":specButton.hide()
        fileWidgetLayout.addWidget(refreshButton)
        self.specButton = specButton
        self.mainLayout.addWidget(self.fileWidget)

    def _reload(self):
        if DEBUG:
            print("_reload called")
        qstring = self.fileCombo.currentText()
        if not len(qstring):
            return

        key = qt.safe_str(qstring)
        ddict = {}
        ddict["event"] = "SourceReloaded"
        ddict["combokey"] = key
        ddict["sourcelist"] = self.mapCombo[key] * 1
        if QTVERSION < '4.0.0':
            self.emit(qt.PYSIGNAL("SourceSelectorSignal"), (ddict,))    
        else:
            self.emit(qt.SIGNAL("SourceSelectorSignal"), ddict)

    def _openFileSlot(self):
        self.openFile(None, None)
    
    def openSource(self, sourcename, specsession=None):
        if specsession is None:
            if sourcename in sps.getspeclist():
                specsession=True 
            else:
                specsession=False
        self.openFile(sourcename, specsession=specsession)

    def openFile(self, filename=None,justloaded=None, specsession = False):
        if DEBUG:
            print("openfile = ",filename)
        staticDialog = False
        if not specsession:
            if justloaded is None: justloaded = True
            if filename is None:
                if self.lastInputDir is not None:
                    if not os.path.exists(self.lastInputDir):
                        self.lastInputDir = None
                wdir = self.lastInputDir
                if wdir is None:wdir = os.getcwd()
                if QTVERSION < '4.0.0':
                    filetypes = ""
                    for filetype in self.fileTypeList:
                        filetypes += filetype+"\n"
                    if sys.platform == 'win32':
                        filelist = qt.QFileDialog.getOpenFileNames(
                                    filetypes,
                                    wdir,
                                    self,
                                    "openFile",
                                    "Open a new EdfFile")
                    else:
                        filedialog = qt.QFileDialog(self,
                                                    "Open new EdfFile(s)",
                                                    1)
                        if self.lastInputDir is not None:
                            filedialog.setDir(self.lastInputDir)
                        filedialog.setMode(filedialog.ExistingFiles)
                        filedialog.setFilters(filetypes)           
                        if filedialog.exec_loop() == qt.QDialog.Accepted:
                            filelist= filedialog.selectedFiles()
                        else:
                            return              
                else:
                    #if sys.platform == 'win32':
                    #if (QTVERSION < '4.3.0') and (sys.platform != 'darwin'):
                    if (sys.version < '3.0') and PyMcaDirs.nativeFileDialogs:
                        filetypes = ""
                        for filetype in self.fileTypeList:
                            filetypes += filetype+"\n"
                        if sys.version < '3.0':
                            filelist = qt.QFileDialog.getOpenFileNames(self,
                                    "Open a new source file",
                                    wdir,
                                    filetypes,
                                    self.lastFileFilter)
                        else:
                            filelist, self.lastFileFilter =\
                                    qt.QFileDialog.getOpenFileNamesAndFilter(\
                                    self,
                                    "Open a new source file",
                                    wdir,
                                    filetypes,
                                    self.lastFileFilter)
                        staticDialog = True
                    else:
                        fdialog = qt.QFileDialog(self)
                        fdialog.setModal(True)
                        fdialog.setWindowTitle("Open a new source file")
                        if hasattr(qt, "QStringList"):
                            strlist = qt.QStringList()
                        else:
                            strlist = []
                        for filetype in self.fileTypeList:
                            strlist.append(filetype)
                        fdialog.setFilters(strlist)
                        fdialog.selectFilter(self.lastFileFilter)
                        fdialog.setFileMode(fdialog.ExistingFiles)
                        fdialog.setDirectory(wdir)
                        ret = fdialog.exec_()
                        if ret == qt.QDialog.Accepted:
                            filelist = fdialog.selectedFiles()
                            self.lastFileFilter = qt.safe_str(\
                                                    fdialog.selectedFilter())
                            fdialog.close()
                            del fdialog                        
                        else:
                            fdialog.close()
                            del fdialog
                            return            
                #filelist.sort()
                filename=[]
                for f in filelist:
                    filename.append(qt.safe_str(f))
                if not len(filename):
                    return
                if len(filename):
                    self.lastInputDir  = os.path.dirname(filename[0])
                    PyMcaDirs.inputDir = os.path.dirname(filename[0])
                    if staticDialog:
                        if len(filename[0]) > 3:
                            #figure out the selected filter
                            extension = filename[0][-3:]
                            self.lastFileFilter = self.fileTypeList[-1]
                            for fileFilter in self.fileTypeList:
                                if extension == fileFilter[-4:-1]:
                                    self.lastFileFilter = fileFilter
                                    break
                justloaded = True
            if justloaded:
                if type(filename) != type([]):
                    filename = [filename]
            if not os.path.exists(filename[0]):
                if '%' not in filename[0]:
                    raise IOError("File %s does not exist" % filename[0])

            #check if it is a stack
            if len(filename) > 1:
                key = "STACK from %s to %s" % (filename[0], filename[-1])
            else:
                key = os.path.basename(filename[0])
        else:
            key = filename
            if key not in sps.getspeclist():
                qt.QMessageBox.critical(self,
                                    "SPS Error", 
                                    "No shared memory source named %s" % key)
                return
        ddict = {}
        ddict["event"] = "NewSourceSelected"
        if key in self.mapCombo.keys():
            if self.mapCombo[key] == filename:
                #Reloaded event
                ddict["event"] = "SourceReloaded"
            else:
                i = 0
                while key in self.mapCombo.keys():
                    key += "_%d" % i
        ddict["combokey"]   = key
        ddict["sourcelist"] = filename
        self.mapCombo[key] = filename
        if ddict["event"] =="NewSourceSelected":
            nitems = self.fileCombo.count()
            if QTVERSION < '4.0.0':
                self.fileCombo.insertItem(key)
                self.fileCombo.setCurrentItem(nitems)
            else:
                self.fileCombo.insertItem(nitems, key)
                self.fileCombo.setCurrentIndex(nitems)
        else:
            if hasattr(qt, "QString"):
                nitem = self.fileCombo.findText(qt.QString(key))
            else:
                nitem = self.fileCombo.findText(key)
            self.fileCombo.setCurrentIndex(nitem)
        if QTVERSION < '4.0.0':
            self.emit(qt.PYSIGNAL("SourceSelectorSignal"), (ddict,))        
        else:
            self.emit(qt.SIGNAL("SourceSelectorSignal"), ddict)

    def closeFile(self):
        if DEBUG:
            print("closeFile called")
        #get current combobox key
        qstring = self.fileCombo.currentText()
        if not len(qstring):
            return
        key = qt.safe_str(qstring)
        ddict = {}
        ddict["event"] = "SourceClosed"
        ddict["combokey"] = key
        ddict["sourcelist"] = self.mapCombo[key] * 1
        if hasattr(qt, "QString"):
            nitem = self.fileCombo.findText(qt.QString(key))
        else:
            nitem = self.fileCombo.findText(key)
        self.fileCombo.removeItem(nitem)
        del self.mapCombo[key]
        if QTVERSION < '4.0.0':
            self.emit(qt.PYSIGNAL("SourceSelectorSignal"), (ddict,))    
        else:
            self.emit(qt.SIGNAL("SourceSelectorSignal"), ddict)

    def openSpec(self):
        speclist = sps.getspeclist()
        if not len(speclist):
            qt.QMessageBox.information(self,
                    "No SPEC Shared Memory Found", 
                    "No shared memory source available")
            return
        if QTVERSION < '4.0.0':
            print("should I keep Qt3 version?")
            return
        menu = qt.QMenu()
        for spec in speclist:
            if hasattr(qt, "QString"):
                menu.addAction(qt.QString(spec), 
                        lambda i=spec:self.openFile(i, specsession=True))
            else:
                menu.addAction(spec, 
                        lambda i=spec:self.openFile(i, specsession=True))
        menu.exec_(self.cursor().pos())

    def _fileSelection(self, qstring):
        if DEBUG:
            print("file selected ", qstring)
        key = str(qstring)
        ddict = {}
        ddict["event"] = "SourceSelected"
        ddict["combokey"] = key
        ddict["sourcelist"] = self.mapCombo[key]
        if QTVERSION < '4.0.0':
            self.emit(qt.PYSIGNAL("SourceSelectorSignal"), (ddict,))    
        else:
            self.emit(qt.SIGNAL("SourceSelectorSignal"), ddict)

def test():
    a = qt.QApplication(sys.argv)
    #new access
    from PyMca import QDataSource
    w= QSourceSelector()
    def mySlot(ddict):
        print(ddict)
        if ddict["event"] == "NewSourceSelected":
            d = QDataSource.QDataSource(ddict["sourcelist"][0])
            w.specfileWidget.setDataSource(d)
            if QTVERSION < '4.0.0':
                a.connect(w, qt.PYSIGNAL("SourceSelectorSignal"),
                      mySlot)
            else:
                a.connect(w, qt.SIGNAL("SourceSelectorSignal"),
                       mySlot)

        
    qt.QObject.connect(a, qt.SIGNAL("lastWindowClosed()"),
              a, qt.SLOT("quit()"))

    if QTVERSION < '4.0.0':
        w.show()
        a.exec_loop()
    else:
        w.show()
        a.exec_()


if __name__=="__main__":
    test()
