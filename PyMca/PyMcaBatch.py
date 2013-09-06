#!/usr/bin/env python
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
__revision__ = "$Revision: 4.12$"
import sys
import os
import time
import subprocess

from PyMca import PyMcaQt as qt

QTVERSION = qt.qVersion()
if QTVERSION > '4.0.0':
    try:
        import h5py
        from PyMca import NexusDataSource
        from PyMca import QNexusWidget
        from PyMca import HDF5Selection
        HDF5SUPPORT = True
    except ImportError:
        HDF5SUPPORT = False
else:
    HDF5SUPPORT = False
from PyMca import ConfigDict
from PyMca import McaAdvancedFitBatch
from PyMca import EdfFileLayer
from PyMca import SpecFileLayer
from PyMca.PyMca_Icons import IconDict
from PyMca import McaCustomEvent
from PyMca import EdfFileSimpleViewer
from PyMca import QtMcaAdvancedFitReport
from PyMca import HtmlIndex
from PyMca import PyMcaDirs
from PyMca import PyMcaBatchBuildOutput

ROIWIDTH = 100.
DEBUG = 0

class McaBatchGUI(qt.QWidget):
    def __init__(self,parent=None,name="PyMca batch fitting",fl=None,
                filelist=None,config=None,outputdir=None, actions=0):
        if QTVERSION < '4.0.0':
            if fl is None:
                fl = qt.Qt.WDestructiveClose
            qt.QWidget.__init__(self,parent,name,fl)
            self.setIcon(qt.QPixmap(IconDict['gioconda16']))
            self.setCaption(name)
        else:
            qt.QWidget.__init__(self, parent)
            self.setWindowTitle(name)
            self.setWindowIcon(qt.QIcon(qt.QPixmap(IconDict['gioconda16'])))
        self._layout = qt.QVBoxLayout(self)
        self._layout.setMargin(0)
        self._layout.setSpacing(0)
        self._edfSimpleViewer = None
        self._timer = None
        self._processList = []
        self._selection = None
        self.__build(actions)               
        if filelist is None: filelist = []
        self.inputDir   = None
        self.outputDir  = None
        if outputdir is not None:
            if os.path.exists(outputdir):
                self.outputDir = outputdir
            else:
                qt.QMessageBox.information(self, "INFO",
                    "Directory %s does not exists\nUsing %s"% (outputdir, self.outputDir))
                
        self.configFile = None
        self.setFileList(filelist)
        self.setConfigFile(config)
        self.setOutputDir(self.outputDir)
    
    def __build(self,actions):
        self.__grid= qt.QWidget(self)
        self._layout.addWidget(self.__grid)
        #self.__grid.setGeometry(qt.QRect(30,30,288,156))
        if QTVERSION < '4.0.0':
            grid       = qt.QGridLayout(self.__grid,3,3,11,6)
            grid.setColStretch(0,0)
            grid.setColStretch(1,1)
            grid.setColStretch(2,0)
        else:
            grid       = qt.QGridLayout(self.__grid)
            grid.setMargin(11)
            grid.setSpacing(6)
        #input list
        listrow  = 0
        listlabel   = qt.QLabel(self.__grid)
        listlabel.setText("Input File list:")
        if QTVERSION < '4.0.0':
            listlabel.setAlignment(qt.QLabel.WordBreak | qt.QLabel.AlignVCenter)
            self.__listView   = qt.QTextView(self.__grid)
            self.__listView.setMaximumHeight(30*listlabel.sizeHint().height())
        else:
            self.__listView   = qt.QTextEdit(self.__grid)
            self.__listView.setMaximumHeight(30*listlabel.sizeHint().height())
        self.__listButton = qt.QPushButton(self.__grid)
        self.__listButton.setText('Browse')
        self.connect(self.__listButton,qt.SIGNAL('clicked()'),self.browseList) 
        grid.addWidget(listlabel,        listrow, 0, qt.Qt.AlignTop|qt.Qt.AlignLeft)
        grid.addWidget(self.__listView,  listrow, 1)
        grid.addWidget(self.__listButton,listrow, 2, qt.Qt.AlignTop|qt.Qt.AlignRight)

        if HDF5SUPPORT:
            self._hdf5Widget = HDF5Selection.HDF5Selection(self)
            grid.addWidget(self._hdf5Widget, listrow+1, 0, 1, 3)
            row_offset = 1
            self._hdf5Widget.hide()
        else:
            row_offset = 0        
        #config file
        configrow = 1 + row_offset
        configlabel = qt.QLabel(self.__grid)
        configlabel.setText("Fit Configuration File:")
        if QTVERSION < '4.0.0':
            configlabel.setAlignment(qt.QLabel.WordBreak | qt.QLabel.AlignVCenter)
        self.__configLine = qt.QLineEdit(self.__grid)
        self.__configLine.setReadOnly(True)
        self.__configButton = qt.QPushButton(self.__grid)
        self.__configButton.setText('Browse')
        self.connect(self.__configButton,qt.SIGNAL('clicked()'),self.browseConfig) 
        grid.addWidget(configlabel,         configrow, 0, qt.Qt.AlignLeft)
        grid.addWidget(self.__configLine,   configrow, 1)
        grid.addWidget(self.__configButton, configrow, 2, qt.Qt.AlignLeft)


        #output dir
        outrow    = 2 + row_offset
        outlabel   = qt.QLabel(self.__grid)
        outlabel.setText("Output dir:")
        if QTVERSION < '4.0.0':
            outlabel.setAlignment(qt.QLabel.WordBreak | qt.QLabel.AlignVCenter)
        self.__outLine = qt.QLineEdit(self.__grid)
        self.__outLine.setReadOnly(True)
        #self.__outLine.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Maximum, qt.QSizePolicy.Fixed))
        self.__outButton = qt.QPushButton(self.__grid)
        self.__outButton.setText('Browse')
        self.connect(self.__outButton,qt.SIGNAL('clicked()'),self.browseOutputDir) 
        grid.addWidget(outlabel,         outrow, 0, qt.Qt.AlignLeft)
        grid.addWidget(self.__outLine,   outrow, 1)
        grid.addWidget(self.__outButton, outrow, 2, qt.Qt.AlignLeft)


        box  = qt.QWidget(self)
        box.l = qt.QHBoxLayout(box)
        box.l.setMargin(11)
        box.l.setSpacing(0)
        vbox1 = qt.QWidget(box)
        vbox1.l = qt.QVBoxLayout(vbox1)
        vbox1.l.setMargin(0)
        vbox1.l.setSpacing(0)
        box.l.addWidget(vbox1)
        vbox2 = qt.QWidget(box)
        vbox2.l = qt.QVBoxLayout(vbox2)
        vbox2.l.setMargin(0)
        vbox2.l.setSpacing(0)

        box.l.addWidget(vbox2)
        vbox3 = qt.QWidget(box)
        vbox3.l = qt.QVBoxLayout(vbox3)
        vbox3.l.setMargin(0)
        vbox3.l.setSpacing(0)

        box.l.addWidget(vbox3)
        self.__fitBox = qt.QCheckBox(vbox1)
        self.__fitBox.setText('Generate .fit Files')
        palette = self.__fitBox.palette()
        #if QTVERSION < '4.0.0':
        #    palette.setDisabled(palette.active())
        #else:
        #    print("palette set disabled")
        self.__fitBox.setChecked(False)
        self.__fitBox.setEnabled(True)
        vbox1.l.addWidget(self.__fitBox)

        self.__imgBox = qt.QCheckBox(vbox2)
        self.__imgBox.setText('Generate Peak Images')
        palette = self.__imgBox.palette()
        if QTVERSION < '4.0.0':
            palette.setDisabled(palette.active())
        else:
            if DEBUG:
                print("palette set disabled")
        self.__imgBox.setChecked(True)
        self.__imgBox.setEnabled(False)
        vbox2.l.addWidget(self.__imgBox)
        """
        self.__specBox = qt.QCheckBox(box)
        self.__specBox.setText('Generate Peak Specfile')
        palette = self.__specBox.palette()
        palette.setDisabled(palette.active())
        self.__specBox.setChecked(False)
        self.__specBox.setEnabled(False)
        """
        self.__htmlBox = qt.QCheckBox(vbox3)
        self.__htmlBox.setText('Generate Report (SLOW!)')
        #palette = self.__htmlBox.palette()
        #palette.setDisabled(palette.active())
        self.__htmlBox.setChecked(False)
        self.__htmlBox.setEnabled(True)
        vbox3.l.addWidget(self.__htmlBox)
        
        #report options
        #reportBox = qt.QHBox(self)
        self.__tableBox = qt.QCheckBox(vbox1)
        self.__tableBox.setText('Table in Report')
        palette = self.__tableBox.palette()
        #if QTVERSION < '4.0.0':
        #    palette.setDisabled(palette.active())
        #else:
        #    print("palette set disabled")
        self.__tableBox.setChecked(True)
        self.__tableBox.setEnabled(False)
        vbox1.l.addWidget(self.__tableBox)

        self.__extendedTable = qt.QCheckBox(vbox2)
        self.__extendedTable.setText('Extended Table')
        self.__extendedTable.setChecked(True)
        self.__extendedTable.setEnabled(False)
        vbox2.l.addWidget(self.__extendedTable)
        
        self.__concentrationsBox = qt.QCheckBox(vbox3)
        self.__concentrationsBox.setText('Concentrations')
        self.__concentrationsBox.setChecked(False)
        self.__concentrationsBox.setEnabled(True)
        #self.__concentrationsBox.setEnabled(False)
        vbox3.l.addWidget(self.__concentrationsBox)
        self._layout.addWidget(box)
        
        # other stuff
        bigbox   = qt.QWidget(self)
        bigbox.l = qt.QHBoxLayout(bigbox)
        bigbox.l.setMargin(0)
        bigbox.l.setSpacing(0)
        
        vBox = qt.QWidget(bigbox)
        vBox.l = qt.QVBoxLayout(vBox)
        vBox.l.setMargin(0)
        vBox.l.setSpacing(2)
        bigbox.l.addWidget(vBox)

        if 0:
            #These options are obsolete now
            self.__overwrite = qt.QCheckBox(vBox)
            self.__overwrite.setText('Overwrite Fit Files')
            self.__overwrite.setChecked(True)
            vBox.l.addWidget(self.__overwrite)
            
            self.__useExisting = qt.QCheckBox(vBox)
            self.__useExisting.setText('Use Existing Fit Files')
            self.__useExisting.setChecked(False)
            vBox.l.addWidget(self.__useExisting)
        
            self.connect(self.__overwrite,   qt.SIGNAL("clicked()"),
                                             self.__clickSignal0)
            self.connect(self.__useExisting, qt.SIGNAL("clicked()"),
                                             self.__clickSignal1)
        self.connect(self.__concentrationsBox, qt.SIGNAL("clicked()"),
                                               self.__clickSignal2)
        self.connect(self.__htmlBox, qt.SIGNAL("clicked()"),
                                               self.__clickSignal3)
        
        boxStep0   = qt.QWidget(bigbox)
        boxStep0.l = qt.QVBoxLayout(boxStep0)
        boxStep  = qt.QWidget(boxStep0)
        boxStep.l= qt.QHBoxLayout(boxStep)
        boxStep.l.setMargin(0)
        boxStep.l.setSpacing(0)
        boxStep0.l.addWidget(boxStep)
        bigbox.l.addWidget(boxStep0)

        if 0:
            self.__boxFStep   = qt.QWidget(boxStep)
            boxFStep = self.__boxFStep
            boxFStep.l = qt.QHBoxLayout(boxFStep)
            boxFStep.l.setMargin(0)
            boxFStep.l.setSpacing(0)
            boxStep.l.addWidget(boxFStep)
            label= qt.QLabel(boxFStep)
            label.setText("File Step:")
            self.__fileSpin = qt.QSpinBox(boxFStep)
            if QTVERSION < '4.0.0':
                self.__fileSpin.setMinValue(1)
                self.__fileSpin.setMaxValue(10)
            else:
                self.__fileSpin.setMinimum(1)
                self.__fileSpin.setMaximum(10)
            self.__fileSpin.setValue(1)
            boxFStep.l.addWidget(label)
            boxFStep.l.addWidget(self.__fileSpin)


            self.__boxMStep   = qt.QWidget(boxStep0)
            boxMStep = self.__boxMStep
            boxMStep.l = qt.QHBoxLayout(boxMStep)
            boxMStep.l.setMargin(0)
            boxMStep.l.setSpacing(0)
            boxStep0.l.addWidget(boxMStep)
            
            label= qt.QLabel(boxMStep)
            label.setText("MCA Step:")
            self.__mcaSpin = qt.QSpinBox(boxMStep)
            if QTVERSION < '4.0.0':
                self.__mcaSpin.setMinValue(1)
                self.__mcaSpin.setMaxValue(10)
            else:
                self.__mcaSpin.setMinimum(1)
                self.__mcaSpin.setMaximum(10)
            self.__mcaSpin.setValue(1)

            boxMStep.l.addWidget(label)
            boxMStep.l.addWidget(self.__mcaSpin)
        

        
        #box2 = qt.QHBox(self)
        self.__roiBox = qt.QCheckBox(vBox)
        self.__roiBox.setText('ROI Fitting Mode')
        vBox.l.addWidget(self.__roiBox)
        #box3 = qt.QHBox(box2)
        self.__box3 = qt.QWidget(boxStep0)
        box3 = self.__box3
        box3.l = qt.QHBoxLayout(box3)
        box3.l.setMargin(0)
        box3.l.setSpacing(0)
        boxStep0.l.addWidget(box3)

        
        label= qt.QLabel(box3)
        label.setText("ROI Width (eV):")
        self.__roiSpin = qt.QSpinBox(box3)
        if QTVERSION < '4.0.0':
            self.__roiSpin.setMinValue(10)
            self.__roiSpin.setMaxValue(1000)
        else:
            self.__roiSpin.setMinimum(10)
            self.__roiSpin.setMaximum(1000)
        self.__roiSpin.setValue(ROIWIDTH)
        box3.l.addWidget(label)
        box3.l.addWidget(self.__roiSpin)


        #BATCH SPLITTING
        self.__splitBox = qt.QCheckBox(vBox)
        self.__splitBox.setText('EXPERIMENTAL: Use several processes')
        vBox.l.addWidget(self.__splitBox)
        #box3 = qt.QHBox(box2)
        self.__box4 = qt.QWidget(boxStep0)
        box4 = self.__box4
        box4.l = qt.QHBoxLayout(box4)
        box4.l.setMargin(0)
        box4.l.setSpacing(0)
        boxStep0.l.addWidget(box4)
        
        label= qt.QLabel(box4)
        label.setText("Number of processes:")
        self.__splitSpin = qt.QSpinBox(box4)
        if QTVERSION < '4.0.0':
            self.__splitSpin.setMinValue(1)
            self.__splitSpin.setMaxValue(1000)
        else:
            self.__splitSpin.setMinimum(1)
            self.__splitSpin.setMaximum(1000)
        self.__splitSpin.setValue(2)
        box4.l.addWidget(label)
        box4.l.addWidget(self.__splitSpin)
        
        self._layout.addWidget(bigbox)

        if actions: self.__buildActions()


    def __clickSignal0(self):
        if self.__overwrite.isChecked():
            self.__useExisting.setChecked(0)
        else:
            self.__useExisting.setChecked(1)

    def __clickSignal1(self):
        if self.__useExisting.isChecked():
            self.__overwrite.setChecked(0)
        else:
            self.__overwrite.setChecked(1)

    def __clickSignal2(self):
        #self.__tableBox.setEnabled(True)
        pass

    def __clickSignal3(self):
        if self.__htmlBox.isChecked():
            self.__tableBox.setEnabled(True)
            #self.__concentrationsBox.setEnabled(True)
            self.__fitBox.setChecked(True)
            self.__fitBox.setEnabled(False)
        else:
            self.__tableBox.setEnabled(False)
            #self.__concentrationsBox.setEnabled(False)
            self.__fitBox.setChecked(False)
            self.__fitBox.setEnabled(True)


    def __buildActions(self):
        box = qt.QWidget(self)
        box.l = qt.QHBoxLayout(box)
        box.l.addWidget(qt.HorizontalSpacer(box))
        self.__dismissButton = qt.QPushButton(box)
        box.l.addWidget(self.__dismissButton)
        box.l.addWidget(qt.HorizontalSpacer(box))
        self.__dismissButton.setText("Close")
        self.__startButton   = qt.QPushButton(box)
        box.l.addWidget(self.__startButton)
        box.l.addWidget(qt.HorizontalSpacer(box))
        self.__startButton.setText("Start")
        self.connect(self.__dismissButton,qt.SIGNAL("clicked()"),self.close)
        self.connect(self.__startButton,qt.SIGNAL("clicked()"),self.start)
        self._layout.addWidget(box)

    def close(self):
        if self._edfSimpleViewer is not None:
            self._edfSimpleViewer.close()
            self._edfSimpleViewer = None
        qt.QWidget.close(self)

    def setFileList(self,filelist=None):
        if filelist is None:filelist = []
        if True or self.__goodFileList(filelist):
            text = ""
            oldtype = None
            #do not sort the file list
            #respect user choice
            #filelist.sort()
            for file in filelist:
                filetype = self.__getFileType(file)
                if filetype is None:return
                if oldtype  is None: oldtype = filetype
                if oldtype != filetype:
                    qt.QMessageBox.critical(self, "ERROR",
                        "Type %s does not match type %s on\n%s"% (filetype,oldtype,file))
                    return           
                text += "%s\n" % file

            if len(filelist):
                if HDF5SUPPORT:
                    if not h5py.is_hdf5(filelist[0]):
                        self._hdf5Widget.hide()
                        self._selection=None
                    else:
                        dialog = qt.QDialog(self)
                        dialog.setWindowTitle('Select your data set')
                        dialog.mainLayout = qt.QVBoxLayout(dialog)
                        dialog.mainLayout.setMargin(0)
                        dialog.mainLayout.setSpacing(0)
                        datasource = NexusDataSource.NexusDataSource(filelist[0])                        
                        nexusWidget = QNexusWidget.QNexusWidget(dialog)
                        nexusWidget.buttons.hide()
                        nexusWidget.setDataSource(datasource)
                        button = qt.QPushButton(dialog)
                        button.setText("Done")
                        button.setAutoDefault(True)
                        qt.QObject.connect(button,
                                           qt.SIGNAL('clicked()'),
                                           dialog.accept)
                        dialog.mainLayout.addWidget(nexusWidget)
                        dialog.mainLayout.addWidget(button)
                        ret = dialog.exec_()
                        cntSelection = nexusWidget.cntTable.getCounterSelection()
                        cntlist = cntSelection['cntlist']
                        if not len(cntlist):
                            text = "No dataset selection"
                            self.showMessage(text)
                            self.__listView.clear()
                            return
                        if not len(cntSelection['y']):
                            text = "No dataset selected as y"
                            self.showMessage(text)
                            self.__listView.clear()
                            return
                        datasource = None
                        selection = {}
                        selection['x'] = []
                        selection['y'] = []
                        selection['m'] = []
                        for key in ['x', 'y', 'm']:
                            if len(cntSelection[key]):
                                for idx in cntSelection[key]:
                                    selection[key].append(cntlist[idx])
                        self._selection = selection
                        self._hdf5Widget.setSelection(selection)
                        #they are not used yet
                        self._hdf5Widget.selectionWidgetsDict['x'].hide()
                        self._hdf5Widget.selectionWidgetsDict['m'].hide()
                        if self._hdf5Widget.isHidden():
                            self._hdf5Widget.show()
                elif filelist[0][-3:].lower() in ['.h5', 'nxs', 'hdf']:
                            text = "Warning, this looks as an HDF5 file "
                            text += "but you do not have HDF5 support."
                            self.showMessage(text)
            self.fileList = filelist
            if len(self.fileList):
                self.inputDir = os.path.dirname(self.fileList[0])
                PyMcaDirs.inputDir = os.path.dirname(self.fileList[0])
            if QTVERSION < '4.0.0':
                self.__listView.setText(text)
            else:
                self.__listView.clear()
                self.__listView.insertPlainText(text)

    def showMessage(self, text):
        msg = qt.QMessageBox(self)
        msg.setWindowTitle("PyMcaBatch Message")
        msg.setIcon(qt.QMessageBox.Information)
        msg.setText(text)
        msg.exec_()
        
    def setConfigFile(self,configfile=None):
        if configfile is None:return
        if self.__goodConfigFile(configfile):
            self.configFile = configfile
            if type(configfile) == type([]):
                #do not sort file list
                #self.configFile.sort()
                self.__configLine.setText(self.configFile[0])
                self.lastInputDir = os.path.dirname(self.configFile[0])
            else:
                self.__configLine.setText(configfile)
                self.lastInputDir = os.path.dirname(self.configFile)
        
    def setOutputDir(self,outputdir=None):
        if outputdir is None:
            return
        if self.__goodOutputDir(outputdir):
            self.outputDir = outputdir
            PyMcaDirs.outputDir = outputdir
            self.__outLine.setText(outputdir)
        else:
            qt.QMessageBox.critical(self, "ERROR",
                "Cannot use output directory:\n%s"% (outputdir))

    def __goodFileList(self,filelist):
        if not len(filelist):
            return True
        for ffile in filelist:
            if not os.path.exists(ffile):
                qt.QMessageBox.critical(self, "ERROR",
                                    'File %s\ndoes not exists' % ffile)
                if QTVERSION < '4.0.0':
                    self.raiseW()
                else:
                    self.raise_()
                return False
        return True
        
    def __goodConfigFile(self,configfile0):
        if type(configfile0) != type([]):
            configfileList = [configfile0]
        else:
            configfileList = configfile0
        for configfile in configfileList:
            if not os.path.exists(configfile):
                qt.QMessageBox.critical(self,
                             "ERROR",'File %s\ndoes not exists' % configfile)
                if QTVERSION < '4.0.0':
                    self.raiseW()
                else:
                    self.raise_()
                return False
            elif len(configfile.split()) > 1:
                if sys.platform != 'win32':
                    qt.QMessageBox.critical(self,
                                 "ERROR",'Configuration File:\n %s\ncontains spaces in the path' % configfile)
                    if QTVERSION < '4.0.0':
                        self.raiseW()
                    else:
                        self.raise_()
                    return False
        return True

    def __goodOutputDir(self,outputdir):
        if not os.path.isdir(outputdir):
            return False
        elif len(outputdir.split()) > 1:
            if sys.platform != 'win32':
                qt.QMessageBox.critical(self,
                    "ERROR",
                    'Output Directory:\n %s\ncontains spaces in the path' % outputdir)
                if QTVERSION < '4.0.0':
                    self.raiseW()
                else:
                    self.raise_()
                return False
        return True

    def __getFileType(self,inputfile):
        try:
            ffile = None
            try:
                ffile   = EdfFileLayer.EdfFileLayer(fastedf=0)
                ffile.SetSource(inputfile)
                fileinfo = ffile.GetSourceInfo()
                if fileinfo['KeyList'] == []:ffile=None
                return "EdfFile"
            except:
                pass
            if (ffile is None):
                ffile   = SpecFileLayer.SpecFileLayer()
                ffile.SetSource(inputfile)
            del ffile
            return "Specfile" 
        except:
            qt.QMessageBox.critical(self,
                                    sys.exc_info()[0],
                                    'I do not know what to do with file\n %s' % ffile)
            if QTVERSION < '4.0.0':
                self.raiseW()
            else:
                self.raise_()
            return None

    def browseList(self):
        self.inputDir = PyMcaDirs.inputDir
        if not os.path.exists(self.inputDir):
            self.inputDir =  os.getcwd()
        wdir = self.inputDir
        if QTVERSION < '4.0.0':
            filedialog = qt.QFileDialog(self,"Open a set of files",1)
            filedialog.setMode(filedialog.ExistingFiles)
            filedialog.setDir(wdir)
        else:
            filedialog = qt.QFileDialog(self)
            filedialog.setWindowTitle("Open a set of files")
            filedialog.setDirectory(wdir)
            filedialog.setModal(1)
            filedialog.setFileMode(filedialog.ExistingFiles)
        
        filetypes  = "McaFiles (*.mca)\nEdfFiles (*.edf)\n"
        if HDF5SUPPORT:
            filetypes += "HDF5 (*.nxs *.h5 *.hdf)\n"
        filetypes += "SpecFiles (*.spec)\nSpecFiles (*.dat)\nAll files (*)"
        #if (QTVERSION < '4.3.0') and sys.platform == "win32":
        if False and PyMcaDirs.nativeFileDialogs:
            #windows file dialogs have difficulties when dealing with
            #thousands of files, I do not know if other platforms too
            if QTVERSION < '4.0.0':
                filelist= filedialog.getOpenFileNames(qt.QString(filetypes),
                        wdir,
                        self,"openFile", "Open a set of files")
            else:
                filelist = qt.QFileDialog.getOpenFileNames(self,
                                "Open a set of files",
                                wdir,
                                filetypes,
                                None)    #This should be the last file filter used
        else:
            if QTVERSION < '4.0.0':
                filedialog.setFilters(filetypes)
                ret = filedialog.exec_loop()
            else:
                filedialog.setFilters(filetypes.split("\n"))
                ret = filedialog.exec_()
            if  ret == qt.QDialog.Accepted:
                filelist=filedialog.selectedFiles()
            else:
                if QTVERSION < '4.0.0':
                    self.raiseW()
                else:
                    self.raise_()
                return
        if len(filelist):
            filelist = [qt.safe_str(x) for x in filelist]
            self.setFileList(filelist)
        if QTVERSION < '4.0.0':
            self.raiseW()
        else:
            self.raise_()

    def browseConfig(self):
        self.inputDir = PyMcaDirs.inputDir
        if not os.path.exists(self.inputDir):
            self.inputDir =  os.getcwd()
        wdir = self.inputDir
        if QTVERSION < '4.0.0':
            filename = qt.QFileDialog(self,"Open a new fit config file",1)
            filename.setMode(filename.ExistingFiles)
            filename.setDir(wdir)
        else:
            filename = qt.QFileDialog(self)
            filename.setWindowTitle("Open a new fit config file")
            filename.setModal(1)
            filename.setFileMode(filename.ExistingFiles)
            filename.setDirectory(wdir)
        filetypes = "Config Files (*.cfg)\nAll files (*)"
        if PyMcaDirs.nativeFileDialogs:
            if QTVERSION < '4.0.0':
                filenameList= filename.getOpenFileNames(qt.QString(filetypes),
                            wdir,
                            self,"openFile", "Open a new fit config file")
            else:
                filenameList = qt.QFileDialog.getOpenFileNames(self,
                                    "Open a new fit config file",
                                    wdir,
                                    filetypes,
                                    None)    #This should be the filter used
        else:
            if QTVERSION < '4.0.0':
                filename.setFilters(filetypes)
                ret = filename.exec_loop() 
            else:
                filename.setFilters(["Config Files (*.cfg)", "All files (*)"])
                ret = filename.exec_()
                
            if  ret == qt.QDialog.Accepted:
                filenameList = filename.selectedFiles()
            else:
                if QTVERSION < '4.0.0':
                    self.raiseW()
                else:
                    self.raise_()
                return

        filename = []
        for f in filenameList:
            filename.append(qt.safe_str(f))

        if len(filename) == 1:
            self.setConfigFile(qt.safe_str(filename[0]))
        elif len(filenameList):
            self.setConfigFile(filename)
        if QTVERSION < '4.0.0':
            self.raiseW()
        else:
            self.raise_()

    def browseOutputDir(self):
        self.outputDir = PyMcaDirs.outputDir
        if not os.path.exists(self.outputDir):
            self.outputDir =  os.getcwd()
        wdir = self.outputDir
        if QTVERSION < '4.0.0':
            outfile = qt.QFileDialog(self,"Output Directory Selection",1)
            outfile.setMode(outfile.DirectoryOnly)
            outfile.setDir(wdir)
            ret = outfile.exec_loop()
        else:
            outfile = qt.QFileDialog(self)
            outfile.setWindowTitle("Output Directory Selection")
            outfile.setModal(1)
            outfile.setDirectory(wdir)
            outfile.setFileMode(outfile.DirectoryOnly)
            ret = outfile.exec_()
        if ret:
            if QTVERSION < '4.0.0':
                outdir=qt.safe_str(outfile.selectedFile())
            else:
                outdir=qt.safe_str(outfile.selectedFiles()[0])
            outfile.close()
            del outfile
            self.setOutputDir(outdir)
        else:
            outfile.close()
            del outfile
        if QTVERSION < '4.0.0':
            self.raiseW()
        else:
            self.raise_()
            
    def start(self):
        if not len(self.fileList):
            qt.QMessageBox.critical(self, "ERROR",'Empty file list')
            if QTVERSION < '4.0.0':
                self.raiseW()
            else:
                self.raise_()
            return

        if self.__splitBox.isChecked():
            if sys.platform == 'darwin':
                if ".app" in os.path.dirname(__file__):
                    text  = 'Multiple processes only supported on MacOS X when built from source\n'
                    text += 'and not when running the frozen binary.'  
                    qt.QMessageBox.critical(self, "ERROR",text)
                    if QTVERSION < '4.0.0':
                        self.raiseW()
                    else:
                        self.raise_()
                    return
            if len(self.fileList) == 1:
                if int(qt.safe_str(self.__splitSpin.text())) > 1:
                    allowSingleFileSplitProcesses = False
                    if HDF5SUPPORT:
                        if h5py.is_hdf5(self.fileList[0]):
                            if DEBUG:
                                print("Disallow single HDF5 process split")
                                print("due to problems with concurrent access")
                            #allowSingleFileSplitProcesses = True
                            allowSingleFileSplitProcesses = False
                    if not allowSingleFileSplitProcesses:
                        text = "Multiple processes can only be used with multiple input files."
                        qt.QMessageBox.critical(self, "ERROR",text)
                        if QTVERSION < '4.0.0':
                            self.raiseW()
                        else:
                            self.raise_()
                        return
                    
        if (self.configFile is None) or (not self.__goodConfigFile(self.configFile)):
            qt.QMessageBox.critical(self, "ERROR",'Invalid fit configuration file')
            if QTVERSION < '4.0.0':
                self.raiseW()
            else:
                self.raise_()
            return
        if type(self.configFile) == type([]):
            if len(self.configFile) != len(self.fileList):
                qt.QMessageBox.critical(self, "ERROR",
      'Number of config files should be either one or equal to number of files')
                if QTVERSION < '4.0.0':
                    self.raiseW()
                else:
                    self.raise_()
                return    
        if (self.outputDir is None) or (not self.__goodOutputDir(self.outputDir)):
            qt.QMessageBox.critical(self, "ERROR",'Invalid output directory')
            if QTVERSION < '4.0.0':
                self.raiseW()
            else:
                self.raise_()
            return
        name = "Batch from %s to %s " % (os.path.basename(self.fileList[ 0]),
                                          os.path.basename(self.fileList[-1]))
        roifit  = self.__roiBox.isChecked()
        html    = self.__htmlBox.isChecked()
        #if html:
        concentrations = self.__concentrationsBox.isChecked()
        #else:
        #    concentrations = 0
        if self.__tableBox.isChecked():
            if self.__extendedTable.isChecked():
                table = 2
            else:
                table = 1
        else:   table =0
        
        #htmlindex = qt.safe_str(self.__htmlIndex.text())
        htmlindex = "index.html"
        if html:
            if  len(htmlindex)<5:
                htmlindex+=".html"
            if  len(htmlindex) == 5:
                htmlindex = "index.html" 
            if htmlindex[-5:] != "html":
                htmlindex+=".html"
        roiwidth = float(qt.safe_str(self.__roiSpin.text()))
        if 0:
            overwrite= self.__overwrite.isChecked()
            filestep = int(qt.safe_str(self.__fileSpin.text()))
            mcastep  = int(qt.safe_str(self.__mcaSpin.text()))
        else:
            overwrite= 1
            filestep = 1
            mcastep = 1
            if len(self.fileList) == 1:
                if self.__splitBox.isChecked():
                    nbatches = int(qt.safe_str(self.__splitSpin.text()))
                    mcastep = nbatches
                
        fitfiles = self.__fitBox.isChecked()
        selection = self._selection
        if selection is None:
            selectionFlag = False
        else:
            selectionFlag = True

        if self._edfSimpleViewer is not None:
            self._edfSimpleViewer.close()
            self._edfSimpleViewer = None
        if roifit:
            window =  McaBatchWindow(name="ROI"+name,actions=1, outputdir=self.outputDir,
                                     html=html, htmlindex=htmlindex, table = 0)
            b = McaBatch(window,self.configFile,self.fileList,self.outputDir,roifit=roifit,
                         roiwidth=roiwidth,overwrite=overwrite,filestep=1,mcastep=1,
                         concentrations=0, fitfiles=fitfiles, selection=selection)
            def cleanup():
                b.pleasePause = 0
                b.pleaseBreak = 1
                b.wait()
                qt.qApp.processEvents()

            def pause():
                if b.pleasePause:
                    b.pleasePause=0
                    window.pauseButton.setText("Pause")
                else:
                    b.pleasePause=1
                    window.pauseButton.setText("Continue") 
            qt.QObject.connect(window.pauseButton,qt.SIGNAL("clicked()"),pause)
            qt.QObject.connect(window.abortButton,qt.SIGNAL("clicked()"),window.close)
            qt.QObject.connect(qt.qApp,qt.SIGNAL("aboutToQuit()"),cleanup)
            self.__window = window
            self.__b      = b
            window.show()
            b.start()
        elif (sys.platform == 'darwin') and\
             ((".app" in os.path.dirname(__file__)) or (not self.__splitBox.isChecked())):
            #almost identical to batch
            window =  McaBatchWindow(name="ROI"+name,actions=1,outputdir=self.outputDir,
                                     html=html,htmlindex=htmlindex, table = table)
            b = McaBatch(window,self.configFile,self.fileList,self.outputDir,roifit=roifit,
                         roiwidth=roiwidth,overwrite=overwrite,filestep=filestep,
                         mcastep=mcastep, concentrations=concentrations, fitfiles=fitfiles, selection=selection)
            def cleanup():
                b.pleasePause = 0
                b.pleaseBreak = 1
                b.wait()
                qt.qApp.processEvents()

            def pause():
                if b.pleasePause:
                    b.pleasePause=0
                    window.pauseButton.setText("Pause")
                else:
                    b.pleasePause=1
                    window.pauseButton.setText("Continue") 
            qt.QObject.connect(window.pauseButton,qt.SIGNAL("clicked()"),pause)
            qt.QObject.connect(window.abortButton,qt.SIGNAL("clicked()"),window.close)
            qt.QObject.connect(qt.qApp,qt.SIGNAL("aboutToQuit()"),cleanup)
            window._rootname = "%s"% b._rootname
            self.__window = window
            self.__b      = b
            window.show()
            b.start()
        elif sys.platform == 'win32':
            listfile = os.path.join(self.outputDir, "tmpfile")
            self.genListFile(listfile, config=False)
            dirname = os.path.dirname(McaAdvancedFitBatch.__file__)
            tmpDirname = os.path.dirname(dirname)
            if tmpDirname.lower().endswith("exe") or\
               tmpDirname.lower().endswith(".zip"):
                frozen = True
                dirname  = os.path.dirname(tmpDirname)
                myself   = os.path.join(dirname, "PyMcaBatch.exe")
                viewer   = os.path.join(dirname, "EdfFileSimpleViewer.exe")
                rgb    = os.path.join(dirname, "PyMcaPostBatch.exe")
                if not os.path.exists(viewer):
                    viewer = None
                if not os.path.exists(rgb):
                    rgb = None
            else:
                frozen = False
                myself = os.path.join(dirname, "PyMcaBatch.py")
                if not os.path.exists(myself):
                    dirname = os.path.dirname(McaAdvancedFitBatch.__file__)
                myself = os.path.join(dirname, "PyMcaBatch.py")
                viewer = os.path.join(dirname, "EdfFileSimpleViewer.py")
                rgb    = os.path.join(dirname, "PyMcaPostBatch.py")
                if not os.path.exists(viewer):
                    viewer = None
                else:
                    viewer = '%s "%s"' % (sys.executable, viewer)
                if not os.path.exists(rgb):
                    rgb    = None
                else:
                    rgb = '%s "%s"' % (sys.executable, rgb)
            self._rgb = rgb
            if type(self.configFile) == type([]):
                cfglistfile = os.path.join(self.outputDir, "tmpfile.cfg")
                self.genListFile(cfglistfile, config=True)
                dirname  = os.path.dirname(dirname)
                if frozen:
                    cmd = '"%s" --cfglistfile="%s" --outdir="%s" --overwrite=%d --filestep=%d --mcastep=%d --html=%d --htmlindex=%s --listfile="%s" --concentrations=%d --table=%d --fitfiles=%d --selection=%d' %\
                                                                  (myself,
                                                                    cfglistfile,
                                                                    self.outputDir, overwrite,
                                                                    filestep, mcastep,
                                                                    html,htmlindex,
                                                                    listfile,concentrations,
                                                                    table, fitfiles, selectionFlag)
                else:
                    cmd = '%s "%s" --cfglistfile="%s" --outdir="%s" --overwrite=%d --filestep=%d --mcastep=%d --html=%d --htmlindex=%s --listfile="%s" --concentrations=%d --table=%d --fitfiles=%d --selection=%d' %\
                                                                  (sys.executable,myself,
                                                                    cfglistfile,
                                                                    self.outputDir, overwrite,
                                                                    filestep, mcastep,
                                                                    html,htmlindex,
                                                                    listfile,concentrations,
                                                                    table, fitfiles, selectionFlag)
            else:
                if not frozen:
                    cmd = '%s "%s" --cfg="%s" --outdir="%s" --overwrite=%d --filestep=%d --mcastep=%d --html=%d --htmlindex=%s --listfile="%s" --concentrations=%d --table=%d --fitfiles=%d --selection=%d' % \
                                                                  (sys.executable, myself,
                                                                    self.configFile,
                                                                    self.outputDir, overwrite,
                                                                    filestep, mcastep,
                                                                    html,htmlindex,
                                                                    listfile,concentrations,
                                                                    table, fitfiles, selectionFlag)
                else:
                    cmd = '"%s" --cfg="%s" --outdir="%s" --overwrite=%d --filestep=%d --mcastep=%d --html=%d --htmlindex=%s --listfile="%s" --concentrations=%d --table=%d --fitfiles=%d --selection=%d' % \
                                                                  (myself,
                                                                    self.configFile,
                                                                    self.outputDir, overwrite,
                                                                    filestep, mcastep,
                                                                    html,htmlindex,
                                                                    listfile,concentrations,
                                                                    table, fitfiles, selectionFlag)
            self.hide()
            qt.qApp.processEvents()
            if DEBUG:
                print("cmd = %s" % cmd)
            if self.__splitBox.isChecked():
                nbatches = int(qt.safe_str(self.__splitSpin.text()))
                if len(self.fileList) > 1:
                    filechunk = int(len(self.fileList)/nbatches)
                    processList = []
                    for i in range(nbatches):
                        beginoffset = filechunk * i
                        endoffset   = len(self.fileList) - filechunk * (i+1)
                        if (i+1) == nbatches:endoffset = 0
                        cmd1        = cmd + " --filebeginoffset=%d --fileendoffset=%d --chunk=%d" % \
                                              (beginoffset, endoffset, i)
                        try:
                            processList.append(subprocess.Popen(cmd1,
                                                            cwd=os.getcwd()))
                        except UnicodeEncodeError:
                            processList.append(\
                                subprocess.Popen(cmd1.encode(sys.getfilesystemencoding()),
                                                 cwd=os.getcwd()))
                            
                        if DEBUG:
                            print("cmd = %s" % cmd1)
                else:
                    #f = open("CMD", 'wb')                        
                    processList = []
                    for i in range(nbatches):
                        #the mcastep has been dealt with above
                        cmd1        = cmd + " --mcaoffset=%d --chunk=%d" % (i, i)
                        processList.append(subprocess.Popen(cmd1, cwd=os.getcwd()))
                        if DEBUG:
                            print("CMD = %s" % cmd1)
                        #f.write(cmd1+"\n")
                    #f.close()
                self._processList = processList
                if self._timer is None:
                    self._timer = qt.QTimer(self)
                    qt.QObject.connect(self._timer,
                                       qt.SIGNAL('timeout()'),
                                       self._pollProcessList)
                if not self._timer.isActive():
                    self._timer.start(1000)
                else:
                    print("timer was already active")
                return
            else:
                try:
                    subprocess.call(cmd)
                except UnicodeEncodeError:
                    try:
                        subprocess.call(cmd.encode(sys.getfilesystemencoding()))
                    except:
                        # be ready for any weird error like missing that encoding
                        try:
                            subprocess.call(cmd.encode('utf-8'))
                        except UnicodeEncodeError:
                            subprocess.call(cmd.encode('latin-1'))
            self.show()
        else:
            listfile = os.path.join(self.outputDir, "tmpfile")
            self.genListFile(listfile, config=False)
            dirname = os.path.dirname(McaAdvancedFitBatch.__file__)
            tmpDirname = os.path.dirname(dirname)
            if tmpDirname.lower().endswith("exe") or\
               tmpDirname.lower().endswith(".zip"):
                frozen = True
                dirname  = os.path.dirname(tmpDirname)
                myself   = os.path.join(dirname, "PyMcaBatch") 
                viewer   = os.path.join(dirname, "EdfFileSimpleViewer")
                rgb    = os.path.join(dirname, "PyMcaPostBatch")
                if not os.path.exists(viewer):
                    viewer = None
                if not os.path.exists(rgb):
                    rgb = None
            else:
                myself = os.path.join(dirname, "PyMcaBatch.py")
                if not os.path.exists(os.path.join(dirname, myself)):
                    dirname = os.path.dirname(McaAdvancedFitBatch.__file__)
                    if not os.path.exists(os.path.join(dirname, myself)):
                        text  = 'Cannot locate PyMcaBatch.py file.\n'
                        qt.QMessageBox.critical(self, "ERROR",text)
                        if QTVERSION < '4.0.0':
                            self.raiseW()
                        else:
                            self.raise_()
                        return
                myself  = sys.executable+" "+ os.path.join(dirname, myself)
                viewer = os.path.join(dirname, "EdfFileSimpleViewer.py")
                rgb    = os.path.join(dirname, "PyMcaPostBatch.py")
                if not os.path.exists(viewer):
                    viewer = None
                else:
                    viewer = '%s "%s"' % (sys.executable, viewer)
                if not os.path.exists(rgb):
                    rgb    = None
                else:
                    rgb = '%s "%s"' % (sys.executable, rgb)

            self._rgb = rgb
            if type(self.configFile) == type([]):
                cfglistfile = os.path.join(self.outputDir, "tmpfile.cfg")
                self.genListFile(cfglistfile, config=True)
                cmd = "%s --cfglistfile=%s --outdir=%s --overwrite=%d --filestep=%d --mcastep=%d --html=%d --htmlindex=%s --listfile=%s  --concentrations=%d --table=%d --fitfiles=%d --selection=%d &" % \
                                                    (myself,
                                                    cfglistfile,
                                                    self.outputDir, overwrite,
                                                    filestep, mcastep, html, htmlindex, listfile,
                                                    concentrations, table, fitfiles, selectionFlag)
            else:
                cmd = "%s --cfg=%s --outdir=%s --overwrite=%d --filestep=%d --mcastep=%d --html=%d --htmlindex=%s --listfile=%s  --concentrations=%d --table=%d --fitfiles=%d  --selection=%d &" % \
                                                   (myself, self.configFile,
                                                    self.outputDir, overwrite,
                                                    filestep, mcastep, html, htmlindex,
                                                    listfile, concentrations, table, fitfiles, selectionFlag)
            if DEBUG:
                print("cmd = %s" % cmd)
            if self.__splitBox.isChecked():
                qt.qApp.processEvents()
                nbatches = int(qt.safe_str(self.__splitSpin.text()))
                filechunk = int(len(self.fileList)/nbatches)
                processList = []
                for i in range(nbatches):
                    beginoffset = filechunk * i
                    endoffset   = len(self.fileList) - filechunk * (i+1)
                    if (i+1) == nbatches:endoffset = 0
                    cmd1 = cmd.replace("&", "") + \
                           " --filebeginoffset=%d --fileendoffset=%d --chunk=%d" % \
                            (beginoffset, endoffset, i)
                    # unfortunately I have to set shell = True
                    # otherways I get a file not found error in the
                    # chold process
                    processList.append(\
                        subprocess.Popen(cmd1,
                                         cwd=os.getcwd(),
                                         shell=True,
                                         close_fds=True))
                self._processList = processList
                self.hide()
                self._pollProcessList()
                if self._timer is None:
                    self._timer = qt.QTimer(self)
                    qt.QObject.connect(self._timer,
                                       qt.SIGNAL('timeout()'),
                                       self._pollProcessList)
                if not self._timer.isActive():
                    self._timer.start(1000)
                else:
                    print("timer was already active")
                return
            else:
                os.system(cmd)
                msg = qt.QMessageBox(self)
                msg.setIcon(qt.QMessageBox.Information)
                text = "Your batch has been started as an independent process."
                msg.setText(text)
                if QTVERSION < '4.0.0':
                    msg.exec_loop()
                else:
                    msg.exec_()
            
    def genListFile(self, listfile, config=None):
        if os.path.exists(listfile):
            try:
                os.remove(listfile)
            except:
                print("Cannot delete file %s" % listfile)
                raise
        if config is None:
            lst = self.fileList
        elif config:
            lst = self.configFile
        else:
            lst = self.fileList

        if lst == self.fileList:
            if self._selection is not None:
                ddict = ConfigDict.ConfigDict()
                ddict['PyMcaBatch'] = {}
                ddict['PyMcaBatch']['filelist'] = lst
                ddict['PyMcaBatch']['selection'] = self._selection
                ddict.write(listfile)
                return
        fd=open(listfile, 'wb')
        for filename in lst:
            # only the file system encoding makes sense here
            fd.write(('%s\n' % filename).encode(sys.getfilesystemencoding()))
        fd.close()

    def _pollProcessList(self):
        processList = self._processList
        rgb         = self._rgb
        if QTVERSION < '4.0.0':rgb = None
        n = 0
        for process in processList:
            if process.poll() is None:
                n += 1
        if n > 0: return
        self._timer.stop()
        self.show()
        if QTVERSION < '4.0.0':
            self.raiseW()
        else:
            self.raise_()

        work = PyMcaBatchBuildOutput.PyMcaBatchBuildOutput(os.path.join(self.outputDir, "IMAGES"))
        if DEBUG:
            a, b, c = work.buildOutput(delete=False)
        else:
            a, b, c = work.buildOutput(delete=True)
        if len(a):
            if self._edfSimpleViewer is None:
                self._edfSimpleViewer = EdfFileSimpleViewer.EdfFileSimpleViewer()
            self._edfSimpleViewer.setFileList(a)
            self._edfSimpleViewer.show()
        if rgb is not None:
            if len(b):
                if sys.platform == "win32":
                    try:
                        subprocess.Popen('%s "%s"' % (rgb, b[0]),
                                         cwd = os.getcwd())
                    except UnicodeEncodeError:
                        subprocess.Popen(('%s "%s"' % (rgb, b[0])).encode(sys.getfilesystemencoding()),
                                         cwd = os.getcwd())
                else:
                    os.system("%s %s &" % (rgb, b[0]))
        work = PyMcaBatchBuildOutput.PyMcaBatchBuildOutput(self.outputDir)
        if DEBUG:
            work.buildOutput(delete=False)
        else:
            work.buildOutput(delete=True)


class McaBatch(qt.QThread,McaAdvancedFitBatch.McaAdvancedFitBatch):
    def __init__(self, parent, configfile, filelist=None, outputdir = None,
                     roifit = None, roiwidth=None, overwrite=1,
                     filestep=1, mcastep=1, concentrations=0,
                     fitfiles=0, filebeginoffset=0, fileendoffset=0,
                     mcaoffset=0, chunk=None,
                     selection=None, lock=None):
        McaAdvancedFitBatch.McaAdvancedFitBatch.__init__(self, configfile, filelist, outputdir,
                                                         roifit=roifit, roiwidth=roiwidth,
                                                         overwrite=overwrite, filestep=filestep,
                                                         mcastep=mcastep,
                                                         concentrations=concentrations,
                                                         fitfiles=fitfiles,
                                                         filebeginoffset = filebeginoffset,
                                                         fileendoffset = fileendoffset,
                                                         mcaoffset  = mcaoffset,
                                                         chunk=chunk,
                                                         selection=selection,
                                                         lock=lock) 
        qt.QThread.__init__(self)
        self.parent = parent
        self.pleasePause = 0
        
    def run(self):
        self.processList()

    def onNewFile(self, ffile, filelist):
        self.__lastOnNewFile = ffile
        ddict = {'file':ffile,
                 'filelist':filelist,
                 'filestep':self.fileStep,
                 'filebeginoffset':self.fileBeginOffset,
                 'fileendoffset':self.fileEndOffset,
                 'event':'onNewFile'}
        
        if QTVERSION < '4.0.0':
            self.postEvent(self.parent, McaCustomEvent.McaCustomEvent(ddict))
        else:
            qt.QApplication.postEvent(self.parent, McaCustomEvent.McaCustomEvent(ddict))

        if self.pleasePause:self.__pauseMethod()

    def onImage(self, key, keylist):
        ddict = {'key':key, 'keylist':keylist, 'event':'onImage'}
        if QTVERSION < '4.0.0':
            self.postEvent(self.parent, McaCustomEvent.McaCustomEvent(ddict))
        else:
            qt.QApplication.postEvent(self.parent, McaCustomEvent.McaCustomEvent(ddict))

    def onMca(self, mca, nmca, filename=None, key=None, info=None):
        if DEBUG:
            print("onMca key = %s" % key)
        ddict = {'mca':mca,
                 'nmca':nmca,
                 'mcastep':self.mcaStep,
                 'filename':filename,
                 'key':key,
                 'info':info,
                 'outputdir':self._outputdir,
                 'useExistingFiles':self.useExistingFiles,
                 'roifit':self.roiFit,
                 'event':'onMca'}
        
        if QTVERSION < '4.0.0':
            self.postEvent(self.parent, McaCustomEvent.McaCustomEvent(ddict))
        else:
            qt.QApplication.postEvent(self.parent, McaCustomEvent.McaCustomEvent(ddict))
            
        if self.pleasePause:self.__pauseMethod()
                                                                   
    def onEnd(self):
        if DEBUG:
            print("onEnd")

        ddict = {'event':'onEnd',
                 'filestep':self.fileStep,
                 'mcastep':self.mcaStep,
                 'chunk':self.chunk,
                 'savedimages':self.savedImages}

        if QTVERSION < '4.0.0':
            self.postEvent(self.parent, McaCustomEvent.McaCustomEvent(ddict))
        else:
            qt.QApplication.postEvent(self.parent, McaCustomEvent.McaCustomEvent(ddict))

        if self.pleasePause:
            self.__pauseMethod()
        
    def __pauseMethod(self):
        if QTVERSION < '4.0.0':
            self.postEvent(self.parent, McaCustomEvent.McaCustomEvent({'event':'batchPaused'}))
        else:
            qt.QApplication.postEvent(self.parent, McaCustomEvent.McaCustomEvent({'event':'batchPaused'}))

        while(self.pleasePause):
            time.sleep(1)

        if QTVERSION < '4.0.0':
            self.postEvent(self.parent, McaCustomEvent.McaCustomEvent({'event':'batchResumed'}))
        else:
            qt.QApplication.postEvent(self.parent, McaCustomEvent.McaCustomEvent({'event':'batchResumed'}))
            

class McaBatchWindow(qt.QWidget):
    def __init__(self,parent=None, name="BatchWindow", fl=0, actions = 0, outputdir=None, html=0,
                    htmlindex = None, table=2, chunk = None):
        if QTVERSION < '4.0.0':
            qt.QWidget.__init__(self, parent, name, fl)
            self.setCaption(name)
        else:
            qt.QWidget.__init__(self, parent)
            self.setWindowTitle(name)
        self.chunk = chunk
        self.l = qt.QVBoxLayout(self)
        #self.l.setAutoAdd(1)
        self.bars =qt.QWidget(self)
        self.l.addWidget(self.bars)
        if QTVERSION < '4.0.0':
            self.barsLayout = qt.QGridLayout(self.bars,2,3)
        else:
            self.barsLayout = qt.QGridLayout(self.bars)
            self.barsLayout.setMargin(2)
            self.barsLayout.setSpacing(3)
        self.progressBar   = qt.QProgressBar(self.bars)
        self.progressLabel = qt.QLabel(self.bars)
        self.progressLabel.setText('File Progress:')
        self.imageBar   = qt.QProgressBar(self.bars)
        self.imageLabel = qt.QLabel(self.bars)
        self.imageLabel.setText('Image in File:')
        self.mcaBar   = qt.QProgressBar(self.bars)
        self.mcaLabel = qt.QLabel(self.bars)
        self.mcaLabel.setText('MCA in Image:')
        
        self.barsLayout.addWidget(self.progressLabel,0,0)        
        self.barsLayout.addWidget(self.progressBar,0,1)
        self.barsLayout.addWidget(self.imageLabel,1,0)        
        self.barsLayout.addWidget(self.imageBar,1,1)
        self.barsLayout.addWidget(self.mcaLabel,2,0)        
        self.barsLayout.addWidget(self.mcaBar,2,1)
        self.status      = qt.QLabel(self)
        self.status.setText(" ")
        self.timeLeft      = qt.QLabel(self)
        self.l.addWidget(self.status)
        self.l.addWidget(self.timeLeft)
        
        self.timeLeft.setText("Estimated time left = ???? min")
        self.time0 = None
        self.html  = html
        if htmlindex is None:htmlindex="index.html"
        self.htmlindex = htmlindex
        self.outputdir = outputdir
        self.table = table
        self.__ended         = False
        self.__writingReport = False

        if actions: self.addButtons()
        self.show()
        if QTVERSION < '4.0.0':
            self.raiseW()
        else:
            self.raise_()


    def addButtons(self):
        self.actions = 1
        self.buttonsBox = qt.QWidget(self)
        l = qt.QHBoxLayout(self.buttonsBox)
        l.addWidget(qt.HorizontalSpacer(self.buttonsBox))        
        self.pauseButton = qt.QPushButton(self.buttonsBox)
        l.addWidget(self.pauseButton)
        l.addWidget(qt.HorizontalSpacer(self.buttonsBox))
        self.pauseButton.setText("Pause")
        self.abortButton   = qt.QPushButton(self.buttonsBox)
        l.addWidget(self.abortButton)
        l.addWidget(qt.HorizontalSpacer(self.buttonsBox))
        self.abortButton.setText("Abort")
        self.l.addWidget(self.buttonsBox)
        self.update()

    def customEvent(self,event):
        if   event.dict['event'] == 'onNewFile':self.onNewFile(event.dict['file'],
                                                               event.dict['filelist'],
                                                               event.dict['filestep'],
                                                               event.dict['filebeginoffset'],
                                                               event.dict['fileendoffset'])
        elif event.dict['event'] == 'onImage':  self.onImage  (event.dict['key'],
                                                               event.dict['keylist'])
        elif event.dict['event'] == 'onMca':    self.onMca    (event.dict)
                                                               #event.dict['mca'],
                                                               #event.dict['nmca'],
                                                               #event.dict['mcastep'],
                                                               #event.dict['filename'],
                                                               #event.dict['key'])
        elif event.dict['event'] == 'onEnd':    self.onEnd(event.dict)

        elif event.dict['event'] == 'batchPaused': self.onPause()
        
        elif event.dict['event'] == 'batchResumed':self.onResume()

        elif event.dict['event'] == 'reportWritten':self.onReportWritten()

        else:
            print("Unhandled event %s" % event) 
                                                

    def onNewFile(self, file, filelist, filestep, filebeginoffset =0, fileendoffset = 0):
        if DEBUG:
            print("onNewFile: %s" % file)
        indexlist = list(range(0,len(filelist),filestep))
        index  = indexlist.index(filelist.index(file)) - filebeginoffset
        #print index + filebeginoffset
        if index == 0:
            self.report= None
            if self.html:
                self.htmlindex = os.path.join(self.outputdir, 'HTML')
                htmlindex = os.path.join(os.path.basename(file)+"_HTMLDIR", 
                            "index.html")
                self.htmlindex = os.path.join(self.htmlindex,htmlindex)
                if os.path.exists(self.htmlindex):
                    try:
                        os.remove(self.htmlindex)
                    except:
                        print("cannot delete file %s" % self.htmlindex)
        nfiles = len(indexlist)-filebeginoffset-fileendoffset
        self.status.setText("Processing file %s" % file)
        e = time.time()
        if QTVERSION < '4.0.0':
            self.progressBar.setTotalSteps(nfiles)
            self.progressBar.setProgress(index)
        else:
            self.progressBar.setMaximum(nfiles)
            self.progressBar.setValue(index)
        if self.time0 is not None:
            t = (e - self.time0) * (nfiles - index)
            self.time0 =e
            if t < 120:
                self.timeLeft.setText("Estimated time left = %d sec" % (t))
            else:
                self.timeLeft.setText("Estimated time left = %d min" % (int(t / 60.)))
        else:
            self.time0 = e
        if sys.platform == 'darwin':
            qt.qApp.processEvents()

    def onImage(self,key,keylist):
        if DEBUG:
            print("onImage %s" % key)
        i = keylist.index(key) + 1
        n = len(keylist)
        if QTVERSION < '4.0.0':
            self.imageBar.setTotalSteps(n)
            self.imageBar.setProgress(i)
            self.mcaBar.setTotalSteps(1)
            self.mcaBar.setProgress(0)
        else:
            self.imageBar.setMaximum(n)
            self.imageBar.setValue(i)
            self.mcaBar.setMaximum(1)
            self.mcaBar.setValue(0)
            

    #def onMca(self, mca, nmca, mcastep):
    def onMca(self, ddict):
        if DEBUG:
            print("onMca ", ddict['mca'])
        mca  = ddict['mca']
        nmca = ddict['nmca']
        mcastep  = ddict['mcastep']
        filename = ddict['filename']
        key = ddict['key']
        info = ddict['info']
        outputdir = ddict['outputdir']
        useExistingFiles = ddict['useExistingFiles']
        self.roiFit = ddict['roifit']
        if self.html:
            try:
                if not self.roiFit:
                    if mca == 0:
                        self.__htmlReport(filename, key, outputdir,
                                          useExistingFiles, info, firstmca = True)
                    else:
                        self.__htmlReport(filename, key, outputdir,
                                          useExistingFiles, info, firstmca = False)
            except:
                print("ERROR on REPORT",sys.exc_info())
                print(sys.exc_info()[1])
                print("filename = %s key =%s " % (filename, key))
                print("If your batch is stopped, please report this")
                print("error sending the above mentioned file and the")
                print("associated fit configuration file.")
        if QTVERSION < '4.0.0':
            self.mcaBar.setTotalSteps(nmca)
            self.mcaBar.setProgress(mca)
        else:
            self.mcaBar.setMaximum(nmca)
            self.mcaBar.setValue(mca)
        if sys.platform == 'darwin':
            qt.qApp.processEvents()

    def __htmlReport(self, filename, key, outputdir, useExistingFiles, info=None, firstmca = True): 
        """
        file=self.file
        fileinfo = file.GetSourceInfo()
        nimages = nscans = len(fileinfo['KeyList'])

        filename = os.path.basename(info['SourceName'])
        """
        fitdir = os.path.join(outputdir,"HTML")
        if not os.path.exists(fitdir):
            try:
                os.mkdir(fitdir)
            except:
                print("I could not create directory %s" % fitdir)
                return
        fitdir = os.path.join(fitdir,filename+"_HTMLDIR")
        if not os.path.exists(fitdir):
            try:
                os.mkdir(fitdir)
            except:
                print("I could not create directory %s" % fitdir)
                return
        localindex = os.path.join(fitdir, "index.html")
        if not os.path.isdir(fitdir):
            print("%s does not seem to be a valid directory" % fitdir)
        else:
            outfile = filename +"_"+key+".html" 
            outfile = os.path.join(fitdir,  outfile)
        useExistingResult = useExistingFiles
        if os.path.exists(outfile):
            if not useExistingFiles:
                try:
                    os.remove(outfile)
                except:
                    print("cannot delete file %s" % outfile)
                useExistingResult = 0
        else:
            useExistingResult = 0    
        outdir = fitdir
        fitdir = os.path.join(outputdir,"FIT")
        fitdir = os.path.join(fitdir,filename+"_FITDIR")
        fitfile= os.path.join(fitdir,  filename +"_"+key+".fit")
        if not os.path.exists(fitfile):
            print("fit file %s does not exists!" % fitfile)
            return
        if self.report is None:
            #first file
            self.forcereport = 0        
            self._concentrationsFile = os.path.join(outputdir,
                                self._rootname + "_concentrations.txt")
            if os.path.exists(self._concentrationsFile):
                """
                #code removed, concentrations in McaAdvancedFitBatch.py
                try:
                    os.remove(self._concentrationsFile)
                except:
                    pass
                """
                pass
            else:
                #this is to generate the concentrations file
                #from an already existing set of fitfiles
                self.forcereport = 1                        
        if self.forcereport or (not useExistingResult):
            self.report = QtMcaAdvancedFitReport.QtMcaAdvancedFitReport(fitfile = fitfile,
                        outfile = outfile, table = self.table)
            self.__writingReport = True
            a=self.report.writeReport()
            """
            #The code below has been moved to McaAdvancedFitBatch.py
            if len(self.report._concentrationsTextASCII) > 1:
                text  = ""
                text += "SOURCE: "+ filename +"\n"
                text += "KEY: "+key+"\n"
                text += self.report._concentrationsTextASCII + "\n"
                f=open(self._concentrationsFile,"a")
                f.write(text)
                f.close()
            """
            self.__writingReport = False
            #qt.QApplication.postEvent(self, McaCustomEvent.McaCustomEvent({'event':'reportWritten'}))
            self.onReportWritten()
            
    def onEnd(self,dict):
        self.__ended = True
        if QTVERSION < '4.0.0':
            n = self.progressBar.progress()
            self.progressBar.setProgress(n+dict['filestep'])
            n = self.mcaBar.progress()
            self.mcaBar.setProgress(n+dict['mcastep'])
        else:
            n = self.progressBar.value()
            self.progressBar.setValue(n+dict['filestep'])
            n = self.mcaBar.value()
            self.mcaBar.setValue(n+dict['mcastep'])
        self.status.setText  ("Batch Finished")
        self.timeLeft.setText("Estimated time left = 0 sec")
        if self.actions:
            self.pauseButton.hide()
            self.abortButton.setText("OK")
        if self.chunk is None:
            if 'savedimages' in dict:
                self.plotImages(dict['savedimages'])
        if self.html:
            if not self.__writingReport:
                directory = os.path.join(self.outputdir,"HTML")
                a = HtmlIndex.HtmlIndex(directory)
                a.buildRecursiveIndex()
        if dict['chunk'] is not None:
            if 0:
                #this was giving troubles using HDF5 files as input
                sys.exit(0)
            else:
                #this seems to work properly
                self.close()
    
    def onReportWritten(self):
        if self.__ended:
            directory = os.path.join(self.outputdir,"HTML")
            a = HtmlIndex.HtmlIndex(directory)
            a.buildRecursiveIndex()

    def onPause(self):    
        pass

    def onResume(self):    
        pass
        
        
    def plotImages(self,imagelist):
        if (sys.platform == 'win32') or (sys.platform == 'darwin'):
            filelist = " "
            for ffile in imagelist:
                filelist+=" %s" % ffile
            try:
                dirname = os.path.dirname(__file__)
            except:
                dirname = os.path.dirname(\
                        os.path.dirname(McaAdvancedFitBatch.__file__))
            if (dirname[-3:] == "exe") or\
               (dirname.lower().endswith(".zip")):
                myself  = os.path.dirname(dirname) 
                myself  = os.path.join(myself, "EdfFileSimpleViewer.exe")
            else:
                myself  = os.path.join(dirname, "EdfFileSimpleViewer.py")
            cmd = '"%s" %s ' % (myself, filelist)
            if 0:
                self.hide()
                qt.qApp.processEvents()
                os.system(cmd)
                self.show()                
            else:
                self.__viewer = EdfFileSimpleViewer.EdfFileSimpleViewer()
                self.__viewer.setFileList(imagelist)
                self.__viewer.show()
        else:
            filelist = " "
            for ffile in imagelist:
                filelist+=" %s" % ffile
            try:
                dirname = os.path.dirname(__file__)
            except:
                dirname = os.path.dirname(\
                        os.path.dirname(McaAdvancedFitBatch.__file__))
            if DEBUG:
                print("final dirname = %s" % dirname)
            if (dirname[-3:] == "exe") or\
               (dirname.lower().endswith(".zip")):
                myself  = os.path.dirname(dirname) 
                myself  = os.path.join(myself, "EdfFileSimpleViewer")
            else:
                myself  = sys.executable+" "+os.path.join(dirname, "EdfFileSimpleViewer.py")
            cmd = "%s %s &" % (myself, filelist)
            if DEBUG:
                print("cmd = %s" % cmd)
            os.system(cmd)
                              
def main():
    import getopt
    options     = 'f'
    longoptions = ['cfg=','outdir=','roifit=','roi=','roiwidth=',
                   'overwrite=', 'filestep=', 'mcastep=', 'html=','htmlindex=',
                   'listfile=','cfglistfile=', 'concentrations=', 'table=', 'fitfiles=',
                   'filebeginoffset=','fileendoffset=','mcaoffset=', 'chunk=',
                   'nativefiledialogs=','selection=']
    filelist = None
    outdir   = None
    cfg      = None
    listfile = None
    cfglistfile = None
    selection = False
    roifit   = 0
    roiwidth = ROIWIDTH
    overwrite= 1
    filestep = 1
    html = 0
    htmlindex= None
    mcastep  = 1
    table    = 2
    fitfiles = 1
    concentrations = 0
    filebeginoffset = 0
    fileendoffset = 0
    mcaoffset = 0
    chunk = None
    opts, args = getopt.getopt(
                    sys.argv[1:],
                    options,
                    longoptions)
    for opt,arg in opts:
        if opt in ('--cfg'):
            cfg = arg
        elif opt in ('--outdir'):
            outdir = arg
        elif opt in ('--roi','--roifit'):
            roifit   = int(arg)
        elif opt in ('--roiwidth'):
            roiwidth = float(arg)
        elif opt in ('--overwrite'):
            overwrite= int(arg)
        elif opt in ('--filestep'):
            filestep = int(arg)
        elif opt in ('--mcastep'):
            mcastep  = int(arg)
        elif opt in ('--html'):
            html  = int(arg)
        elif opt in ('--htmlindex'):
            htmlindex  = arg
        elif opt in ('--listfile'):
            listfile  = arg
        elif opt in ('--cfglistfile'):
            cfglistfile  = arg
        elif opt in ('--concentrations'):
            concentrations  = int(arg)
        elif opt in ('--table'):
            table  = int(arg)
        elif opt in ('--fitfiles'):
            fitfiles  = int(arg)
        elif opt in ('--filebeginoffset'):
            filebeginoffset = int(arg)
        elif opt in ('--fileendoffset'):
            fileendoffset   = int(arg)
        elif opt in ('--mcaoffset'):
            mcaoffset  = int(arg)
        elif opt in ('--chunk'):
            chunk  = int(arg)
        elif opt in ('--selection'):
            selection  = int(arg)
            if selection:
                selection = True
            else:
                selection = False
        elif opt in ('--nativefiledialogs'):
            if int(arg):
                PyMcaDirs.nativeFileDialogs = True
            else:
                PyMcaDirs.nativeFileDialogs = False
    if listfile is None: 
        filelist=[]
        for item in args:
            filelist.append(item)
        selection = None
    else:
        if selection:
            tmpDict = ConfigDict.ConfigDict()
            tmpDict.read(listfile)
            tmpDict = tmpDict['PyMcaBatch']
            filelist = tmpDict['filelist']
            if type(filelist) == type(""):
                filelist = [filelist]
            selection = tmpDict['selection']
        else:
            fd = open(listfile, 'rb')
            filelist = fd.readlines()
            fd.close()
            for i in range(len(filelist)):
                filelist[i]=filelist[i].decode(sys.getfilesystemencoding()).replace('\n','')
            selection = None
    if cfglistfile is not None:
        fd = open(cfglistfile, 'rb')
        cfg = fd.readlines()
        fd.close()
        for i in range(len(cfg)):
            cfg[i]=cfg[i].decode(sys.getfilesystemencoding()).replace('\n','')
    app=qt.QApplication(sys.argv) 
    winpalette = qt.QPalette(qt.QColor(230,240,249),qt.QColor(238,234,238))
    app.setPalette(winpalette)       
    if len(filelist) == 0:
        qt.QObject.connect(app,qt.SIGNAL("lastWindowClosed()"),app, qt.SLOT("quit()"))
        w = McaBatchGUI(actions=1)
        if QTVERSION < '4.0.0':
            app.setMainWidget(w)
            w.show()
            app.exec_loop()
        else:
            w.show()
            w.raise_()
            app.exec_()
    else:
        qt.QObject.connect(app,qt.SIGNAL("lastWindowClosed()"),app, qt.SLOT("quit()"))
        text = "Batch from %s to %s" % (os.path.basename(filelist[0]), os.path.basename(filelist[-1]))
        window =  McaBatchWindow(name=text,actions=1,
                                outputdir=outdir,html=html, htmlindex=htmlindex, table=table,
                                chunk = chunk)
                                
        if html:fitfiles=1
        try:
            b = McaBatch(window,cfg,filelist,outdir,roifit=roifit,roiwidth=roiwidth,
                     overwrite = overwrite, filestep=filestep, mcastep=mcastep,
                      concentrations=concentrations, fitfiles=fitfiles,
                      filebeginoffset=filebeginoffset,fileendoffset=fileendoffset,
                      mcaoffset=mcaoffset, chunk=chunk, selection=selection)
        except:
            msg = qt.QMessageBox()
            msg.setIcon(qt.QMessageBox.Critical)
            msg.setText("%s" % sys.exc_info()[1])
            if QTVERSION < '4.0.0':
                msg.exec_loop()
            else:
                msg.exec_()
            return

            
        def cleanup():
            b.pleasePause = 0
            b.pleaseBreak = 1
            b.wait()
            qt.qApp.processEvents()

        def pause():
            if b.pleasePause:
                b.pleasePause=0
                window.pauseButton.setText("Pause")
            else:
                b.pleasePause=1
                window.pauseButton.setText("Continue") 
        qt.QObject.connect(window.pauseButton,qt.SIGNAL("clicked()"),pause)
        qt.QObject.connect(window.abortButton,qt.SIGNAL("clicked()"),window.close)
        qt.QObject.connect(app,qt.SIGNAL("aboutToQuit()"),cleanup)        
        window._rootname = "%s"% b._rootname
        window.show()
        b.start()
        if QTVERSION < '4.0.0':
            app.setMainWidget(window)
            app.exec_loop()
        else:
            app.exec_()
 
if __name__ == "__main__":
    main()
 
 
# PyMcaBatch.py --cfg=/mntdirect/_bliss/users/sole/COTTE/WithLead.cfg --outdir=/tmp/   /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0007.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0008.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0009.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0010.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0011.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0012.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0013.edf &
# PyMcaBatch.exe --cfg=E:/COTTE/WithLead.cfg --outdir=C:/tmp/   E:/COTTE/ch09/ch09__mca_0003_0000_0007.edf E:/COTTE/ch09/ch09__mca_0003_0000_0008.edf
