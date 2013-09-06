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
__revision__ = "$Revision: 1.13 $"
import sys
import os
import numpy
import time

from PyMca import PyMcaQt as qt
QTVERSION = qt.qVersion()
if QTVERSION >= '4.0.0':
    qt.Qt.WDestructiveClose = "TO BE DONE"
from PyMca.PyMca_Icons import IconDict
from PyMca import McaCustomEvent
from PyMca import EdfFile
from PyMca import SpecFileLayer
from PyMca import PyMcaDirs


class Mca2EdfGUI(qt.QWidget):
    def __init__(self,parent=None,name="Mca to Edf Conversion",fl=qt.Qt.WDestructiveClose,
                filelist=None,outputdir=None, actions=0):
        if qt.qVersion() < '4.0.0':
            qt.QWidget.__init__(self,parent,name,fl)
            self.setIcon(qt.QPixmap(IconDict['gioconda16']))
            self.setCaption(name)
        else:
            qt.QWidget.__init__(self, parent)
            self.setWindowTitle(name)
            self.setWindowIcon(qt.QIcon(qt.QPixmap(IconDict['gioconda16'])))
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(0)
        #layout.setAutoAdd(1)
        self.__build(actions)               
        if filelist is None: filelist = []
        self.outputDir = None
        self.inputDir  = None
        self.setFileList(filelist)
        self.setOutputDir(outputdir)
    
    def __build(self,actions):
        self.__grid= qt.QWidget(self)
        #self.__grid.setGeometry(qt.QRect(30,30,288,156))
        if QTVERSION < '4.0.0':
            grid       = qt.QGridLayout(self.__grid,3,3,11,6)
            grid.setColStretch(0,0)
            grid.setColStretch(1,1)
            grid.setColStretch(2,0)
        else:
            grid  = qt.QGridLayout(self.__grid)
            grid.setMargin(11)
            grid.setSpacing(6)
        #input list
        listrow  = 0
        listlabel   = qt.QLabel(self.__grid)
        listlabel.setText("Input File list:")
        if QTVERSION < '4.0.0':
            listlabel.setAlignment(qt.QLabel.WordBreak | qt.QLabel.AlignVCenter)
            self.__listView   = qt.QTextView(self.__grid)
        else:
            self.__listView   = qt.QTextEdit(self.__grid)
        self.__listView.setMaximumHeight(30*listlabel.sizeHint().height())
        self.__listButton = qt.QPushButton(self.__grid)
        self.__listButton.setText('Browse')
        self.connect(self.__listButton,qt.SIGNAL('clicked()'),self.browseList) 
        grid.addWidget(listlabel,        listrow, 0, qt.Qt.AlignTop|qt.Qt.AlignLeft)
        grid.addWidget(self.__listView,  listrow, 1)
        grid.addWidget(self.__listButton,listrow, 2, qt.Qt.AlignTop|qt.Qt.AlignRight)
        
        #output dir
        outrow    = 1
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

        #step
        filesteprow =2
        filesteplabel = qt.QLabel(self.__grid)
        filesteplabel.setText("New EDF file each")
        filesteplabel2 = qt.QLabel(self.__grid)
        self.__fileSpin = qt.QSpinBox(self.__grid)
        if QTVERSION < '4.0.0':
            self.__fileSpin.setMinValue(1)
            self.__fileSpin.setMaxValue(999999)
        else:
            self.__fileSpin.setMinimum(1)
            self.__fileSpin.setMaximum(999999)
        self.__fileSpin.setValue(1)
        
        filesteplabel2.setText("mca")
        grid.addWidget(filesteplabel,  filesteprow, 0, qt.Qt.AlignLeft)
        grid.addWidget(self.__fileSpin,filesteprow, 1)
        grid.addWidget(filesteplabel2, filesteprow, 2, qt.Qt.AlignLeft)

        self.mainLayout.addWidget(self.__grid)
        if actions: self.__buildActions()

    def __buildActions(self):
        box = qt.QWidget(self)
        boxLayout = qt.QHBoxLayout(box)
        boxLayout.addWidget(qt.HorizontalSpacer(box))
        self.__dismissButton = qt.QPushButton(box)
        boxLayout.addWidget(self.__dismissButton)
        boxLayout.addWidget(qt.HorizontalSpacer(box))
        self.__dismissButton.setText("Close")
        self.__startButton   = qt.QPushButton(box)
        boxLayout.addWidget(self.__startButton)
        boxLayout.addWidget(qt.HorizontalSpacer(box))
        self.__startButton.setText("Start")
        self.mainLayout.addWidget(box)
        self.connect(self.__dismissButton,qt.SIGNAL("clicked()"),self.close)
        self.connect(self.__startButton,qt.SIGNAL("clicked()"),self.start)

    def setFileList(self,filelist=None):
        if filelist is None:filelist = []
        if True or self.__goodFileList(filelist):
            text = ""
            #respect initial file list choice
            #filelist.sort()
            for ffile in filelist:
                text += "%s\n" % ffile
            self.fileList = filelist
            self.__listView.setText(text)
            if len(filelist):
                PyMcaDirs.inputDir = os.path.dirname(filelist[0])
        
        
    def setOutputDir(self,outputdir=None):
        if outputdir is None:return
        if self.__goodOutputDir(outputdir):
            self.outputDir = outputdir
            self.__outLine.setText(outputdir)
            PyMcaDirs.outputDir = self.outputDir
        else:
            qt.QMessageBox.critical(self, "ERROR",
            "Cannot use output directory:\n%s"% (outputdir))

    def __goodFileList(self,filelist):
        if not len(filelist):return True
        for file in filelist:
            if not os.path.exists(file):
                qt.QMessageBox.critical(self, "ERROR",'File %s\ndoes not exists' % file)
                self.raiseW()
                return False
        return True
        
    def __goodOutputDir(self,outputdir):
        if os.path.isdir(outputdir):return True
        else:return False

    def browseList(self):
        if self.inputDir is None:self.inputDir = PyMcaDirs.inputDir
        if not os.path.exists(self.inputDir):
            self.inputDir =  os.getcwd()
        wdir = self.inputDir
        if QTVERSION < '4.0.0':
            filedialog = qt.QFileDialog(self,"Open a set of files",1)
            filedialog.setMode(filedialog.ExistingFiles)
            filedialog.setDir(wdir)
            filedialog.setFilters("Mca Files (*.mca)\nSpec Files (*.dat)\nAll Files (*)\n")
            if filedialog.exec_loop() == qt.QDialog.Accepted:
                filelist0=filedialog.selectedFiles()
            else:
                self.raiseW()
                return
        else:
            filedialog = qt.QFileDialog(self)
            filedialog.setWindowTitle("Open a set of files")
            filedialog.setDirectory(wdir)
            filedialog.setFilters(["Mca Files (*.mca)",
                                   "Spec Files (*.dat)",
                                   "All Files (*)"])
            filedialog.setModal(1)
            filedialog.setFileMode(filedialog.ExistingFiles)
            ret = filedialog.exec_()
            if  ret == qt.QDialog.Accepted:
                filelist0=filedialog.selectedFiles()
            else:
                self.raise_()
                return
            
        filelist = []
        for f in filelist0:
            filelist.append(qt.safe_str(f))
        if len(filelist):
            self.setFileList(filelist)
        if QTVERSION < '4.0.0':
            self.raiseW()
        else:
            self.raise_()

    def browseConfig(self):
        filename = qt.QFileDialog(self,"Open a new fit config file",1)
        filename.setMode(filename.ExistingFiles)
        filename.setFilters("Config Files (*.cfg)\nAll files (*)")
        if filename.exec_loop() == qt.QDialog.Accepted:
            filename = filename.selectedFile()
        else:
            self.raiseW()
            return
        filename= qt.safe_str(filename)
        if len(filename):
            self.setConfigFile(filename)  
        self.raiseW()

    def browseOutputDir(self):
        if self.outputDir is None:
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
                outdir = qt.safe_str(outfile.selectedFile())
            else:
                outdir = qt.safe_str(outfile.selectedFiles()[0])
            outfile.close()
            del outfile
            self.setOutputDir(outdir)
        else:
            # pyflakes http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=666494
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
        if (self.outputDir is None) or (not self.__goodOutputDir(self.outputDir)):
            qt.QMessageBox.critical(self, "ERROR",'Invalid output directory')
            if QTVERSION < '4.0.0':
                self.raiseW()
            else:
                self.raise_()
            return
        name = "Batch from %s to %s " % (os.path.basename(self.fileList[ 0]),
                                          os.path.basename(self.fileList[-1]))

        window =  Mca2EdfWindow(name="Mca 2 Edf "+name,actions=1)
        self.fileStep = int(qt.safe_str(self.__fileSpin.text()))
        b = Mca2EdfBatch(window,self.fileList,self.outputDir, self.fileStep)
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


class Mca2EdfBatch(qt.QThread):
    def __init__(self, parent, filelist=None, outputdir = None, filestep = 1):
        self._filelist  = filelist
        self.outputdir = outputdir
        self.filestep  = filestep
        qt.QThread.__init__(self)
        self.parent = parent
        self.pleasePause = 0

    def processList(self):
        self.__ncols = None
        self.__nrows = self.filestep
        counter = 0
        ffile   = SpecFileLayer.SpecFileLayer()
        for fitfile in self._filelist:
            self.onNewFile(fitfile, self._filelist) 
            ffile.SetSource(fitfile)
            fileinfo = ffile.GetSourceInfo()
            # nscans = len(fileinfo['KeyList'])
            for scankey in  fileinfo['KeyList']:
                scan,order = scankey.split(".")
                info,data  = ffile.LoadSource(scankey)
                scan_obj = ffile.Source.select(scankey)
                if info['NbMca'] > 0:
                    for i in range(info['NbMca']):
                        point = int(i/info['NbMcaDet']) + 1
                        mca   = (i % info['NbMcaDet'])  + 1 
                        key = "%s.%s.%05d.%d" % (scan,order,point,mca)
                        if i == 0:
                            mcainfo,mcadata = ffile.LoadSource(key)
                        mcadata = scan_obj.mca(i+1)
                        y0 = numpy.array(mcadata, numpy.float)
                        if counter == 0:
                            key0 = "%s key %s" % (os.path.basename(fitfile), key)
                            self.__ncols = len(y0)
                            image = numpy.zeros((self.__nrows,self.__ncols), \
                                                numpy.float)
                        if self.__ncols !=  len(y0):
                            print("spectrum has different number of columns")
                            print("skipping it")
                        else:
                            image[counter,:] = y0[:]
                            if (counter+1) == self.filestep:
                                if self.filestep > 1:
                                    key1    = "%s key %s" % (os.path.basename(fitfile), key)
                                    title = "%s to %s" % (key0, key1)
                                else:
                                    title = key0
                                if 1:
                                    ddict={}
                                    if 'Channel0' in mcainfo:
                                        ddict['MCA start ch'] =\
                                                   int(mcainfo['Channel0'])
                                    if 'McaCalib' in mcainfo:
                                        ddict['MCA a'] = mcainfo['McaCalib'][0]
                                        ddict['MCA b'] = mcainfo['McaCalib'][1]
                                        ddict['MCA c'] = mcainfo['McaCalib'][2]
                                else:
                                    ddict = mcainfo
                                ddict['Title'] = title
                                edfname = os.path.join(self.outputdir,title.replace(" ","_")+".edf")
                                edfout  = EdfFile.EdfFile(edfname)
                                edfout.WriteImage (ddict , image, Append=0)
                                counter = 0
                            else:
                                counter += 1

        self.onEnd()

    def run(self):
        self.processList()

    def onNewFile(self, file, filelist):
        qt.QApplication.postEvent(self.parent, McaCustomEvent.McaCustomEvent({'file':file,
                                                                   'filelist':filelist,
                                                                   'event':'onNewFile'}))
        if self.pleasePause:self.__pauseMethod()

    def onEnd(self):
        qt.QApplication.postEvent(self.parent, McaCustomEvent.McaCustomEvent({'event':'onEnd'}))
        if self.pleasePause:self.__pauseMethod()
        

    def __pauseMethod(self):
        qt.QApplication.postEvent(self.parent, McaCustomEvent.McaCustomEvent({'event':'batchPaused'}))
        while(self.pleasePause):
            time.sleep(1)
        qt.QApplication.postEvent(self.parent, McaCustomEvent.McaCustomEvent({'event':'batchResumed'}))
            

class Mca2EdfWindow(qt.QWidget):
    def __init__(self,parent=None, name="BatchWindow", fl=0, actions = 0):
        if qt.qVersion() < '4.0.0':
            qt.QWidget.__init__(self, parent, name, fl)
            self.setCaption(name)
        else:
            qt.QWidget.__init__(self, parent)
            self.setWindowTitle(name)
        self.l = qt.QVBoxLayout(self)
        self.l.setMargin(0)
        self.l.setSpacing(0)
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
        
        self.barsLayout.addWidget(self.progressLabel,0,0)        
        self.barsLayout.addWidget(self.progressBar,0,1)
        self.status      = qt.QLabel(self)
        self.l.addWidget(self.status)
        self.status.setText(" ")
        self.timeLeft      = qt.QLabel(self)
        self.l.addWidget(self.timeLeft)
        self.timeLeft.setText("Estimated time left = ???? min")
        self.time0 = None
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
                                                               event.dict['filelist'])
        elif event.dict['event'] == 'onEnd':    self.onEnd(event.dict)

        elif event.dict['event'] == 'batchPaused': self.onPause()
        
        elif event.dict['event'] == 'batchResumed':self.onResume()

        else:
            print("Unhandled event %s" % event)
                                

    def onNewFile(self, file, filelist):
        indexlist = range(0,len(filelist))
        index  = indexlist.index(filelist.index(file))
        nfiles = len(indexlist)
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

    def onEnd(self,ddict):
        if QTVERSION < '4.0.0':
            n = self.progressBar.progress()
            self.progressBar.setProgress(n+1)
        else:
            n = self.progressBar.value()
            self.progressBar.setValue(n+1)            
        self.status.setText  ("Batch Finished")
        self.timeLeft.setText("Estimated time left = 0 sec")
        if self.actions:
            self.pauseButton.hide()
            self.abortButton.setText("OK")

    def onPause(self):    
        pass

    def onResume(self):    
        pass

def main():
    import getopt
    options     = 'f'
    longoptions = ['outdir=', 'listfile=', 'mcastep=']
    filelist = None
    outdir   = None
    listfile = None
    mcastep  = 1
    opts, args = getopt.getopt(
                    sys.argv[1:],
                    options,
                    longoptions)
    for opt,arg in opts:
        if opt in ('--outdir'):
            outdir = arg
        elif opt in  ('--listfile'):
            listfile  = arg
        elif opt in  ('--mcastep'):
            mcastep  = int(arg)
    if listfile is None: 
        filelist=[]
        for item in args:
            filelist.append(item)
    else:
        fd = open(listfile)
        filelist = fd.readlines()
        fd.close()
        for i in range(len(filelist)):
            filelist[i]=filelist[i].replace('\n','')
    app=qt.QApplication(sys.argv) 
    winpalette = qt.QPalette(qt.QColor(230,240,249),qt.QColor(238,234,238))
    app.setPalette(winpalette)       
    if len(filelist) == 0:
        qt.QObject.connect(app,qt.SIGNAL("lastWindowClosed()"),app, qt.SLOT("quit()"))
        w = Mca2EdfGUI(actions=1)
        if QTVERSION < '4.0.0':
            app.setMainWidget(w)
            w.show()
            app.exec_loop()
        else:
            w.show()
            sys.exit(app.exec_())
    else:
        qt.QObject.connect(app,qt.SIGNAL("lastWindowClosed()"),app, qt.SLOT("quit()"))
        text = "Batch from %s to %s" % (os.path.basename(filelist[0]), os.path.basename(filelist[-1]))
        window =  Mca2EdfWindow(name=text,actions=1)
        b = Mca2EdfBatch(window,filelist,outdir,mcastep)
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
        window.show()
        b.start()
        if QTVERSION < '4.0.0':
            app.setMainWidget(window)
            app.exec_loop()
        else:
            sys.exit(app.exec_())
                
if __name__ == "__main__":
    main()
 
# Mca2Edf.py  --outdir=/tmp --mcastep=1 *.mca
 
# PyMcaBatch.py --cfg=/mntdirect/_bliss/users/sole/COTTE/WithLead.cfg --outdir=/tmp/   /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0007.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0008.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0009.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0010.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0011.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0012.edf /mntdirect/_bliss/users/sole/COTTE/ch09/ch09__mca_0003_0000_0013.edf &
# PyMcaBatch.exe --cfg=E:/COTTE/WithLead.cfg --outdir=C:/tmp/   E:/COTTE/ch09/ch09__mca_0003_0000_0007.edf E:/COTTE/ch09/ch09__mca_0003_0000_0008.edf
