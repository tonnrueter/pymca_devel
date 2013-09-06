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
import os.path

from PyMca import QtBlissGraph
from PyMca import PyMcaQt as qt

if not hasattr(qt, 'QString'):
    QString = qt.safe_str
    QStringList = list
else:
    QString = qt.QString
    QStringList = qt.QStringList
QTVERSION = qt.qVersion()

if QTVERSION > '4.0.0':
    QT4 = True
    try:
        from PyMca import QPyMcaMatplotlibSave
        MATPLOTLIB = True
    except ImportError:
        MATPLOTLIB = False
else:
    QT4 = False
    MATPLOTLIB = False
from PyMca.PyMca_Icons import IconDict
from PyMca import ColormapDialog
from PyMca import PyMcaPrintPreview
from PyMca import ArraySave
from PyMca import PyMcaDirs
from PyMca import SpecFileDataInfo

DEBUG = 0
SOURCE_TYPE = 'EdfFile'
__revision__ = "$Revision: 1.35 $"

class EdfFile_StandardArray(qt.QWidget):
    def __init__(self, parent=None, name="Edf_StandardArray", fl=0, images=None, rows=None, cols=None):
        if images is None:images = 1
        if rows is None:rows = 0
        if cols is None:cols = 0        
        qt.QWidget.__init__(self, parent)
        if qt.qVersion() < '4.0.0':
            layout = qt.QGridLayout(self, 4, 2)
            layout.setColStretch(0,0)
            layout.setColStretch(1,1)
        else:
            layout = qt.QGridLayout(self)
        layout.setMargin(5)

        ilab= qt.QLabel("Image:", self)
        self.plab= qt.QLabel("Plot", self)
        self.ylab= qt.QLabel("Columns :", self)

        layout.addWidget(ilab, 0, 0, qt.Qt.AlignRight)
        layout.addWidget(self.plab, 1, 0, qt.Qt.AlignRight)
        layout.addWidget(self.ylab, 2, 0, qt.Qt.AlignRight|qt.Qt.AlignTop)

        self.iCombo= qt.QComboBox(self)
        self.iCombo.setEditable(0)

        self.plotCombo= qt.QComboBox(self)
        self.plotCombo.setEditable(0)
        if qt.qVersion() < '4.0.0':
            self.plotCombo.insertItem("Rows")
            self.plotCombo.insertItem("Columns")
            self.yList= qt.QListBox(self)
            self.yList.setSelectionMode(qt.QListBox.Multi)
        else:
            self.plotCombo.insertItems(0, ["Rows", "Columns"])
            self.yList= qt.QListWidget(self)
            #self.yList.setSelectionMode(qt.QListBox.Multi)


        layout.addWidget(self.iCombo,   0, 1)
        layout.addWidget(self.plotCombo,1, 1)
        layout.addWidget(self.yList,    2, 1)

        self.connect(self.plotCombo, qt.SIGNAL("activated(int)"), self.__plotChanged)
        self.connect(self.iCombo, qt.SIGNAL("activated(int)"),    self.__iChanged)
        self.setImages(images)

        self.setDataSize(rows, cols)

    def setImages(self,images,info=None):
        self.iCombo.clear()
        if info is None: info = []
        for i in range(images):
            if qt.qVersion() < '4.0.0':
                if len(info) == images:
                    self.iCombo.insertItem("Image %d Key %s" % (i,info[i]))  
                else:
                    self.iCombo.insertItem("Image %d" % i)
            else:
                if len(info) == images:
                    self.iCombo.insertItem(i, "Image %d Key %s" % (i,info[i])) 
                else:
                    self.iCombo.insertItem(i, "Image %d" % i)

    def setCurrentImage(self,image):
        if image < self.iCombo.count():
            if QT4:self.iCombo.setCurrentIndex(image)
            else:  self.iCombo.setCurrentItem(image)
    
    def setDataSize(self, rows, cols):
        self.rows= rows
        self.cols= cols

        idx= self.cols<=self.rows
        if qt.qVersion() < '4.0.0':
            self.plotCombo.setCurrentItem(idx)
        else:
            self.plotCombo.setCurrentIndex(idx)
        self.__plotChanged(idx)

    def __plotChanged(self, index):
        if index==1:        
            self.ylab.setText('Columns')
            txt= "Column"
            val= self.cols
        else:
            self.ylab.setText('Rows')
            txt= "Row"
            val= self.rows
        if qt.qVersion() < '4.0.0': self.yList.clear()
        else:self.yList.clear()
        for x in range(val):
            if QT4:self.yList.addItem("%s %d"%(txt,x))
            else:  self.yList.insertItem("%s %d"%(txt,x))
        dict={}
        dict['event'] = "plotChanged"
        dict['plot']  =  txt+"s"
        if qt.qVersion() < '4.0.0':
            self.emit(qt.PYSIGNAL("widgetSignal"),(dict,))
        else:
            self.emit(qt.SIGNAL("widgetSignal"),(dict))

    def __iChanged(self, index):
        dict={}
        dict['event'] = "imageChanged"
        dict['index'] =  index
        if qt.qVersion() < '4.0.0':
            self.emit(qt.PYSIGNAL("widgetSignal"),(dict,))
        else:
            self.emit(qt.SIGNAL("widgetSignal"),(dict))

    def getSelection(self):
        selection= []

        if QTVERSION < '4.0.0':
            idx = self.plotCombo.currentItem()
        else:
            idx = self.plotCombo.currentIndex()
        if idx==1: plot= "cols"
        else: plot= "rows"

        if qt.qVersion() < '4.0.0':
            idx = self.iCombo.currentItem()
        else:
            idx = self.iCombo.currentIndex()
        if idx==0: image= None
        else: image= idx-1

        if QTVERSION < '4.0.0':
            ylist= [ idx for idx in range(self.yList.count()) if self.yList.isSelected(idx) ]
        else:
            ylist= [ idx for idx in range(self.yList.count()) if self.yList.isItemSelected(self.yList.item(idx)) ]
        for y in ylist:
            selection.append({"plot":plot, "image": image,"x":None, "y":y})
        return selection

    def markImageSelected(self,imagelist=[]):
        if qt.qVersion() < '4.0.0':
            current = self.iCombo.currentItem()
        else:
            current = self.iCombo.currentIndex()
        images  = self.iCombo.count()
        #self.iCombo.clear()
        msg = " (selected)"
        for i in range(images):
            index = "%d" % i
            if qt.qVersion() < '4.0.0':
                text = qt.safe_str(self.iCombo.text(i)).split(msg)[0]
            else:
                text = qt.safe_str(self.iCombo.itemText(i)).split(msg)[0]
            key  = text.split()[-1]
            if qt.qVersion() < '4.0.0':
                if key in imagelist:
                    self.iCombo.changeItem("%s%s" % (text,msg),i)
                else:
                    self.iCombo.changeItem("%s" % (text),i)
            else:
                if key in imagelist:
                    self.iCombo.setItemText(i, "%s%s" % (text,msg))
                else:
                    self.iCombo.setItemText(i, "%s" % (text))
        if qt.qVersion() < '4.0.0':
            self.iCombo.setCurrentItem(current)
        else:
            self.iCombo.setCurrentIndex(current)
        
        
    def markRowSelected(self, rowlist=[]):
        if not qt.safe_str(self.plotCombo.currentText()) == "Rows":
            return
        current = self.yList.currentItem()
        n       = self.yList.count()
        self.yList.clear()
        for index in range(n):
            if qt.qVersion() < '4.0.0':
                if index in rowlist:
                    self.yList.insertItem(" Row %d (selected)" % index)
                else:
                    self.yList.insertItem(" Row %d" % index)
            else:
                if index in rowlist:
                    self.yList.addItem(" Row %d (selected)" % index)
                else:
                    self.yList.addItem(" Row %d" % index)
        #print "asking set"
        #self.yList.setCurrentItem(current)
        #print "DONE"
    
    def markColSelected(self, collist=[]):
        if not qt.safe_str(self.plotCombo.currentText()) == "Columns":
            return
        current = self.yList.currentItem()
        n       = self.yList.count()
        self.yList.clear()
        if QTVERSION < '4.0.0':
            for index in range(n):
                if index in collist:
                    self.yList.insertItem(" Column %d (selected)" % index)
                else:
                    self.yList.insertItem(" Column %d" % index)            
            self.yList.setCurrentItem(current)
        else:
            for index in range(n):
                if index in collist:
                    self.yList.addItem(" Column %d (selected)" % index)
                else:
                    self.yList.addItem(" Column %d" % index)            
            self.yList.setCurrentItem(current)
        

class QEdfFileWidget(qt.QWidget):
    def __init__(self, parent=None, justviewer=False):
        if qt.qVersion() < '4.0.0':
            qt.QWidget.__init__(self, parent)
        else:
            qt.QWidget.__init__(self, parent)
        self.justViewer = justviewer
        self.dataSource= None
        self.oldsource = ""
        self.oldcurrentArray = None
        self.data= None
        self.currentFile= None
        self.currentArray= 0
        self._matplotlibSaveImage = None
        self.selection= None
        self.__plotting = "Columns"
        self._edfstack = None
        self.lastInputDir = None
        self.colormapDialog = None
        self.colormap  = None
        self.printPreview = PyMcaPrintPreview.PyMcaPrintPreview(modal = 0)
        if DEBUG:
            print("printPreview id = %d" % id(self.printPreview))

        #self.selectPixmap= qt.QPixmap(icons.selected)
        #self.unselectPixamp= qt.QPixmap(icons.unselected)
        self.mapComboName= {}

        self.mainLayout= qt.QVBoxLayout(self)
        self.toolBar = None
        self._buildToolBar()

        # --- splitter
        self.splitter= qt.QSplitter(self)
        if QT4:
            self.splitter.setOrientation(qt.Qt.Vertical)        
        else:
            self.splitter.setOrientation(qt.QSplitter.Vertical)
    
        # --- graph
        self.graph=QtBlissGraph.QtBlissGraph(self.splitter)
        self.graph.canvas().setMouseTracking(1)
        self.graph.setTitle('')
        self.graph.xlabel('Columns')
        self.graph.ylabel('Rows')
        if QTVERSION < '4.0.0':
            self.connect(self.graph,qt.PYSIGNAL('QtBlissGraphSignal')  ,
                         self.widgetSignal)
        else:
            self.connect(self.graph,qt.SIGNAL('QtBlissGraphSignal')  ,
                         self.widgetSignal)
        self._x1Limit = self.graph.getx1axislimits()[-1]
        self._y1Limit = self.graph.gety1axislimits()[-1]
        #self.graph.hide()
        # --- array parameter
        self.__dummyW = qt.QWidget(self.splitter)
        self.__dummyW.layout =qt.QVBoxLayout(self.__dummyW)
        self.__dummyW.layout.setMargin(0)
        self.__dummyW.layout.setSpacing(0)
        if not justviewer:
            if not QT4:
                self.applygroup = qt.QHButtonGroup(self.__dummyW,"")
                self.applytoone = qt.QCheckBox(self.applygroup)
                self.applytoone.setText("Apply to seen  image")
                self.applytoone.setChecked(1)
                self.applytoall = qt.QCheckBox(self.applygroup)
                self.applytoall.setText("Apply to all in file")
                self.applygroup.insert(self.applytoone,0)
                self.applygroup.insert(self.applytoall,1)
                self.applygroup.setExclusive(1)
                self.__dummyW.layout.addWidget(self.applygroup)            
                self.applygroup.setFlat(1)
                self.applygroup.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.MinimumExpanding,
                                                             qt.QSizePolicy.Fixed))
                self.connect(self.applygroup,qt.SIGNAL("clicked(int)"),self.groupSignal)
            else:
                self.applygroupContainer = qt.QWidget(self.__dummyW)
                self.applytoone = qt.QCheckBox(self.applygroupContainer)
                self.applytoone.setText("Apply to seen  image")
                self.applytoone.setChecked(1)
                self.applytoall = qt.QCheckBox(self.applygroupContainer)
                self.applytoall.setText("Apply to all images in list")
                self.applygroup = qt.QButtonGroup()
                self.applygroup.addButton(self.applytoone, 0)
                self.applygroup.addButton(self.applytoall, 1)
                self.applygroup.setExclusive(True)
                self.applygroupLayout = qt.QHBoxLayout(self.applygroupContainer)
                self.applygroupLayout.setMargin(0)
                self.applygroupLayout.setSpacing(0)
                self.applygroupLayout.addWidget(self.applytoone)
                self.applygroupLayout.addWidget(self.applytoall)
                self.__dummyW.layout.addWidget(self.applygroupContainer) 
                self.connect(self.applygroup,qt.SIGNAL("buttonClicked(int)"),
                             self.groupSignal)

        self.dataInfoWidgetDict = {}
        self.paramWidget = EdfFile_StandardArray(self.__dummyW)
        self.__dummyW.layout.addWidget(self.paramWidget)
        if QTVERSION < '4.0.0':
            self.connect(self.paramWidget,
                     qt.PYSIGNAL("widgetSignal"),
                     self.widgetSignal)
        else:
            self.connect(self.paramWidget,
                     qt.SIGNAL("widgetSignal"),
                     self.widgetSignal)

        if justviewer:
            self.paramWidget.plab.hide()
            self.paramWidget.plotCombo.hide()
            self.paramWidget.ylab.hide()
            self.paramWidget.yList.hide()
            
        self.allImages = 0
        # --- main layout
        self.mainLayout.setMargin(5)
        self.mainLayout.setSpacing(2)

        #self.mainLayout.addWidget(self.infoBar)
        self.mainLayout.addWidget(self.splitter)
        if not justviewer: self._buildActions()

    def _buildToolBar(self):
        if QTVERSION < '4.0.0':
            if qt.qVersion() < '3.0':
                self.colormapIcon= qt.QIconSet(qt.QPixmap(IconDict["colormap16"]))
            else:
                self.colormapIcon= qt.QIconSet(qt.QPixmap(IconDict["colormap"]))
            self.zoomResetIcon	= qt.QIconSet(qt.QPixmap(IconDict["zoomreset"]))
            self.printIcon	= qt.QIconSet(qt.QPixmap(IconDict["fileprint"]))
            self.saveIcon	= qt.QIconSet(qt.QPixmap(IconDict["filesave"]))
            self.infoIcon = None
        else:
            self.colormapIcon   = qt.QIcon(qt.QPixmap(IconDict["colormap"]))
            self.zoomResetIcon	= qt.QIcon(qt.QPixmap(IconDict["zoomreset"]))
            self.printIcon	= qt.QIcon(qt.QPixmap(IconDict["fileprint"]))
            self.saveIcon	= qt.QIcon(qt.QPixmap(IconDict["filesave"]))
            try:
                self.infoIcon	= qt.QApplication.style().\
                                  standardIcon(qt.QStyle.SP_MessageBoxInformation)
            except:
                self.infoIcon = None

        self.toolBar = qt.QWidget(self)
        self.toolBarLayout = qt.QHBoxLayout(self.toolBar)
        self.toolBarLayout.setMargin(0)
        self.toolBarLayout.setMargin(2)
        self.mainLayout.addWidget(self.toolBar)
        #Autoscale
        self._addToolButton(self.zoomResetIcon,
                            self._zoomReset,
                            'Auto-Scale the Graph')

        #colormap
        self._addToolButton(self.colormapIcon,
                            self.selectColormap,
                            'Color-Scale the Graph')
        
        #info
        if self.infoIcon is not None:
            self._addToolButton(self.infoIcon,
                             self._showInformation,
                            'Show source information')

        #save
        if MATPLOTLIB:
            tb = self._addToolButton(self.saveIcon,
                                 self.__saveIconSignal,
                                 'Export Graph')
            self._saveMenu = qt.QMenu()
            self._saveMenu.addAction(QString("Standard"),    self._saveIconSignal)
            self._saveMenu.addAction(QString("Matplotlib") , self._saveMatplotlibImage)            
        else:
            tb = self._addToolButton(self.saveIcon,
                                 self._saveIconSignal,
                                 'Export Graph')

        #info
        self.infoText = qt.QLabel(self.toolBar)
        self.infoText.setText("    X = ???? Y = ???? Z = ????")
        self.toolBarLayout.addWidget(self.infoText)

        self.toolBarLayout.addWidget(qt.HorizontalSpacer(self.toolBar))

        # ---print
        tb = self._addToolButton(self.printIcon,
                                 self.printGraph,
                                 'Print the Graph')

    def _addToolButton(self, icon, action, tip, toggle=None):
        tb      = qt.QToolButton(self.toolBar)            
        if QTVERSION < '4.0.0':
            tb.setIconSet(icon)
            qt.QToolTip.add(tb,tip) 
            if toggle is not None:
                if toggle:
                    tb.setToggleButton(1)
        else:
            tb.setIcon(icon)
            tb.setToolTip(tip)
            if toggle is not None:
                if toggle:
                    tb.setCheckable(1)
        self.toolBarLayout.addWidget(tb)
        self.connect(tb,qt.SIGNAL('clicked()'), action)
        return tb

    def _showInformation(self):
        if (self.data is None) or \
           (self.currentArray is None):
            qt.QMessageBox.information(self, "No data",\
                                       "No information to be shown")
            return

        #this could be cached because implies a new reading
        infoSource= self.data.getSourceInfo()
        info = self.data.getKeyInfo(infoSource['KeyList']\
                                             [self.currentArray])
        infoWidget = SpecFileDataInfo.SpecFileDataInfo(info, parent=None)
        infoWidget.show()
        infoWidget.notifyCloseEventToWidget(self)
        self.dataInfoWidgetDict[id(infoWidget)] = infoWidget

    def _dataInfoClosed(self, ddict):
        if ddict['event'] == "SpecFileDataInfoClosed":
            key = ddict['id']
            if key in self.dataInfoWidgetDict:
                del self.dataInfoWidgetDict[key]

    def customEvent(self, event):
        if hasattr(event, 'dict'):
            ddict = event.dict
            self._dataInfoClosed(ddict)

    def _zoomReset(self):
        if DEBUG:
            print("_zoomReset")
        self.graph.zoomReset()

    def _saveMatplotlibImage(self):
        if self._matplotlibSaveImage is None:
            if (self.currentArray is None) or \
                (self.data is None):
                self._matplotlibSaveImage = QPyMcaMatplotlibSave.SaveImageSetup(None,
                                                                                None)
            else:
                self._matplotlibSaveImage = QPyMcaMatplotlibSave.SaveImageSetup(None,
                                                                                self.lastData)
        else:
            self._matplotlibSaveImage.setImageData(self.lastData)
        self._matplotlibSaveImage.show()
        self._matplotlibSaveImage.raise_()

    def __saveIconSignal(self):
        self._saveMenu.exec_(self.cursor().pos())        

    def _saveIconSignal(self):
        self.lastInputDir = PyMcaDirs.outputDir

        fileTypeList = ["Data *.dat",
                        "ImageData *.tif",
                        "Image *.png",
                        "Image *.jpg",
                        "ZoomedImage *.png",
                        "ZoomedImage *.jpg",
                        "Widget *.png",
                        "Widget *.jpg"]

        outfile = qt.QFileDialog(self)
        outfile.setModal(1)
        if QTVERSION < '4.0.0':
            outfile.setCaption("Output File Selection")
            filterlist = fileTypeList[0]
            for f in fileTypeList:
                filterlist += "\n%s" % f
            outfile.setFilters(filterlist)
            outfile.setMode(outfile.AnyFile)
            outfile.setDir(self.lastInputDir)
            ret = outfile.exec_loop()
        else:
            outfile.setWindowTitle("Output File Selection")
            strlist = QStringList()
            for f in fileTypeList:
                strlist.append(f)
            outfile.setFilters(strlist)
            outfile.setFileMode(outfile.AnyFile)
            outfile.setAcceptMode(outfile.AcceptSave)
            outfile.setDirectory(self.lastInputDir)
            ret = outfile.exec_()

        if not ret: return
        filterused = qt.safe_str(outfile.selectedFilter()).split()
        filetype = filterused[0]
        extension = filterused[1]
        if QTVERSION < '4.0.0':
            outstr = qt.safe_str(outfile.selectedFile())
        else:
            outstr = qt.safe_str(outfile.selectedFiles()[0])
        try:            
            outputFile = os.path.basename(outstr)
        except:
            outputFile  = outstr
        outputDir  = os.path.dirname(outstr)
        self.lastInputDir = outputDir
        PyMcaDirs.outputDir = outputDir
        
        #always overwrite for the time being
        if len(outputFile) < len(extension[1:]):
            outputFile += extension[1:]
        elif outputFile[-4:] != extension[1:]:
            outputFile += extension[1:]
        outputFile = os.path.join(outputDir, outputFile)
        if os.path.exists(outputFile):
            try:
                os.remove(outputFile)
            except:
                qt.QMessageBox.critical(self, "Save Error", "Cannot overwrite existing file")
                return


        tiff = False
        if filetype.upper() == "IMAGEDATA":
            tiff = True
        if (filetype.upper() == "DATA") or tiff:
            if (self.data is None) or \
               (self.currentArray is None):
                qt.QMessageBox.information(self, "No data",\
                                           "No data to be saved")
                return
            i = 0
            for sname in self.data.sourceName:
                if i == 0:
                    selfdatasourceName = sname
                    i = 1
                else:
                    selfdatasourceName += "|"+sname
            key = self.data.getSourceInfo()['KeyList'][self.currentArray]
            label = selfdatasourceName +"_"+"Key"+"_"+key
            try:
                if tiff:
                    ArraySave.save2DArrayListAsMonochromaticTiff([self.lastData],
                                                                 outputFile,
                                                                 labels = [label],
                                                                 dtype=None)
                else:
                    ArraySave.save2DArrayListAsASCII([self.lastData],
                                                     outputFile,
                                                     labels = [label])
            except:
                qt.QMessageBox.critical(self, "Save Error", "%s" % \
                                        sys.exc_info()[1])
                return
        elif filetype.upper() == "IMAGE":
            self.saveGraphImage(outputFile, original=True)
        elif filetype.upper() == "ZOOMEDIMAGE":
            self.saveGraphImage(outputFile,original=False)
        else:
            self.saveGraphWidget(outputFile)

    def saveGraphImage(self, filename,original=True):
        fformat = filename[-3:].upper()
        if original:
            if QTVERSION < '4.0.0':
                pixmap = qt.QPixmap(self.graph.plotImage.image)
            else:
                pixmap = qt.QPixmap.fromImage(self.graph.plotImage.image)
        else:
            pixmap = qt.QPixmap.grabWidget(self.graph.canvas())
        if pixmap.save(filename, fformat):
            return
        else:
            qt.QMessageBox.critical(self, "Save Error", "%s" % sys.exc_info()[1])
            return

    def saveGraphWidget(self, filename):
        fformat = filename[-3:].upper()
        pixmap = qt.QPixmap.grabWidget(self.graph)
        if pixmap.save(filename, fformat):
            return
        else:
            qt.QMessageBox.critical(self, "Save Error", "%s" % sys.exc_info()[1])
            return

    def setSaveDirectory(self, wdir):
        if os.path.exists(wdir):
            self.lastInputDir = wdir
            return True
        else:
            return False
    
    def printGraph(self):
        pixmap = qt.QPixmap.grabWidget(self.graph.canvas())
        self.printPreview.addPixmap(pixmap)
        if self.printPreview.isHidden():
            self.printPreview.show()
        if QTVERSION < '4.0.0':
            self.printPreview.raiseW()
        else:
            self.printPreview.raise_()

    def _buildActions(self):
        if QTVERSION < '4.0.0':
            return self._buildActionsQt3()
        self.buttonBox = qt.QWidget(self)
        buttonBox = self.buttonBox
        self.buttonBoxLayout = qt.QGridLayout(buttonBox)
        self.buttonBoxLayout.setMargin(2)
        self.buttonBoxLayout.setSpacing(2)
        
        self.add2DButton = qt.QPushButton(buttonBox)
        self.add2DButton.setText("ADD 2D")
        self.remove2DButton = qt.QPushButton(buttonBox)
        self.remove2DButton.setText("REMOVE 2D")
        self.replace2DButton = qt.QPushButton(buttonBox)
        self.replace2DButton.setText("REPLACE 2D")

        self.addButton = qt.QPushButton(buttonBox)
        self.addButton.setText("ADD")
        self.removeButton = qt.QPushButton(buttonBox)
        self.removeButton.setText("REMOVE")
        self.replaceButton = qt.QPushButton(buttonBox)
        self.replaceButton.setText("REPLACE")
        
        self.buttonBoxLayout.addWidget(self.add2DButton, 0, 0)
        self.buttonBoxLayout.addWidget(self.remove2DButton, 0, 1)
        self.buttonBoxLayout.addWidget(self.replace2DButton, 0, 2)

        self.buttonBoxLayout.addWidget(self.addButton, 1, 0)
        self.buttonBoxLayout.addWidget(self.removeButton, 1, 1)
        self.buttonBoxLayout.addWidget(self.replaceButton, 1, 2)
        
        self.mainLayout.addWidget(buttonBox)
        
        self.connect(self.add2DButton, qt.SIGNAL("clicked()"), 
                    self._add2DClicked)
        self.connect(self.remove2DButton, qt.SIGNAL("clicked()"), 
                    self._remove2DClicked)
        self.connect(self.replace2DButton, qt.SIGNAL("clicked()"), 
                    self._replace2DClicked)

        self.connect(self.addButton, qt.SIGNAL("clicked()"), 
                    self._addClicked)
        self.connect(self.removeButton, qt.SIGNAL("clicked()"), 
                    self._removeClicked)
        self.connect(self.replaceButton, qt.SIGNAL("clicked()"), 
                    self._replaceClicked)

    def _buildActionsQt3(self):
        self.buttonBox = qt.QWidget(self)
        buttonBox = self.buttonBox
        self.buttonBoxLayout = qt.QHBoxLayout(buttonBox)
        
        self.addButton = qt.QPushButton(buttonBox)
        self.addButton.setText("ADD")
        self.removeButton = qt.QPushButton(buttonBox)
        self.removeButton.setText("REMOVE")
        self.replaceButton = qt.QPushButton(buttonBox)
        self.replaceButton.setText("REPLACE")
        
        self.buttonBoxLayout.addWidget(self.addButton)
        self.buttonBoxLayout.addWidget(self.removeButton)
        self.buttonBoxLayout.addWidget(self.replaceButton)
        
        self.mainLayout.addWidget(buttonBox)
        
        self.connect(self.addButton, qt.SIGNAL("clicked()"), 
                    self._addClicked)

        self.connect(self.removeButton, qt.SIGNAL("clicked()"), 
                    self._removeClicked)

        self.connect(self.replaceButton, qt.SIGNAL("clicked()"), 
                    self._replaceClicked)


    def groupSignal(self,i):
        self.allImages = i

    def widgetSignal(self,dict={}):
        if 'event' in dict:
            if dict['event']    == 'plotChanged':
                self.__plotting = dict['plot']
                self.__refreshSelection()
            elif dict['event']    == 'MouseAt':
                x = round(dict['y'])
                if x < 0: x = 0
                y = round(dict['x'])
                if y < 0: y = 0
                if (self.data is None) or \
                   (self.currentArray is None):
                    self.infoText.setText("    X = %d Y = %d Z = ????" %\
                                                   (y, x))
                else:
                    limits = self.lastData.shape
                    x = min(int(x), limits[0]-1)
                    y = min(int(y), limits[1]-1)

                    z = self.lastData[x, y]
                    self.infoText.setText("    X = %d Y = %d Z = %.4g" %\
                                                   (y, x, z))
            elif dict['event']    == 'MouseClick':
                if self.justViewer:return
                col = min(int(round(dict['x'])), self._x1Limit - 1)
                row = min(int(round(dict['y'])), self._y1Limit - 1)
                if row < 0: row = 0
                if col < 0: col = 0
                if self.data is None: 
                    self.graph.removeImage()
                    wid = self.__getParamWidget('array')
                    wid.setImages(1)
                    return
                if self.data.sourceName is None:return
                if self.selection is None:
                    self.selection = {}
                nsel = {}
                i = 0
                for sname in self.data.sourceName:
                    if i == 0:
                        selfdatasourceName = sname
                        i = 1
                    else:
                        selfdatasourceName += "|"+sname
                nsel['SourceType'] = self.data.sourceType
                nsel['SourceName'] = selfdatasourceName
                nsel['selection']  = None
                key_list = self.data.getSourceInfo()['KeyList']
                if self.currentArray == len(key_list):
                    key = '0.0'
                else:
                    key = key_list[self.currentArray] 
                if self.allImages:
                    arraynamelist = key_list
                else:
                    arraynamelist = [key]
                for key in arraynamelist:
                    nsel['Key']        = key
                    signalsel = {}
                    signalsel['SourceType'] = self.data.sourceType
                    signalsel['SourceName'] = self.data.sourceName
                    signalsel['selection']  = None
                    signalsel['Key']  = key
                    if self.__plotting == 'Rows':
                        ptype = 'rows'
                        nsel[key] = {'rows':[{'y':row,'x':None}],'cols':[]}
                        signalsel['Key'] += ".r.%d" % row
                    else:
                        nsel[key] = {'rows':[],'cols':[{'y':col,'x':None}]}
                        signalsel['Key'] += ".c.%d" % col
                        ptype = 'cols'
                    name = ""

                    i = int(key.split(".")[0])
                    if i > 0:
                        signalsel['legend'] = os.path.basename(self.data.sourceName[i-1]) +" "+signalsel['Key']
                    else:
                        signalsel['legend'] = "EDF Stack "+ os.path.basename(self.data.sourceName[0])+\
                                              " "+signalsel['Key']
                    if self.selection == {}:
                        self.setSelected([nsel],reset=0)
                        if QTVERSION < '4.0.0':
                            self.emit(qt.PYSIGNAL("addSelection"), ([signalsel],))
                        else:
                            self.emit(qt.SIGNAL("addSelection"), [signalsel])
                    elif not (nsel['SourceName'] in self.selection):
                        self.setSelected([nsel],reset=0)
                        if QTVERSION < '4.0.0':
                            self.emit(qt.PYSIGNAL("addSelection"), ([signalsel],))
                        else:
                            self.emit(qt.SIGNAL("addSelection"), [signalsel])
                    elif not (key in self.selection[nsel['SourceName']]):
                        self.setSelected([nsel],reset=0)
                        if QTVERSION < '4.0.0':
                            self.emit(qt.PYSIGNAL("addSelection"), ([signalsel],))
                        else:
                            self.emit(qt.SIGNAL("addSelection"), [signalsel])
                    elif len(self.selection[nsel['SourceName']][key][ptype]) == 0:
                        self.setSelected([nsel],reset=0)
                        if QTVERSION < '4.0.0':
                            self.emit(qt.PYSIGNAL("addSelection"), ([signalsel],))
                        else:
                            self.emit(qt.SIGNAL("addSelection"), [signalsel])
                    elif nsel[key][ptype][0] not in self.selection[nsel['SourceName']][key][ptype]:
                        self.setSelected([nsel],reset=0)
                        if QTVERSION < '4.0.0':
                            self.emit(qt.PYSIGNAL("addSelection"), ([signalsel],))
                        else:
                            self.emit(qt.SIGNAL("addSelection"), [signalsel])
                    else:
                        self.removeSelection([nsel])
            elif dict['event']  == 'imageChanged':
                if DEBUG:
                    print("Image changed")
                if dict['index'] != self.currentArray: 
                    self.currentArray = dict['index']
                    self.refresh()
                if DEBUG:
                    print("self.currentArray = ",self.currentArray)


    def openFileOLD(self, filename=None):
        if DEBUG:
            print("openfile = ",filename)
        if filename is None:
            filename= qt.QFileDialog(self,"Open a new EdfFile", 1)
            filename.setFilters("EdfFiles (*.edf)\nEdfFiles (*.mca)\nEdfFiles (*ccd)\nAll files (*)")
            if filename.exec_loop() == qt.QDialog.Accepted:
                filename= qt.safe_str(filename.selectedFile())
            else:
                return
            if not len(filename):    return
    
        if filename in self.mapComboName.keys():
            self.selectFile(filename)
        else:
            if not self.data.SetSource(filename):
                qt.QMessageBox.critical(self, "ERROR opening EdfFile",
                        "Cannot open following EdfFile:\n%s"%(filename))
            else:
                filename= self.data.SourceName
                self.mapComboName[filename]= os.path.basename(filename)
                self.fileCombo.insertItem(self.mapComboName[filename])
                self.selectFile(filename)

    def openFile(self, filename=None,justloaded=None):
        if DEBUG:
            print("openfile = %s" % filename)
        if justloaded is None:justloaded = 0
        if filename is None:
            self.lastInputDir = PyMcaDirs.inputDir
            if QT4:
                fdialog = qt.QFileDialog(self)
                fdialog.setModal(True)
                fdialog.setWindowTitle("Open a new EdfFile")
                strlist = QStringList()
                strlist.append("EDF Files *edf")
                strlist.append("EDF Files *ccd")
                strlist.append("All Files *")
                fdialog.setFilters(strlist)
                fdialog.setFileMode(fdialog.ExistingFiles)
                ret = fdialog.exec_()
                if ret == qt.QDialog.Accepted:
                    filelist = fdialog.selectedFiles()
                    fdialog.close()
                    del fdialog                        
                else:
                    fdialog.close()
                    del fdialog
                    return            
            elif sys.platform == 'win32':
                wdir = self.lastInputDir
                if wdir is None:wdir = ""
                filelist = qt.QFileDialog.getOpenFileNames("EdfFiles (*.edf)\nEdfFiles (*mca)\nEdfFiles (*ccd)\nAll files (*)",
                            wdir,
                            self,"openFile", "Open a new EdfFile")
            else:
                filedialog = qt.QFileDialog(self,"Open new EdfFile(s)",1)
                if self.lastInputDir is not None:filedialog.setDir(self.lastInputDir)
                filedialog.setMode(filedialog.ExistingFiles)
                filedialog.setFilters("EdfFiles (*.edf)\nEdfFiles (*.mca)\nEdfFiles (*ccd)\nAll files (*)")           
                if filedialog.exec_loop() == qt.QDialog.Accepted:
                    filelist= filedialog.selectedFiles()
                else:
                    return
            #respect selection choice
            #filelist.sort()
            filename=[]
            for f in filelist:
                filename.append(qt.safe_str(f))
            if not len(filename):    return
            if len(filename):
                self.lastInputDir  = os.path.dirname(filename[0])
                PyMcaDirs.inputDir = os.path.dirname(filename[0])
            justloaded = 1
        if justloaded:
            if type(filename) != type([]):
                filename = [filename]
        if not os.path.exists(filename[0]):
            raise IOError("File %s does not exist" % filename[0])
        if (justloaded) and (filename in self.mapComboName.keys()):
            self.selectFile(filename,justloaded=justloaded)
        elif 1:
            combokey = os.path.basename(filename[0])
            self.mapComboName[combokey]= filename
            self.selectFile(combokey,justloaded=justloaded)
        else:
            if not self.data.SetSource(filename):
                qt.QMessageBox.critical(self, "ERROR opening EdfFile",
                        "Cannot open following EdfFile:\n%s"%(filename))
            else:
                filename= self.data.SourceName.split("|")
                if len(filename) > 1:
                    combokey = 'EDF Stack'
                    self._edfstack = filename
                else:
                    combokey = os.path.basename(filename[0])
                if combokey not in self.mapComboName.keys():
                    self.mapComboName[combokey]= filename[0]
                    if QT4:
                        self.fileCombo.addItem(combokey)                    
                    else:
                        self.fileCombo.insertItem(combokey)
                self.selectFile(combokey,justloaded=justloaded)


    def selectFile(self, filename=None, justloaded=None):
        if justloaded is None:justloaded=0
        if filename is not None:
            #if qt.safe_str(self.fileCombo.currentText()) !=\
            #   self.mapComboName[filename]:
            if filename == 'EDF Stack':
                filename= self._edfstack
            else:
                filename = self.mapComboName[filename]
            if justloaded and (filename==self._edfstack):
                self.currentArray=len(self.data.getSourceInfo()['KeyList'])
            else:
                self.currentArray=0
        self.refresh()
    
    def selectColormap(self):
        if self.colormap is None: return
        if self.colormapDialog.isHidden():
            self.colormapDialog.show()
        if qt.qVersion() < '4.0.0':self.colormapDialog.raiseW()
        else:  self.colormapDialog.raise_()          
        self.colormapDialog.show()

    def updateColormap(self, *var):
        if len(var) > 6:
            self.colormap = [var[0],
                             var[1],
                             var[2],
                             var[3],
                             var[4],
                             var[5],
                             var[6]]
        elif len(var) > 5:
            self.colormap = [var[0],
                             var[1],
                             var[2],
                             var[3],
                             var[4],
                             var[5]]
        else:
            self.colormap = [var[0],
                             var[1],
                             var[2],
                             var[3],
                             var[4],
                             var[5]]
        self.graph.setY1AxisInverted(True)
        self.graph.imagePlot(self.lastData,
                         colormap = self.colormap,
                         xmirror = False,
                         ymirror = False)
        self.graph.replot()

    def closeFile(self, filename=None):
        if filename is None:
            ffile= qt.safe_str(self.fileCombo.currentText())
            #if file != "EDF Stack":
            #    filename = self.mapComboName[file]
            #else:
            #    filename ="EDF Stack"
            filename = ffile

        #print self.selection
        if (self.selection is not None) and filename in self.selection:
            nmca = 0
            for key in self.selection[filename].keys():
                nmca += len(self.selection[filename][key]['rows']) + len(self.selection[filename][key]['cols'])
            if nmca:
                msg= "%d mca are linked to that EdfFile source.\n"% nmca
                msg+="Do you really want to delete all these graphs ??"
                ans= qt.QMessageBox.information(self, "Remove SpecFile %s"%filename, msg,
                        qt.QMessageBox.No, qt.QMessageBox.Yes)
                if ans==qt.QMessageBox.No: return
                try:
                    self.emit(qt.PYSIGNAL("delSelection"), (self.data.SourceName, mcakeys))
                except:
                    print("This is to be implemented")
        
        for idx in range(self.fileCombo.count()):
            if qt.qVersion() < '4.0.0':
                itext = self.fileCombo.text(idx)
            else:
                itext = self.fileCombo.itemText(idx)            
            if filename == "EDF Stack":
                if itext == filename:
                    self.fileCombo.removeItem(idx)
                    del self.mapComboName[filename]
                    break
            elif qt.safe_str(itext) ==\
                 os.path.basename(self.mapComboName[filename]):
                self.fileCombo.removeItem(idx)
                del self.mapComboName[filename]
                break

        if not self.fileCombo.count():
            self.data.sourceName = None
            self._reset()
            #self.selectFile()
        else:
            self.selectFile(self.mapComboName.keys()[0])

    def __fileSelection(self, ffile):
        ffile= qt.safe_str(ffile)
        for filename, comboname in self.mapComboName.items():
            if filename == ffile:
                self.selectFile(filename)
                break

    def _reset(self):
        self.graph.removeImage()
        self.oldsource = None
        self.graph.clearmarkers()
        self.graph.replot()
        wid = self.__getParamWidget('array')
        wid.setImages(1)
        wid.setDataSize(0,0)


    def setDataSource(self,data=None):
        if DEBUG:
            print("setData(self, data) called")
            print("data = ",data)
        self.data= data
        self.refresh()

    def refresh(self):
        if DEBUG:
            print("refresh method called")
        if self.data is None:
            self._reset()
            #wid = self.__getParamWidget('array')
            #wid.setImages(1)
            return
        if self.data.sourceName is None:    return
        self.currentFile = self.data.sourceName
        #this gives the number of images in the file
        infoSource= self.data.getSourceInfo()
        if DEBUG:
            print("info :")
            print(infoSource)
        
        nimages=len(infoSource['KeyList'])
        #print self.data.SourceName,"nimages = ",nimages
        loadsum = 0
        if nimages == 1:
            self.currentArray =  0
        elif self.currentArray > nimages:
            self.currentArray =  0
        elif self.currentArray == nimages:
            loadsum=1            
        #print "SUM = ",loadsum, infoSource['KeyList']
        #print self.currentArray
        if (self.oldsource != self.currentFile) or (self.oldcurrentArray != self.currentArray):
            if DEBUG:
                print("I have to read again ... ")
            if not loadsum:
                if DEBUG:
                    print("Not Loading the sum")
                dataObject = self.data.getDataObject(infoSource['KeyList']\
                                                     [self.currentArray])
                info = dataObject.info
                data = dataObject.data
                imageinfo = infoSource['KeyList']
            else:
                if DEBUG:
                    print("Loading the sum")
                dataObject = self.data.getDataObject('0.0')
                info = dataObject.info
                data = dataObject.data           
                imageinfo = infoSource['KeyList']
            wid= self.__getParamWidget("array")
            if nimages > 1:
                if 'Title' in info:
                    i = 0
                    for key in self.data.getSourceInfo()['KeyList']:
                        source,image = key.split(".")
                        source = int(source)
                        image  = int(image)
                        dataObject = self.data.getDataObject(key)
                        header = dataObject.info
                        if 'Title' in header:
                            imageinfo[i] += "- " + header['Title']
                        i+=1
                if DEBUG:
                    print("NOT ADDING 0.0 - SUM KEY")
                    wid.setImages(nimages+1,info = imageinfo+["0.0 - SUM"])
                wid.setImages(nimages,info = imageinfo)
            else:
                if 'Title' in info:
                    imageinfo [self.currentArray] += info['Title']  
                wid.setImages(nimages,  info = imageinfo)                
            wid.setCurrentImage(self.currentArray)
            #P.B. -> pointer(a,d1,d2,i1,i2) = a+ (i1+i2 * d1) 
            wid.setDataSize(int(info["Dim_2"]), int(info["Dim_1"]))
            if DEBUG:
                print("Image size = %d x %d" % (int(info["Dim_2"]),
                                                int(info["Dim_1"])))
                print("data  size = ", data.shape)

            if self.graph.isHidden():
                self.graph.show()
            ##self.graph.setx1axislimits(0, int(info["Dim_2"]))
            ##self.graph.sety1axislimits(0, int(info["Dim_1"]))
            self._x1Limit = int(info["Dim_1"])
            self._y1Limit = int(info["Dim_2"])
            self.graph.clear()
            minData = data.min()
            maxData = data.max()
            wasnone = 0
            self.lastData = data
            if self.colormapDialog is None:
                wasnone = 1
                self.colormapDialog = ColormapDialog.ColormapDialog()
                self.colormapDialog.colormapIndex  = self.colormapDialog.colormapList.index("Temperature")
                self.colormapDialog.colormapString = "Temperature"
                if QTVERSION < '4.0.0':
                    self.connect(self.colormapDialog,
                                 qt.PYSIGNAL("ColormapChanged"),
                                 self.updateColormap)
                else:
                    self.connect(self.colormapDialog,
                                 qt.SIGNAL("ColormapChanged"),
                                 self.updateColormap)
            self.colormapDialog.setDataMinMax(minData, maxData)
            if wasnone:
                self.colormapDialog.setAutoscale(1)
                self.colormapDialog.setColormap(self.colormapDialog.colormapIndex)
            self.colormap = (self.colormapDialog.colormapIndex,
                             self.colormapDialog.autoscale,
                             self.colormapDialog.minValue, 
                             self.colormapDialog.maxValue,
                             minData, maxData)
            #self.graph.imagePlot(data=data, colormap = self.colormap)
            self.colormapDialog._update()
            self.graph.setY1AxisInverted(True)
            self.graph.imagePlot(data=data,
                                 colormap = self.colormap,
                                 ymirror = False)
            
        self.__refreshSelection()
        self.graph.replot()
        self.oldsource       = "%s" % self.data.sourceName
        self.oldcurrentArray = self.currentArray * 1
        
    def __getParamWidget(self, widtype):
        return self.paramWidget

    def _replaceClicked(self):
        if DEBUG:
            print("replace clicked")
        selkeys= self.__getSelectedKeys()
        if len(selkeys):
            #self.eh.event(self.repEvent, selkeys)
            if DEBUG:
                print("Replace event")
            if self.allImages:
                arraynamelist = self.data.getSourceInfo()['KeyList']
            else:
                arraynamelist = []
                for selection in selkeys:
                    arraynamelist.append(selection['Key'])
            sellist=[]
            signalsellist = []
            for arrayname in arraynamelist:
                sel = {}
                sel['SourceType'] = SOURCE_TYPE            
                for selection in selkeys:
                    signalsel = {}
                    signalsel.update(sel)
                    signalsel['selection']  = None
                    signalsel['SourceName'] = self.data.sourceName
                    if not ('SourceName' in sel):
                        sel['SourceName'] = selection['SourceName']
                    arrayname = selection['Key']
                    if not ('Key' in sel):
                        sel['Key'] = selection['Key']
                    signalsel['Key'] = selection['Key']
                    if not (arrayname in sel):
                        sel[arrayname] = {'rows':[],'cols':[]}
                    if selection['plot'] == 'cols':
                        sel[arrayname]['cols'].append({'x':selection['x'],'y':selection['y']})
                        signalsel['Key'] += ".c.%d" % selection['y']
                    if selection['plot'] == 'rows':
                        sel[arrayname]['rows'].append({'x':selection['x'],'y':selection['y']})                              
                        signalsel['Key'] += ".r.%d" % selection['y']
                    i = int(signalsel['Key'].split(".")[0])
                    if i > 0:
                        signalsel['legend'] = os.path.basename(self.data.sourceName[i-1]) +" "+signalsel['Key']
                    else:
                        signalsel['legend'] = "EDF Stack "+ \
                                              os.path.basename(self.data.sourceName[0])+\
                                              " "+signalsel['Key']

                    """
                    if selection['plot'] == 0:
                         sel[arrayname]['mca'].append({'x':selection['x'],'y':selection['y']})
                    """
                    signalsellist.append(signalsel)
                sellist.append(sel)
            self.setSelected(sellist,reset=1)
            if QTVERSION < '4.0.0':
                self.emit(qt.PYSIGNAL("replaceSelection"), (signalsellist,))
            else:
                self.emit(qt.SIGNAL("replaceSelection"), signalsellist)

    def _add2DClicked(self, replace=False):
        if DEBUG:
            print("ADD 2D clicked")
        if (self.data is None) or \
           (self.currentArray is None):
            return

        #this is not very efficient because it could be cached
        #while this implies a new reading
        infoSource= self.data.getSourceInfo()
        sel = {}
        sel['SourceType'] = infoSource['SourceType'] 
        sel['SourceName'] = self.data.sourceName
        sel['Key'] = infoSource['KeyList'][self.currentArray]
        f, i = sel['Key'].split(".")
        f = int(f) - 1
        sel['legend'] = os.path.basename(self.data.sourceName[f]) +\
                        " "+ ("%s" % self.paramWidget.iCombo.currentText())
        sel['selectiontype'] = '2D' 
        sel['imageselection']  = True
        sel['mcaselection']  = False
        sel['scanselection'] = False
        sel['selection'] = None
        if replace:
            self.emit(qt.SIGNAL("replaceSelection"), [sel])
        else:
            self.emit(qt.SIGNAL("addSelection"), [sel])

    def _remove2DClicked(self):
        if DEBUG:
            print("REMOVE 2D clicked")
        infoSource= self.data.getSourceInfo()
        sel = {}
        sel['SourceType'] = infoSource['SourceType'] 
        sel['SourceName'] = self.data.sourceName
        sel['Key'] = infoSource['KeyList'][self.currentArray]
        f, i = sel['Key'].split(".")
        f = int(f) - 1
        sel['legend'] = os.path.basename(self.data.sourceName[f]) +\
                        " "+ qt.safe_str(self.paramWidget.iCombo.currentText())
        sel['selectiontype'] = '2D' 
        sel['imageselection']  = True
        sel['mcaselection']  = False
        sel['scanselection'] = False
        sel['selection'] = None
        self.emit(qt.SIGNAL("removeSelection"), [sel])

    def _replace2DClicked(self):
        if DEBUG:
            print("REPLACE 2D clicked")
        self._add2DClicked(replace=True)

    def _addClicked(self):
        if DEBUG:
            print("select clicked")
        selkeys= self.__getSelectedKeys()
        if DEBUG:
            print("selected keys = ",selkeys) 
        if len(selkeys):
            #self.eh.event(self.addEvent, selkeys)
            if DEBUG:
                print("Select event")
            if self.allImages:
                arraynamelist = self.data.getSourceInfo()['KeyList']
            else:
                arraynamelist = []
                for selection in selkeys:
                    arraynamelist.append(selection['Key'])
            sellist=[]
            sellistsignal = []
            for arrayname in arraynamelist:
                sel = {}
                sel['SourceType'] = SOURCE_TYPE
                for selection in selkeys:
                    selsignal = {}
                    selsignal['SourceType'] = self.data.sourceType
                    selsignal['SourceName'] = self.data.sourceName
                    selsignal['selection'] = None
                    selsignal['Key'] = arrayname
                    if not ('SourceName' in sel):
                        sel['SourceName'] = selection['SourceName']
                    #arrayname = selection['Key']
                    if not ('Key' in sel):
                        sel['Key'] = arrayname
                    if not (arrayname in sel):
                        sel[arrayname] = {'rows':[],'cols':[]}
                    if selection['plot'] == 'cols':
                        sel[arrayname]['cols'].append({'x':selection['x'],
                                                        'y':selection['y']})
                        selsignal["Key"] += ".c.%d" % int(selection['y'])
                        i = int(selsignal["Key"].split(".")[0])
                        if i > 0:
                            selsignal['legend'] = os.path.basename(self.data.sourceName[i-1]) +" "+selsignal['Key']
                        else:
                            selsignal['legend'] = "EDF Stack "+ os.path.basename(self.data.sourceName[0])+\
                                                  " "+selsignal['Key']
                    if selection['plot'] == 'rows':
                        sel[arrayname]['rows'].append({'x':selection['x'],
                                                        'y':selection['y']})
                        selsignal["Key"] += ".r.%d" % int(selection['y'])
                        i = int(selsignal["Key"].split(".")[0])
                        if i > 0:
                            selsignal['legend'] = os.path.basename(self.data.sourceName[i-1]) +\
                                                  " "+selsignal['Key']
                        else:
                            selsignal['legend'] = "EDF Stack "+\
                                                  os.path.basename(self.data.sourceName[0])+\
                                                  " "+selsignal['Key']
                    sellistsignal.append(selsignal)
                sellist.append(sel)
            if self.selection is None: 
                self.setSelected(sellist,reset=1)
            else:
                self.setSelected(sellist,reset=0)
            if QTVERSION < '4.0.0':
                self.emit(qt.PYSIGNAL("addSelection"), (sellistsignal,))
            else:
                self.emit(qt.SIGNAL("addSelection"), sellistsignal)
            
    def __getSelectedKeys(self):
        selkeys= []
        parwid= self.paramWidget
        #.visibleWidget()
        if self.currentArray is not None:
            for sel in parwid.getSelection():
                sel["SourceName"]= self.currentFile
                sel['SourceType'] = SOURCE_TYPE
                if 0:
                    sel["Key"]= "%d" % self.currentArray
                else:
                    keylist = self.data.getSourceInfo()['KeyList']
                    if self.currentArray == len(keylist):
                        sel["Key"]= "0.0"
                    else:
                        sel["Key"]= keylist[self.currentArray]                
                selkeys.append(sel)
        return selkeys

    def _removeClicked(self):
        if DEBUG:
            print("remove clicked")
        selkeys= self.__getSelectedKeys()
        returnedselection=[]
        signalsellist = []
        if len(selkeys):
            #self.eh.event(self.delEvent, selkeys)
            if DEBUG:
                print("Remove Event")
                print("self.selection before = ",self.selection)
            if self.allImages:
                arraynamelist = self.data.getSourceInfo()['KeyList']
            else:
                arraynamelist = []
                for selection in selkeys:
                    arraynamelist.append(selection['Key'])
            for arrayname in arraynamelist:
                for selection in selkeys:
                    sel = {}
                    i = 0
                    for sname in self.data.sourceName:
                        if i == 0:
                            selfdatasourceName = sname
                            i = 1
                        else:
                            selfdatasourceName += "|"+sname
                    sel['SourceName'] = selfdatasourceName
                    sel['SourceType'] = SOURCE_TYPE            
                    #sel['Key'] = selection['Key']
                    #arrayname = "%s" % selection['Key']
                    sel['Key'] = arrayname
                    sel[arrayname] = {'rows':[],'cols':[]}
                    if selection['plot'] == 'cols':
                         sel[arrayname]['cols'].append({'x':selection['x'],'y':selection['y']})
                    if selection['plot'] == 'rows':
                         sel[arrayname]['rows'].append({'x':selection['x'],'y':selection['y']})
                    if self.selection is not None:
                        if DEBUG:
                            print("step 1")
                        if sel['SourceName'] in self.selection:
                            if DEBUG:
                                print("step 2")
                            if arrayname in self.selection[sel['SourceName']]:
                                if DEBUG:
                                    print("step 3")
                                if 'rows' in self.selection[sel['SourceName']][arrayname]:
                                    if DEBUG:
                                        print("step 4")
                                    for couple in  sel[arrayname]['rows']:
                                        if couple in  self.selection[sel['SourceName']][arrayname]['rows']:
                                            index= self.selection[sel['SourceName']][arrayname]['rows'].index(couple)
                                            del self.selection[sel['SourceName']][arrayname]['rows'][index]
                                            signalsel = {}
                                            signalsel.update(sel)
                                            signalsel['SourceName'] = self.data.sourceName
                                            signalsel['Key'] += ".r.%d" % couple['y']
                                            i = int(signalsel['Key'].split(".")[0])
                                            if i > 0:
                                                signalsel['legend'] = os.path.basename(self.data.sourceName[i-1]) +" "+signalsel['Key']
                                            else:
                                                signalsel['legend'] = "EDF Stack "+ \
                                                                      os.path.basename(self.data.sourceName[0])+\
                                                                      " "+signalsel['Key']
                                            signalsellist.append(signalsel)
                                    for couple in  sel[arrayname]['cols']:
                                        if couple in  self.selection[sel['SourceName']][arrayname]['cols']:
                                            index= self.selection[sel['SourceName']][arrayname]['cols'].index(couple)
                                            del self.selection[sel['SourceName']][arrayname]['cols'][index]
                                            signalsel = {}
                                            signalsel.update(sel)
                                            signalsel['SourceName'] = self.data.sourceName
                                            signalsel['Key'] += ".c.%d" % couple['y']
                                            i = int(signalsel['Key'].split(".")[0])
                                            if i > 0:
                                                signalsel['legend'] = os.path.basename(self.data.sourceName[i-1]) +" "+signalsel['Key']
                                            else:
                                                signalsel['legend'] = "EDF Stack "+ \
                                                                      os.path.basename(self.data.sourceName[0])+\
                                                                      " "+signalsel['Key']
                                            signalsellist.append(signalsel)
                                    seln = {}
                                    seln['SourceName'] = sel['SourceName'] 
                                    seln['SourceType'] = SOURCE_TYPE            
                                    seln['Key']        = sel['Key']
                                    seln[seln['Key']]  = self.selection[seln['SourceName']][seln['Key']]
                                    self.setSelected([seln],reset=0)
                    returnedselection.append(sel)
            if QTVERSION < '4.0.0':
                self.emit(qt.PYSIGNAL("removeSelection"), (signalsellist,))
            else:
                self.emit(qt.SIGNAL("removeSelection"), signalsellist)
            
    def removeSelection(self,selection):
        if type(selection) != type([]):
            selection=[selection]
        signalsellist = []
        for sel in selection:
                arrayname = sel['Key']
                if self.selection is not None:
                    if DEBUG:
                        print("step 1")
                    if sel['SourceName'] in self.selection:
                        if DEBUG:
                            print("step 2")
                        if arrayname in self.selection[sel['SourceName']]:
                            if DEBUG:
                                print("step 3")
                            if 'rows' in self.selection[sel['SourceName']][arrayname]:
                                if DEBUG:
                                    print("step 4")
                                for couple in  sel[arrayname]['rows']:
                                    if couple in  self.selection[sel['SourceName']][arrayname]['rows']:
                                        index= self.selection[sel['SourceName']][arrayname]['rows'].index(couple)
                                        del self.selection[sel['SourceName']][arrayname]['rows'][index]
                                        signalsel = {}
                                        signalsel.update(sel)
                                        signalsel['SourceName'] = self.data.sourceName
                                        signalsel['Key'] += ".r.%d" % couple['y']
                                        i = int(signalsel['Key'].split(".")[0])
                                        if i > 0:
                                            signalsel['legend'] = os.path.basename(self.data.sourceName[i-1]) +" "+signalsel['Key']
                                        else:
                                            signalsel['legend'] = "EDF Stack "+ \
                                                                  os.path.basename(self.data.sourceName[0])+\
                                                                  " "+signalsel['Key']
                                        signalsellist.append(signalsel)
                                for couple in  sel[arrayname]['cols']:
                                    if couple in  self.selection[sel['SourceName']][arrayname]['cols']:
                                        index= self.selection[sel['SourceName']][arrayname]['cols'].index(couple)
                                        del self.selection[sel['SourceName']][arrayname]['cols'][index]
                                        signalsel = {}
                                        signalsel.update(sel)
                                        signalsel['SourceName'] = self.data.sourceName
                                        signalsel['Key'] += ".r.%d" % couple['y']
                                        i = int(signalsel['Key'].split(".")[0])
                                        if i > 0:
                                            signalsel['legend'] = os.path.basename(self.data.sourceName[i-1]) +" "+signalsel['Key']
                                        else:
                                            signalsel['legend'] = "EDF Stack "+ \
                                                                  os.path.basename(self.data.sourceName[0])+\
                                                                  " "+signalsel['Key']
                                        signalsellist.append(signalsel)
                                seln = {}
                                seln['SourceName'] = sel['SourceName'] 
                                seln['SourceType'] = SOURCE_TYPE            
                                seln['Key']        = sel['Key']
                                seln[seln['Key']]  = self.selection[seln['SourceName']][seln['Key']]
                                self.setSelected([seln],reset=0)
        if QTVERSION < '4.0.0':
            self.emit(qt.PYSIGNAL("removeSelection"), (signalsellist,))
        else:
            self.emit(qt.SIGNAL("removeSelection"), signalsellist)


                             
    def setSelected(self,sellist,reset=1):
        if DEBUG:
            print("setSelected(self,sellist,reset=1) called")
            print("sellist = ",sellist)
            print("selection before = ",self.selection)
            print("reset = ",reset)
        if reset:
            self.selection = {}
        elif self.selection is None:
            self.selection = {}
        for sel in sellist:
            specname = sel['SourceName']
            if type(specname) == type([]):
                for i in range(len(sel['SourceName'])):
                    if i == 0:
                        specname = sel['SourceName'][i]
                    else:
                        specname += "|"+sel['SourceName'][i]
            #selkey is the array name what to do if multiple array names?
            if type(sel["Key"]) == type([]):
                selkey = sel["Key"][0]
            else:
                selkey = sel["Key"]
            if not (specname in self.selection):
                self.selection[specname]= {}
            if not (selkey in self.selection[specname]):
                self.selection[specname][selkey] = {'rows':[],'cols':[]}
            if 'rows' in sel[selkey]:
                for rowsel in sel[selkey]['rows']:
                    if rowsel not in self.selection[specname][selkey]['rows']:
                        self.selection[specname][selkey]['rows'].append(rowsel)   
            if 'cols' in sel[selkey]:
                for rowsel in sel[selkey]['cols']:
                    if rowsel not in self.selection[specname][selkey]['cols']:
                        self.selection[specname][selkey]['cols'].append(rowsel)   
        if DEBUG:
            print("self.selection after = ",self.selection)
        self.__refreshSelection()

    def getSelection(self):
        """
        Give the dicionary of dictionaries as an easy to understand list of
        individual selections
        """
        selection = []
        if self.selection is None: return selection
        for sourcekey in self.selection.keys():
            for arraykey in self.selection[sourcekey].keys():
                sel={}
                sel['SourceName']   = sourcekey
                sel['SourceType']   = 'EdfFile'
                sel['Key']          = arraykey
                sel[arraykey]        = self.selection[sourcekey][arraykey]
                selection.append(sel)
        return selection

        
    def __refreshSelection(self):
        if DEBUG:
            print("__refreshSelection(self) called")
            print(self.selection)
            print("self.data.SourceName = ",self.data.sourceName)
        if self.selection is not None:
            if self.data is None:return
            if self.data.sourceName is None: return
            if type(self.data.sourceName) == type([]):
                i = 0
                for sname in self.data.sourceName:
                    if i == 0:
                        selfdatasourceName = sname
                        i = 1
                    else:
                        selfdatasourceName += "|"+sname
            else:
                selfdatasourceName = self.data.sourceName
            if "|" in self.data.sourceName:
                #print "here should be the multiple" 
                #sel = self.selection.get(self.data.SourceName[0], {})
                sel = self.selection.get(selfdatasourceName, {})
            else:
                sel = self.selection.get(selfdatasourceName, {})
            selkeys = []
            for key in sel.keys():
                if (sel[key]['rows'] != []) or (sel[key]['cols'] !=  []):
                    selkeys.append(key)
            if DEBUG:
                print("selected images =",selkeys,"but self.selection = ",self.selection)
                print("and self.selection.get(self.data.SourceName, {}) =",sel)
            
            wid = self.__getParamWidget("array")
            wid.markImageSelected(selkeys)
            #imagedict = sel.get("%d" % self.currentArray, {})
            keylist = self.data.getSourceInfo()['KeyList']
            if self.currentArray == len(keylist):
                imagedict = sel.get("0.0",{})
            else:
                imagedict = sel.get(keylist[self.currentArray],{})
            if not ('rows' in imagedict):
                imagedict['rows'] = []
            if not ('cols' in imagedict):
                imagedict['cols'] = []
            rows = []
            for dict in imagedict['rows']:
                if 'y' in dict:
                    if dict['y'] not in rows:
                        rows.append(dict['y'])
            wid.markRowSelected(rows) 
            cols = []
            for dict in imagedict['cols']:
                if 'y' in dict:
                    if dict['y'] not in cols:
                        cols.append(dict['y'])            
            wid.markColSelected(cols)
            self.graph.clearmarkers()
            for i in rows:
                label = "R%d" % i
                marker=self.graph.inserty1marker(0.1,i,label=label)
                self.graph.setmarkercolor(marker,"white")
            for i in cols:
                label = "C%d" % i
                marker=self.graph.insertx1marker(i, 0.1,label=label)
                self.graph.setmarkercolor(marker,"white")
            self.graph.replot()
            return

    def closeEvent(self, event):
        if self.colormapDialog is not None:
            self.colormapDialog.close()
        if self._matplotlibSaveImage is not None:
            self._matplotlibSaveImage.close()
        qt.QWidget.closeEvent(self, event)


def test2():
    a= qt.QApplication(sys.argv)
    a.connect(a, qt.SIGNAL("lastWindowClosed()"),a,qt.SLOT("quit()"))

    w = EdfFile_StandardArray()
    w.show()
    if qt.qVersion() < '4.0.0':
        a.exec_loop()
    else:
        a.exec_()
        
def test():
    import sys
    from PyMca import EdfFileDataSource
    def repSelection(sel):    print("replaceSelection", sel)
    def removeSelection(sel): print("removeSelection", sel)
    def addSelection(sel):    print("addSelection", sel)

    a= qt.QApplication(sys.argv)
    a.connect(a, qt.SIGNAL("lastWindowClosed()"),a,qt.SLOT("quit()"))

    w = QEdfFileWidget()
    #print w
    if len(sys.argv) > 1:
        d = EdfFileDataSource.EdfFileDataSource([sys.argv[1]])
    elif os.path.exists('test.edf'):
        d = EdfFileDataSource.EdfFileDataSource(['test.edf'])
    else:
        print("Usage:")
        print("python QEdfFileWidget edffile")
        sys.exit(0)
    w.setDataSource(d)
    if QTVERSION < '4.0.0':
        qt.QObject.connect(w,qt.PYSIGNAL("addSelection"),addSelection)
        qt.QObject.connect(w,qt.PYSIGNAL("removeSelection"),removeSelection)
        qt.QObject.connect(w,qt.PYSIGNAL("replaceSelection"),repSelection)
    else:
        qt.QObject.connect(w,qt.SIGNAL("addSelection"),addSelection)
        qt.QObject.connect(w,qt.SIGNAL("removeSelection"),removeSelection)
        qt.QObject.connect(w,qt.SIGNAL("replaceSelection"),repSelection)
    w.show()
    if qt.qVersion() < '4.0.0':
        a.exec_loop()
    else:
        a.exec_()

if __name__=="__main__":
    test()
 


