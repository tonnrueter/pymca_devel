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
import sys
import os

from PyMca import QtBlissGraph
from PyMca import PyMcaQt as qt
from PyMca.PyMca_Icons import IconDict
from PyMca import PyMcaPrintPreview
from PyMca import PyMcaDirs

QTVERSION = qt.qVersion()
DEBUG = 0

class RGBCorrelatorGraph(qt.QWidget):
    def __init__(self, parent = None, selection=False, colormap=False,
                 imageicons=False, standalonesave=True, standalonezoom=True,
                 profileselection=False):
        qt.QWidget.__init__(self, parent)
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(0)
        if QTVERSION < '4.0.0':
            profileselection = False
        self._buildToolBar(selection, colormap, imageicons,
                           standalonesave,
                           standalonezoom=standalonezoom,
                           profileselection=profileselection)
        self.graph = QtBlissGraph.QtBlissGraph(self)
        self.graph.xlabel("Column")
        self.graph.ylabel("Row")
        self.graph.yAutoScale = 1
        self.graph.xAutoScale = 1
        if profileselection:
            if len(self._pickerSelectionButtons):
                self.connect(self.graph,
                         qt.SIGNAL('PolygonSignal'),
                         self._graphPolygonSignalReceived)
                self.connect(self._pickerSelectionWidthValue,
                             qt.SIGNAL('valueChanged(int)'),
                             self.setPickerSelectionWith)
        
        self.saveDirectory = os.getcwd()
        self.mainLayout.addWidget(self.graph)
        self.printPreview = PyMcaPrintPreview.PyMcaPrintPreview(modal = 0)
        if DEBUG:
            print("printPreview id = %d" % id(self.printPreview))

    def sizeHint(self):
        return qt.QSize(1.5 * qt.QWidget.sizeHint(self).width(),
                        qt.QWidget.sizeHint(self).height())

    def _buildToolBar(self, selection=False, colormap=False,
                      imageicons=False, standalonesave=True,
                      standalonezoom=True,profileselection=False):
        if QTVERSION < '4.0.0':
            if qt.qVersion() < '3.0':
                self.colormapIcon= qt.QIconSet(qt.QPixmap(IconDict["colormap16"]))
            else:
                self.colormapIcon= qt.QIconSet(qt.QPixmap(IconDict["colormap"]))
            self.selectionIcon	= qt.QIconSet(qt.QPixmap(IconDict["normal"]))
            self.zoomResetIcon	= qt.QIconSet(qt.QPixmap(IconDict["zoomreset"]))
            self.printIcon	= qt.QIconSet(qt.QPixmap(IconDict["fileprint"]))
            self.saveIcon	= qt.QIconSet(qt.QPixmap(IconDict["filesave"]))
            self.xAutoIcon	= qt.QIconSet(qt.QPixmap(IconDict["xauto"]))
            self.yAutoIcon	= qt.QIconSet(qt.QPixmap(IconDict["yauto"]))
            self.imageIcon      = qt.QIconSet(qt.QPixmap(IconDict["image"]))
            self.eraseSelectionIcon     = qt.QIconSet(qt.QPixmap(IconDict["eraseselect"]))
            self.rectSelectionIcon      = qt.QIconSet(qt.QPixmap(IconDict["boxselect"]))
            self.brushSelectionIcon     = qt.QIconSet(qt.QPixmap(IconDict["brushselect"]))
            self.brushIcon              = qt.QIconSet(qt.QPixmap(IconDict["brush"]))
            self.hFlipIcon	= qt.QIconSet(qt.QPixmap(IconDict["gioconda16mirror"]))
        else:
            self.colormapIcon   = qt.QIcon(qt.QPixmap(IconDict["colormap"]))
            self.selectionIcon	= qt.QIcon(qt.QPixmap(IconDict["normal"]))
            self.zoomResetIcon	= qt.QIcon(qt.QPixmap(IconDict["zoomreset"]))
            self.printIcon	= qt.QIcon(qt.QPixmap(IconDict["fileprint"]))
            self.saveIcon	= qt.QIcon(qt.QPixmap(IconDict["filesave"]))            
            self.xAutoIcon	= qt.QIcon(qt.QPixmap(IconDict["xauto"]))
            self.yAutoIcon	= qt.QIcon(qt.QPixmap(IconDict["yauto"]))
            self.hFlipIcon	= qt.QIcon(qt.QPixmap(IconDict["gioconda16mirror"]))
            self.imageIcon     = qt.QIcon(qt.QPixmap(IconDict["image"]))
            self.eraseSelectionIcon = qt.QIcon(qt.QPixmap(IconDict["eraseselect"]))
            self.rectSelectionIcon  = qt.QIcon(qt.QPixmap(IconDict["boxselect"]))
            self.brushSelectionIcon = qt.QIcon(qt.QPixmap(IconDict["brushselect"]))
            self.brushIcon          = qt.QIcon(qt.QPixmap(IconDict["brush"]))
            self.additionalIcon     = qt.QIcon(qt.QPixmap(IconDict["additionalselect"]))
            self.hLineIcon     = qt.QIcon(qt.QPixmap(IconDict["horizontal"]))
            self.vLineIcon     = qt.QIcon(qt.QPixmap(IconDict["vertical"]))
            self.lineIcon     = qt.QIcon(qt.QPixmap(IconDict["diagonal"]))

        self.toolBar = qt.QWidget(self)
        self.toolBarLayout = qt.QHBoxLayout(self.toolBar)
        self.toolBarLayout.setMargin(0)
        self.toolBarLayout.setSpacing(0)
        self.mainLayout.addWidget(self.toolBar)
        #Autoscale
        if standalonezoom:
            tb = self._addToolButton(self.zoomResetIcon,
                            self._zoomReset,
                            'Auto-Scale the Graph')
        else:
            tb = self._addToolButton(self.zoomResetIcon,
                            None,
                            'Auto-Scale the Graph')
        self.zoomResetToolButton = tb
        #y Autoscale
        tb = self._addToolButton(self.yAutoIcon,
                            self._yAutoScaleToggle,
                            'Toggle Autoscale Y Axis (On/Off)',
                            toggle = True, state=True)
        if qt.qVersion() < '4.0.0':
            tb.setState(qt.QButton.On)
        else:
            tb.setDown(True)
        self.yAutoScaleToolButton = tb
        tb.setDown(True)

        #x Autoscale
        tb = self._addToolButton(self.xAutoIcon,
                            self._xAutoScaleToggle,
                            'Toggle Autoscale X Axis (On/Off)',
                            toggle = True, state=True)
        self.xAutoScaleToolButton = tb
        tb.setDown(True)

        #colormap
        if colormap:
            tb = self._addToolButton(self.colormapIcon,
                                     None,
                                    'Change Colormap')
            self.colormapToolButton = tb

        #flip
        tb = self._addToolButton(self.hFlipIcon,
                                 None,
                                 'Flip Horizontal')
        self.hFlipToolButton = tb


        #save
        if standalonesave:
            tb = self._addToolButton(self.saveIcon,
                                 self._saveIconSignal,
                                 'Save Graph')
        else:
            tb = self._addToolButton(self.saveIcon,
                                 None,
                                 'Save')
        self.saveToolButton = tb

        #Selection
        if selection:
            tb = self._addToolButton(self.selectionIcon,
                                None,
                                'Toggle Selection Mode',
                                toggle = True,
                                state = False)
            if qt.qVersion() < '4.0.0':
                tb.setState(qt.QButton.Off)
            else:
                tb.setDown(False)
            self.selectionToolButton = tb
        #image selection icons
        if imageicons:
            tb = self._addToolButton(self.imageIcon,
                                     None,
                                     'Reset')
            self.imageToolButton = tb

            tb = self._addToolButton(self.eraseSelectionIcon,
                                     None,
                                     'Erase Selection')
            self.eraseSelectionToolButton = tb

            tb = self._addToolButton(self.rectSelectionIcon,
                                     None,
                                     'Rectangular Selection')
            self.rectSelectionToolButton = tb

            tb = self._addToolButton(self.brushSelectionIcon,
                                     None,
                                     'Brush Selection')
            self.brushSelectionToolButton = tb

            tb = self._addToolButton(self.brushIcon,
                                     None,
                                     'Select Brush')
            self.brushToolButton = tb
            if QTVERSION > '4.0.0':
                tb = self._addToolButton(self.additionalIcon,
                                     None,
                                     'Additional Selections Menu')
                self.additionalSelectionToolButton = tb            
        else:
            self.imageToolButton = None

        #picker selection
        self._pickerSelectionButtons = []
        if profileselection:
            self._profileSelection = True
            self._polygonSelection = False
            self._pickerSelectionButtons = []
            if self._profileSelection:
                tb = self._addToolButton(self.hLineIcon,
                                     self._hLineProfileClicked,
                                     'Horizontal Profile Selection',
                                     toggle=True,
                                     state=False)
                self.hLineProfileButton = tb
                self._pickerSelectionButtons.append(tb)
    
                tb = self._addToolButton(self.vLineIcon,
                                     self._vLineProfileClicked,
                                     'Vertical Profile Selection',
                                     toggle=True,
                                     state=False)
                self.vLineProfileButton = tb
                self._pickerSelectionButtons.append(tb)

                tb = self._addToolButton(self.lineIcon,
                                     self._lineProfileClicked,
                                     'Line Profile Selection',
                                     toggle=True,
                                     state=False)
                self.lineProfileButton = tb
                self._pickerSelectionButtons.append(tb)

                self._pickerSelectionWidthLabel = qt.QLabel(self.toolBar)
                self._pickerSelectionWidthLabel.setText("W:")
                self.toolBar.layout().addWidget(self._pickerSelectionWidthLabel)
                self._pickerSelectionWidthValue = qt.QSpinBox(self.toolBar)
                self._pickerSelectionWidthValue.setMinimum(1)
                self._pickerSelectionWidthValue.setMaximum(1000)
                self.toolBar.layout().addWidget(self._pickerSelectionWidthValue)
                #tb = self._addToolButton(None,
                #                     self._lineProfileClicked,
                #                     'Line Profile Selection',
                #                     toggle=True,
                #                     state=False)
                #tb.setText = "W:"
                #self.lineWidthProfileButton = tb
                #self._pickerSelectionButtons.append(tb)
            if self._polygonSelection:
                print("Polygon selection not implemented yet")
        #hide profile selection buttons
        if imageicons:
            for button in self._pickerSelectionButtons:
                button.hide()

        self.infoWidget = qt.QWidget(self.toolBar)
        self.infoWidget.mainLayout = qt.QHBoxLayout(self.infoWidget)
        self.infoWidget.mainLayout.setMargin(0)
        self.infoWidget.mainLayout.setSpacing(0)
        self.infoWidget.label = qt.QLabel(self.infoWidget)
        self.infoWidget.label.setText("X = ???? Y = ???? Z = ????")
        self.infoWidget.mainLayout.addWidget(self.infoWidget.label)
        self.toolBarLayout.addWidget(self.infoWidget)
        self.infoWidget.hide()

        self.toolBarLayout.addWidget(qt.HorizontalSpacer(self.toolBar))

        # ---print
        tb = self._addToolButton(self.printIcon,
                                 self.printGraph,
                                 'Prints the Graph')

    def showInfo(self):
        self.infoWidget.show()

    def hideInfo(self):
        self.infoWidget.hide()

    def setInfoText(self, text):
        self.infoWidget.label.setText(text)

    def infoText(self):
        return self.infoWidget.label.text()

    def setXLabel(self, label="Column"):
        return self.graph.x1Label(label)

    def setYLabel(self, label="Row"):
        return self.graph.ylabel(label)

    def getXLabel(self):
        return self.graph.x1Label()

    def getYLabel(self):
        return self.graph.y1Label()

    def hideImageIcons(self):
        if self.imageToolButton is None:return
        self.imageToolButton.hide()
        self.eraseSelectionToolButton.hide()
        self.rectSelectionToolButton.hide()
        self.brushSelectionToolButton.hide()
        self.brushToolButton.hide()
        if QTVERSION > '4.0.0':
            self.additionalSelectionToolButton.hide()

    def showImageIcons(self):
        if self.imageToolButton is None:return
        self.imageToolButton.show()
        self.eraseSelectionToolButton.show()
        self.rectSelectionToolButton.show()
        self.brushSelectionToolButton.show()
        self.brushToolButton.show()
        if QTVERSION > '4.0.0':
            self.additionalSelectionToolButton.show()

    def _hLineProfileClicked(self):
        for button in self._pickerSelectionButtons:
            if button != self.hLineProfileButton:
                button.setChecked(False)

        if self.hLineProfileButton.isChecked():
            self._setPickerSelectionMode("HORIZONTAL")
        else:
            self._setPickerSelectionMode(None)

    def _vLineProfileClicked(self):
        for button in self._pickerSelectionButtons:
            if button != self.vLineProfileButton:
                button.setChecked(False)
        if self.vLineProfileButton.isChecked():
            self._setPickerSelectionMode("VERTICAL")
        else:
            self._setPickerSelectionMode(None)

    def _lineProfileClicked(self):
        for button in self._pickerSelectionButtons:
            if button != self.lineProfileButton:
                button.setChecked(False)
        if self.lineProfileButton.isChecked():
            self._setPickerSelectionMode("LINE")
        else:
            self._setPickerSelectionMode(None)

    def setPickerSelectionWith(self, intValue):
        self._pickerSelectionWidthValue.setValue(intValue)
        #get the current mode
        mode = "NONE"
        for button in self._pickerSelectionButtons:
            if button.isChecked():
                if button == self.hLineProfileButton:
                    mode = "HORIZONTAL"
                elif button == self.vLineProfileButton:
                    mode = "VERTICAL"
                elif button == self.lineProfileButton:
                    mode = "LINE"        
        ddict = {}
        ddict['event'] = "PolygonWidthChanged"
        ddict['pixelwidth'] = self._pickerSelectionWidthValue.value()
        ddict['mode'] = mode
        self.emit(qt.SIGNAL('PolygonSignal'), ddict)

    def hideProfileSelectionIcons(self):
        if not len(self._pickerSelectionButtons):
            return
        for button in self._pickerSelectionButtons:
            button.setChecked(False)
            button.hide()
        self._pickerSelectionWidthLabel.hide()
        self._pickerSelectionWidthValue.hide()
        self.graph.setPickerSelectionModeOff()

    def showProfileSelectionIcons(self):
        if not len(self._pickerSelectionButtons):
            return
        for button in self._pickerSelectionButtons:
            button.show()
        self._pickerSelectionWidthLabel.show()
        self._pickerSelectionWidthValue.show()

    def _setPickerSelectionMode(self, mode=None):
        if mode is None:
            self.graph.setPickerSelectionModeOff()
            self.graph.enableZoom(True)
        else:
            try:
                self.graph.enableZoom(False)
                self.graph.setPickerSelectionModeOn(mode)
            except:
                self.graph.enableZoom(True)
                qt.QMessageBox.critical(self, "Cannot set picker mode %s" % mode,
                                        "%s" % sys.exc_info()[1])
                if DEBUG:
                    raise
        ddict = {}
        if mode is None:
            mode = "NONE"
        ddict['event'] = "PolygonModeChanged"
        ddict['mode'] = mode
        self.emit(qt.SIGNAL('PolygonSignal'), ddict)

    def _graphPolygonSignalReceived(self, ddict):
        if DEBUG:
            print("PolygonSignal Received")
            for key in ddict.keys():
                print(key, ddict[key])
        ddict['pixelwidth'] = self._pickerSelectionWidthValue.value()
        self.emit(qt.SIGNAL('PolygonSignal'), ddict)


    def _addToolButton(self, icon, action, tip, toggle=None, state=None, position=None):
        tb      = qt.QToolButton(self.toolBar)            
        if QTVERSION < '4.0.0':
            tb.setIconSet(icon)
            qt.QToolTip.add(tb,tip) 
            if toggle is not None:
                if toggle:
                    tb.setToggleButton(1)
                    if state is not None:
                        if state:
                            tb.setState(qt.QButton.On)
        else:
            if icon is not None:
                tb.setIcon(icon)
            tb.setToolTip(tip)
            if toggle is not None:
                if toggle:
                    tb.setCheckable(1)
                    if state is not None:
                        if state:
                            tb.setChecked(state)
                    else:
                        tb.setChecked(False)
        if position is not None:
            self.toolBarLayout.insertWidget(position, tb)
        else:
            self.toolBarLayout.addWidget(tb)
        if action is not None:
            self.connect(tb,qt.SIGNAL('clicked()'), action)
        return tb

    def _zoomReset(self, replot=None):
        if DEBUG:
            print("_zoomReset")
        if replot is None:
            replot = True
        if self.graph is not None:
            self.graph.zoomReset()
            if self.graph.yAutoScale:
                if hasattr(self, '_y1Limit'):
                    self.graph.sety1axislimits(0, self._y1Limit)
            if self.graph.xAutoScale:
                if hasattr(self, '_x1Limit'):
                    self.graph.setx1axislimits(0, self._x1Limit)
            if replot:
                self.graph.replot()

    def _yAutoScaleToggle(self):
        if self.graph is not None:
            if self.graph.yAutoScale:
                self.graph.yAutoScale = False
                self.yAutoScaleToolButton.setDown(False)
            else:
                self.graph.yAutoScale = True
                self.yAutoScaleToolButton.setDown(True)
            
    def _xAutoScaleToggle(self):
        if self.graph is not None:
            if self.graph.xAutoScale:
                self.graph.xAutoScale = False
                self.xAutoScaleToolButton.setDown(False)
            else:
                self.graph.xAutoScale = True
                self.xAutoScaleToolButton.setDown(True)

    def _saveIconSignal(self):
        self.saveDirectory = PyMcaDirs.outputDir

        fileTypeList = ["Image *.png",
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
            outfile.setDir(self.saveDirectory)
            ret = outfile.exec_loop()
        else:
            outfile.setWindowTitle("Output File Selection")
            if hasattr(qt, "QStringList"):
                strlist = qt.QStringList()
            else:
                strlist = []
            for f in fileTypeList:
                strlist.append(f)
            outfile.setFilters(strlist)
            outfile.setFileMode(outfile.AnyFile)
            outfile.setAcceptMode(qt.QFileDialog.AcceptSave)
            outfile.setDirectory(self.saveDirectory)
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
            outputFile = outstr
        outputDir  = os.path.dirname(outstr)
        self.saveDirectory = outputDir
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
                qt.QMessageBox.critical(self, "Save Error",
                                        "Cannot overwrite existing file")
                return

        if filetype.upper() == "IMAGE":
            self.saveGraphImage(outputFile, original = True)
        elif filetype.upper() == "ZOOMEDIMAGE":
            self.saveGraphImage(outputFile, original = False)
        else:
            self.saveGraphWidget(outputFile)

    def saveGraphImage(self, filename, original = False):
        format_ = filename[-3:].upper()
        if original:
            #This is the whole image, not the zoomed one ...
            if QTVERSION < '4.0.0':
                pixmap = qt.QPixmap(self.graph.plotImage.image)
            else:
                pixmap = qt.QPixmap.fromImage(self.graph.plotImage.image)
        else:
            pixmap = qt.QPixmap.grabWidget(self.graph.canvas())
        if pixmap.save(filename, format_):
            return
        else:
            qt.QMessageBox.critical(self, "Save Error",
                                    "%s" % sys.exc_info()[1])
            return

    def saveGraphWidget(self, filename):
        format_ = filename[-3:].upper()
        pixmap = qt.QPixmap.grabWidget(self.graph)
        if pixmap.save(filename, format_):
            return
        else:
            qt.QMessageBox.critical(self, "Save Error", "%s" % sys.exc_info()[1])
            return

    def setSaveDirectory(self, wdir):
        if os.path.exists(wdir):
            self.saveDirectory = wdir
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

    def selectColormap(self):
        qt.QMessageBox.information(self, "Open", "Not implemented (yet)")  


class MyQLabel(qt.QLabel):
    def __init__(self,parent=None,name=None,fl=0,bold=True, color= qt.Qt.red):
        qt.QLabel.__init__(self,parent)
        if qt.qVersion() <'4.0.0':
            self.color = color
            self.bold  = bold
        else:
            palette = self.palette()
            role = self.foregroundRole()
            palette.setColor(role,color)
            self.setPalette(palette)
            self.font().setBold(bold)


    if qt.qVersion() < '4.0.0':
        def drawContents(self, painter):
            painter.font().setBold(self.bold)
            pal =self.palette()
            pal.setColor(qt.QColorGroup.Foreground,self.color)
            self.setPalette(pal)
            qt.QLabel.drawContents(self,painter)
            painter.font().setBold(0)

def test():
    app = qt.QApplication([])
    qt.QObject.connect(app,
                       qt.SIGNAL("lastWindowClosed()"),
                       app,
                       qt.SLOT('quit()'))

    container = RGBCorrelatorGraph()
    container.show()
    if QTVERSION < '4.0.0':
        app.setMainWidget(container)
        app.exec_loop()
    else:
        app.exec_()

if __name__ == "__main__":
    test()
        
