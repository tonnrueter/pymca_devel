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
__revision__ = "$Revision: 1.3 $"
import sys
from PyMca import spswrap as sps
from PyMca import PyMcaQt as qt
from PyMca import SpecFileCntTable
from PyMca import MaskImageWidget
QTVERSION = qt.qVersion()
import PyMca.PyMca_Icons as icons
DEBUG = 0
SOURCE_TYPE = 'SPS'
SCAN_MODE = True

if QTVERSION > '4.0.0':
    class QGridLayout(qt.QGridLayout):
        def addMultiCellWidget(self, w, r0, r1, c0, c1, *var):
            self.addWidget(w, r0, c0, 1 + r1 - r0, 1 + c1 - c0)

class SPSFramesMcaWidget(qt.QWidget):
    def __init__(self, parent=None):
        qt.QWidget.__init__(self, parent)
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(0)
        self.graphWidget = MaskImageWidget.MaskImageWidget(self,
                                    imageicons=False,
                                    selection=False)
        self.graph = self.graphWidget.graphWidget.graph 
        self.mainLayout.addWidget(self.graphWidget)

    def setInfo(self, info):
        self.setDataSize(info["rows"], info["cols"])
        self.setTitle(info["Key"])
        self.info=info

    def setDataSource(self, data):
        self.data = data
        self.connect(self.data,
                     qt.SIGNAL("updated"),
                     self._update)
        
        dataObject = self._getDataObject()
        self.graphWidget.setImageData(dataObject.data)
        self.lastDataObject = dataObject

    def _update(self, ddict):
        targetwidgetid = ddict.get('targetwidgetid', None)
        if targetwidgetid not in [None, id(self)]:
            return
        dataObject = self._getDataObject(ddict['Key'],
                                        selection=None)
        if dataObject is not None:
            self.graphWidget.setImageData(dataObject.data)
            self.lastDataObject = dataObject

    def _getDataObject(self, key=None, selection=None):
        if key is None:
            key = self.info['Key']
        dataObject = self.data.getDataObject(key,
                                             selection=None,
                                             poll=False)
        if dataObject is not None:
            dataObject.info['legend'] = self.info['Key']
            dataObject.info['imageselection'] = False
            dataObject.info['scanselection'] = False
            dataObject.info['targetwidgetid'] = id(self)
            self.data.addToPoller(dataObject)
        return dataObject
        

    def setDataSize(self,rows,cols,selsize=None):
        self.rows= rows
        self.cols= cols
        if self.cols<=self.rows:
            self.idx='cols'
        else:
            self.idx='rows'

    def setTitle(self, title):
        self.graph.setTitle("%s"%title)

    def getSelection(self):
        keys = {"plot":self.idx,"x":0,"y":1}
        return [keys]
        

class SPSScanArrayWidget(SpecFileCntTable.SpecFileCntTable):
    def setInfo(self, info):
        if DEBUG: 
            print("info = ", info)
        if "envdict" in info:
            if len(info["envdict"].keys()):
                #We have environment information
                if "datafile" in info["envdict"]:
                    if info["envdict"]["datafile"] != "/dev/null":
                        if DEBUG:
                            print("I should send a signal, either from here or from the parent to the dispatcher")
                        if DEBUG:
                            print("SPEC data file = %s" % datafile)
                #usefull keys = ["datafile", "scantype", "axistitles","plotlist", "xlabel", "ylabel"]
                #
                #info = self.data.getKeyInfo(sel[0])
                #except:
                #    info, data = self.data.LoadSource(sel[0])
                cntList = info.get("LabelNames", [])
                ycntidx = info["envdict"].get('plotlist', "")
                if len(ycntidx):
                    ycntidx   = ycntidx.split(',')
                self.build(cntList)
                #self.cntTable.setCounterSelection(self._oldCntSelection)
                return
        if info['cols'] > 0:
            #arrayname = info['Key']
            arrayname = 'Column'
            cntList = []
            for i in range(info['cols']):
                cntList.append('%s_%03d' % (arrayname, i))
        self.build(cntList)
            
    def getSelection(self):
        #get selected counter keys
        cnt_sel = self.getCounterSelection()
        sel_list = []
        
        #build the appropriate selection for mca's
        if len(cnt_sel['cntlist']):
            if len(cnt_sel['y']): #if there is something to plot
                for index in cnt_sel['y']:
                    sel = {}
                    sel['selection'] = {}
                    sel['plot'] = 'scan'
                    sel['scanselection']  = True
                    sel['selection']['x'] = cnt_sel['x'] 
                    sel['selection']['y'] = [index] 
                    sel['selection']['m'] = cnt_sel['m']
                    sel['selection']['cntlist'] = cnt_sel['cntlist']
                    sel_list.append(sel)
        return sel_list


class SPSMcaArrayWidget(qt.QWidget):
    def __init__(self, parent=None, name="SPS_MCA_DATA", fl=0, title="MCA", size=(0,8192)):
        if QTVERSION < '4.0.0':
            qt.QWidget.__init__(self, parent, name, fl)
            layout= qt.QGridLayout(self, 5, 2)
        else:
            qt.QWidget.__init__(self, parent)
            self.setWindowTitle(name)
            layout= QGridLayout(self)
        layout.setMargin(5)
        layout.setSpacing(0)

        self.title= qt.QLabel(self)
        font= self.title.font()
        font.setBold(1)
        self.title.setFont(font)
        if QTVERSION < '4.0.0':
            layout.addMultiCellWidget(self.title, 0, 0, 0, 1, qt.Qt.AlignCenter)
            layout.addRowSpacing(0, 40)
        else:
            #layout.addMultiCellWidget(self.title, 0, 0, 0, 1, qt.Qt.AlignCenter)
            layout.addWidget(self.title, 0, 0)
            layout.setAlignment(self.title, qt.Qt.AlignCenter)
        self.setTitle(title)
        
    def setInfo(self, info):
        self.setDataSize(info["rows"], info["cols"])
        self.setTitle(info["Key"])

    def setDataSize(self,rows,cols,selsize=None):
        self.rows= rows
        self.cols= cols
        if self.cols<=self.rows:
            self.idx='cols'
        else:
            self.idx='rows'

    def setTitle(self, title):
        self.title.setText("%s"%title)

    def getSelection(self):
        keys = {"plot":self.idx,"x":0,"y":1}
        return [keys]

class SPSXiaArrayWidget(qt.QWidget):
    def __init__(self, parent=None, name="SPS_XIA_DATA", fl=0, title="XIA", size=(0,8192)):
        if QTVERSION < '4.0.0':       
            qt.QWidget.__init__(self, parent, name, fl)
            layout= qt.QGridLayout(self, 2, 1)
        else:
            qt.QWidget.__init__(self, parent)
            layout= qt.QGridLayout(self)
        layout.setMargin(5)
        layout.setSpacing(0)

        self.title= qt.QLabel(self)
        font= self.title.font()
        font.setBold(1)
        self.title.setFont(font)
        self.title.setText(title)

        if QTVERSION < '4.0.0':
            self.detList= qt.QListBox(self)
            self.detList.setSelectionMode(qt.QListBox.Multi)

            layout.addWidget(self.title, 0, 0, qt.Qt.AlignCenter)
            layout.addRowSpacing(0, 40)
            layout.addWidget(self.detList, 1, 0)
        else:
            self.detList= qt.QListWidget(self)
            self.detList.setSelectionMode(qt.QAbstractItemView.ExtendedSelection)

            layout.addWidget(self.title, 0, 0)
            layout.setAlignment(self.title, qt.Qt.AlignCenter)
            ##layout.addRowSpacing(0, 40)
            if DEBUG:
                print("row spacing")
            layout.addWidget(self.detList, 1, 0)

    def setTitle(self, title):
        self.title.setText("%s"%title)

    def setInfo(self, info):
        self.setDataSize(info["rows"], info["cols"], info.get("Detectors", None))
        self.setTitle(info["Key"])

    def setDataSize(self, rows, cols, dets=None):
        self.rows= rows
        self.cols= cols

        if dets is None or (len(dets)!=(rows-1)):
            dets= range(self.rows)

        self.detList.clear()
        if QTVERSION < '4.0.0':
            for idx in range(1, self.rows):
                self.detList.insertItem("Detector %d"%dets[idx-1])
        else:
            for idx in range(1, self.rows):
                self.detList.addItem("Detector %d"%dets[idx-1])

    def getSelection(self):
        selection= []
        if QTVERSION < '4.0.0':
            ylist= [ (idx+1) for idx in range(self.detList.count()) if self.detList.isSelected(idx) ]
        else:
            itemlist = self.detList.selectedItems()
            ylist = [int(str(item.text()).split()[-1]) for item in itemlist]
        for y in ylist:
            selection.append({"plot":"XIA", "x":0, "y":y})
        return selection

class SPS_ImageArray(qt.QWidget):
    def __init__(self, parent=None, name="SPS_ImageArray", fl=0, title="MCA", size=(0,8192)):
        if QTVERSION < '4.0.0':
            qt.QWidget.__init__(self, parent, name, fl)
            layout= qt.QGridLayout(self, 5, 2)
        else:
            qt.QWidget.__init__(self, parent)
            self.setWindowTitle(name)
            layout= QGridLayout(self)
        layout.setMargin(5)
        layout.setSpacing(0)

        self.title= qt.QLabel(self)
        font= self.title.font()
        font.setBold(1)
        self.title.setFont(font)
        if QTVERSION < '4.0.0':
            layout.addMultiCellWidget(self.title, 0, 0, 0, 1, qt.Qt.AlignCenter)
            layout.addRowSpacing(0, 40)
        else:
            #layout.addMultiCellWidget(self.title, 0, 0, 0, 1, qt.Qt.AlignCenter)
            layout.addWidget(self.title, 0, 0)
            layout.setAlignment(self.title, qt.Qt.AlignCenter)
        self.setTitle(title)
        
    def setInfo(self, info):
        self.setDataSize(info["rows"], info["cols"])
        self.setTitle(info["Key"])

    def setDataSize(self,rows,cols,selsize=None):
        self.rows= rows
        self.cols= cols

    def setTitle(self, title):
        self.title.setText("%s"%title)

    def getSelection(self):
        #get selected counter keys
        sel_list = []
        sel = {}
        sel['selection'] = {}
        sel['plot'] = 'image'
        sel['scanselection']  = False
        sel['selection'] = None
        sel_list.append(sel)
        return sel_list
 
class SPS_StandardArray(qt.QWidget):
    def __init__(self, parent=None, name="SPS_StandardArray", fl=0, rows=0, cols=0):
        if QTVERSION < '4.0.0':
            qt.QWidget.__init__(self, parent, name, fl)
            layout= qt.QGridLayout(self, 4, 2)
        else:
            qt.QWidget.__init__(self, parent)
            layout= qt.QGridLayout(self)
        layout.setMargin(5)
        layout.setSpacing(0)

        plab= qt.QLabel("Plot", self)
        xlab= qt.QLabel("X :", self)
        ylab= qt.QLabel("Y :", self)

        layout.addWidget(plab, 0, 0, qt.Qt.AlignRight)
        layout.addWidget(xlab, 1, 0, qt.Qt.AlignRight)
        layout.addWidget(ylab, 2, 0, qt.Qt.AlignRight|qt.Qt.AlignTop)

        self.plotCombo= qt.QComboBox(self)
        self.plotCombo.setEditable(0)
        if QTVERSION < '4.0.0':
            self.plotCombo.insertItem("Rows")
            self.plotCombo.insertItem("Columns")
        else:
            self.plotCombo.addItem("Rows")
            self.plotCombo.addItem("Columns")

        self.xCombo= qt.QComboBox(self)
        self.xCombo.setEditable(0)

        if QTVERSION < '4.0.0':
            self.yList= qt.QListBox(self)
            self.yList.setSelectionMode(qt.QListBox.Multi)
        else:
            self.yList= qt.QListWidget(self)
            self.yList.setSelectionMode(qt.QAbstractItemView.ExtendedSelection)

        layout.addWidget(self.plotCombo, 0, 1)
        layout.addWidget(self.xCombo, 1, 1)
        layout.addWidget(self.yList, 2, 1)

        self.connect(self.plotCombo, qt.SIGNAL("activated(int)"), self.__plotChanged)

        self.setDataSize(rows, cols)

    def setDataSize(self, rows, cols):
        self.rows= rows
        self.cols= cols

        idx= self.cols<=self.rows
        if QTVERSION < '4.0.0':
            self.plotCombo.setCurrentItem(idx)
        else:
            self.plotCombo.setCurrentIndex(idx)
        self.__plotChanged(idx)

    def __plotChanged(self, index):
        if index==1:
            txt= "Column"
            val= self.cols
        else:
            txt= "Row"
            val= self.rows
        self.xCombo.clear()
        if QTVERSION < '4.0.0':
            self.xCombo.insertItem("Array Index")
            self.yList.clear()
            for x in range(val):
                self.xCombo.insertItem("%s %d"%(txt,x))
                self.yList.insertItem("%s %d"%(txt,x))
            if val==2:
                self.xCombo.setCurrentItem(0)
                self.__xChanged(0)
        else:
            self.xCombo.addItem("Array Index")
            self.yList.clear()
            for x in range(val):
                self.xCombo.addItem("%s %d"%(txt,x))
                self.yList.addItem("%s %d"%(txt,x))
            if val==2:
                self.xCombo.setCurrentIndex(0)
                self.__xChanged(0)

    def __xChanged(self, index):
        pass

    def getSelection(self):
        selection= []
        if QTVERSION < '4.0.0':
            idx= self.plotCombo.currentItem()
        else:
            idx= self.plotCombo.currentIndex()
        if idx==1: plot= "cols"
        else: plot= "rows"

        if QTVERSION < '4.0.0':
            idx= self.xCombo.currentItem()
        else:
            idx= self.xCombo.currentIndex()
        if idx==0: x= None
        else: x= idx-1

        if QTVERSION < '4.0.0':
            ylist= [ idx for idx in range(self.yList.count()) if self.yList.isSelected(idx) ]
        else:
            itemlist = self.yList.selectedItems()
            ylist = [int(str(item.text()).split()[-1]) for item in itemlist]
        for y in ylist:
            selection.append({"plot":plot, "x":x, "y":y})
        return selection


class QSpsWidget(qt.QWidget):
    HiddenArrays= ["MCA_DATA_PARAM", "XIA_STAT", "XIA_DET"]
    WidgetArrays= {"scan":SPSScanArrayWidget, 
                   "xia": SPSXiaArrayWidget, 
                   "mca": SPSMcaArrayWidget, 
                   "array": SPS_StandardArray,
                   "image": SPS_ImageArray,
                   "frames_mca":SPSFramesMcaWidget,
                   "frames_image":qt.QWidget,
                   "empty": qt.QWidget}
    TypeArrays= {"MCA_DATA": "mca", "XIA_PLOT": "mca",
                 "XIA_DATA": "xia", "XIA_BASELINE":"xia",
                 "SCAN_D": "scan", "image_data":"image" }

    def __init__(self, parent=None, name="SPSSelector", fl=0):
        if QTVERSION < '4.0.0':
            qt.QWidget.__init__(self, parent, name, fl)
        else:
            qt.QWidget.__init__(self, parent)

        self.dataSource= None

        self.data= None
        self.currentSpec= None
        self.currentArray= None
        self.selection= None
        self.openFile = self.refreshSpecList

        self.selectPixmap= qt.QPixmap(icons.selected)
        self.unselectPixamp= qt.QPixmap(icons.unselected)

        mainLayout= qt.QVBoxLayout(self)

        # --- spec name selection
        specWidget= qt.QWidget(self)
        self.specCombo= qt.QComboBox(specWidget)
        self.specCombo.setEditable(0)
        if QTVERSION < '4.0.0':
            self.reload= qt.QIconSet(qt.QPixmap(icons.reload))
            refreshButton= qt.QToolButton(specWidget)
            refreshButton.setIconSet(self.reload)
            self.closeIcon= qt.QIconSet(qt.QPixmap(icons.fileclose))
            closeButton= qt.QToolButton(specWidget)
            closeButton.setIconSet(self.closeIcon)
        else:
            self.reload= qt.QIcon(qt.QPixmap(icons.reload))
            refreshButton= qt.QToolButton(specWidget)
            refreshButton.setIcon(self.reload)
            self.closeIcon= qt.QIcon(qt.QPixmap(icons.fileclose))
            closeButton= qt.QToolButton(specWidget)
            closeButton.setIcon(self.closeIcon)
        refreshButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        closeButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        specLayout= qt.QHBoxLayout(specWidget)
        specLayout.addWidget(self.specCombo)
        specLayout.addWidget(refreshButton)
        specLayout.addWidget(closeButton)

        self.connect(refreshButton, qt.SIGNAL("clicked()"), self.refreshSpecList)
        self.connect(closeButton,qt.SIGNAL("clicked()"), self.closeCurrentSpec)
        self.connect(self.specCombo, qt.SIGNAL("activated(const QString &)"), 
                    self.refreshArrayList)

        # --- splitter
        self.splitter= qt.QSplitter(self)
        if QTVERSION < '4.0.0':
            self.splitter.setOrientation(qt.QSplitter.Vertical)
        else:
            self.splitter.setOrientation(qt.Qt.Vertical)

        # --- shm array list
        if QTVERSION < '4.0.0':
            self.arrayList= qt.QListView(self.splitter, "ShmArrayList")
            self.arrayList.addColumn("")
            self.arrayList.addColumn("Array Name")
            self.arrayList.addColumn("Rows")
            self.arrayList.addColumn("Cols")
            self.arrayList.setSorting(-1)
            self.arrayList.header().setClickEnabled(0,-1)
            self.arrayList.setAllColumnsShowFocus(1)
            self.arrayList.setSelectionMode(qt.QListView.Single)
            self.connect(self.arrayList,
                         qt.SIGNAL("selectionChanged()"),
                         self.__arraySelection)
        else:
            self.arrayList= qt.QTreeWidget(self.splitter)
            labels = ["","Array Name", "Rows","Cols"]
            self.arrayList.setColumnCount(len(labels))
            self.arrayList.setHeaderLabels(labels)
            self.arrayList.setSelectionMode(qt.QAbstractItemView.SingleSelection)
            self.connect(self.arrayList,
                         qt.SIGNAL("itemSelectionChanged()"),
                         self.__arraySelection)

        # --- array parameter
        if QTVERSION < '4.0.0':
            self.paramIndex= {}
            self.paramWidget= qt.QWidgetStack(self.splitter)
            for wtype in self.WidgetArrays.keys():
                widclass= self.WidgetArrays[wtype]
                wid= widclass(self.paramWidget)
                self.paramWidget.addWidget(wid)
                self.paramIndex[wtype]= self.paramWidget.id(wid)
        else:
            self.paramIndex= {}
            self.paramWidget= qt.QStackedWidget(self.splitter)
            for wtype in self.WidgetArrays.keys():
                widclass= self.WidgetArrays[wtype]
                wid= widclass(self.paramWidget)
                self.paramWidget.addWidget(wid)
                self.paramIndex[wtype]= self.paramWidget.indexOf(wid)

        # --- command buttons
        butWidget= qt.QWidget(self)
        butWidget.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Minimum,
                                               qt.QSizePolicy.Minimum))
        addButton= qt.QPushButton("Add", butWidget)
        removeButton= qt.QPushButton("Remove", butWidget)
        replaceButton= qt.QPushButton("Replace", butWidget)

        butLayout= qt.QHBoxLayout(butWidget)
        butLayout.addWidget(addButton)
        butLayout.addWidget(removeButton)
        butLayout.addWidget(replaceButton)
        butLayout.setMargin(5)

        self.connect(addButton, qt.SIGNAL("clicked()"), self.__addClicked)
        self.connect(replaceButton, qt.SIGNAL("clicked()"), self.__replaceClicked)
        self.connect(removeButton, qt.SIGNAL("clicked()"), self.__removeClicked)

        # --- main layout
        mainLayout.setMargin(5)
        mainLayout.setSpacing(5)
        mainLayout.addWidget(specWidget)
        if __name__ != "__main__":
            specWidget.hide()
        mainLayout.addWidget(self.splitter)
        mainLayout.addWidget(butWidget)

    def setData(self,data=None):
        if DEBUG:
            print("setData(self, data) called")
            print("spec data = ",data)
        self.data= data
        self.refreshSpecList()
        self.refreshDataSelection()

    def setDataSource(self,data=None):
        if DEBUG:
            print("setDataSource(self, data) called")
            print("spec data = ",data)
        self.data= data
        self.refreshSpecList()
        self.refreshDataSelection()
        if data is None:
            self.arrayList.clear()
        else:
            self.refreshArrayList(data.sourceName)

    def refreshSpecList(self):
        speclist= sps.getspeclist()
        if self.specCombo.count():
            selected= str(self.specCombo.currentText())
        else:    selected= None
        self.specCombo.clear()
        if len(speclist):
            for spec in speclist:
                if QTVERSION < '4.0.0':
                    self.specCombo.insertItem(spec)
                else:
                    self.specCombo.addItem(spec)
            self.selectSpec(selected or speclist[0])

    def selectSpec(self, specname=None):
        if QTVERSION < '4.0.0':
            for idx in range(self.specCombo.count()):
                if str(self.specCombo.text(idx))==specname:
                    self.specCombo.setCurrentItem(idx)
        else:
            for idx in range(self.specCombo.count()):
                if str(self.specCombo.itemText(idx))==specname:
                    self.specCombo.setCurrentIndex(idx)

    def __getCurrentSpec(self):
        if self.specCombo.count():
            return str(self.specCombo.currentText())
        else:    return None

    def refreshDataSelection(self, source=None):
        spec= self.__getCurrentSpec()
        if spec is not None and self.dataSource is not None:
            arraylist= self.dataSource.GetDataList(spec)
            item= self.arrayList.firstChild()
            while item is not None:
                name= str(item.text(1))
                if name in arraylist:    item.setPixmap(0, self.selectPixmap)
                else:            item.setPixmap(0, self.unselectPixmap)
            item= item.nextSibling()

    def closeCurrentSpec(self):
        spec= self.__getCurrentSpec()
        if spec is not None and self.dataSource is not None:
            arraylist= self.DataSource.GetDataList(spec)
            if len(arraylist):
                msg= "%d spectrums are linked to that SPEC source.\n"%(len(arraylist))
                msg+= "Do you really want to delete all these spectrums ??"
                ans= qt.QMessageBox.information(self, "Remove Spec Shared %s"%spec, msg, \
                        qt.QMessageBox.No, qt.QMessageBox.Yes)
                if ans.qt.QMessageBox.Yes:
                    self.dataSource.RemoveData(spec)

    def refreshArrayList(self,qstring):
        self.arrayList.clear()
        #spec= self.__getCurrentSpec()
        self.currentSpec = str(qstring)
        spec = self.currentSpec
        if spec is not None:
            arraylist= {}
            for array in sps.getarraylist(spec):
                if array not in self.HiddenArrays:
                    info= sps.getarrayinfo(spec, array)
                    rows= info[0]
                    cols= info[1]
                    type= info[2]
                    flag= info[3]
                    if DEBUG:
                        print(" array = ", array)
                        print(" flag = ", flag)
                        print(" type = ", type)
                    if type!=sps.STRING:
                        if (flag & sps.TAG_ARRAY) == sps.TAG_ARRAY:
                            arraylist[array]= (rows, cols)
            if len(arraylist.keys()):
                arrayorder= arraylist.keys()
                arrayorder.sort()
                arrayorder.reverse()
                if QTVERSION < '4.0.0':
                    for name in arrayorder:
                        self.arrayList.insertItem(qt.QListViewItem(self.arrayList, 
                            "", name, str(arraylist[name][0]), str(arraylist[name][1])))
                else:
                    for name in arrayorder:
                        item = (qt.QTreeWidgetItem(self.arrayList, 
                            ["", name, str(arraylist[name][0]), str(arraylist[name][1])]))
                self.refreshDataSelection()
            
        self.__getParamWidget("empty")

    def __arraySelection(self):
        if QTVERSION < '4.0.0':
            item= self.arrayList.selectedItem()
            if item is not None:
                self.currentArray= str(item.text(1))
            else:
                #click on empty space
                return
        else:
            item= self.arrayList.selectedItems()
            if len(item):
                item = item[0]
                self.currentArray= str(item.text(1))
            else:
                #click on empty space
                return
        #self.data.SetSource(self.currentSpec)
        #self.data.LoadSource(self.currentArray)
        info= self.data.getKeyInfo(self.currentArray)
        wid= None
        atype = None
        if 0 and ((info['flag'] & sps.TAG_FRAMES) == sps.TAG_FRAMES) and\
           ((info['flag'] & sps.TAG_IMAGE) == sps.TAG_IMAGE):
            atype = "frames_image"
        elif ((info['flag'] & sps.TAG_FRAMES) == sps.TAG_FRAMES) and\
           ((info['flag'] & sps.TAG_MCA) == sps.TAG_MCA):
            atype = "frames_mca"
        elif (info['flag'] & sps.TAG_IMAGE) == sps.TAG_IMAGE:
            atype = "image"
        elif (info['flag'] & sps.TAG_MCA)  == sps.TAG_MCA:
            atype = "mca"
        elif (info['flag'] & sps.TAG_SCAN) == sps.TAG_SCAN:
            atype = "scan"
        elif (info['rows'] > 100) and (info['cols'] > 100):
            atype = "image"
        if atype is not None:
            wid= self.__getParamWidget(atype)
            wid.setInfo(info)
            if hasattr(wid, "setDataSource"):
                wid.setDataSource(self.data)
        else:
            for (array, atype) in self.TypeArrays.items():
                if self.currentArray[0:len(array)]==array:
                    wid= self.__getParamWidget(atype)
                    wid.setInfo(info)
                    if hasattr(wid, "setDataSource"):
                        wid.setDataSource(self.data)
                    break
        if wid is None:
            arrayType = "ARRAY"
            wid= self.__getParamWidget("array")
            wid.setDataSize(info["rows"], info["cols"])
        else:
            arrayType = atype.upper()

        #emit a selection to inform about the change
        ddict = {}
        ddict['SourceName'] = self.data.sourceName
        ddict['SourceType'] = self.data.sourceType
        ddict['event'] = "SelectionTypeChanged"
        if arrayType in ["IMAGE"]:
            ddict['SelectionType'] = self.data.sourceName +" "+self.currentArray
        elif arrayType in ["MCA", "XIA"]:
            ddict['SelectionType'] = "MCA"
        elif arrayType in ["ARRAY"]:
            ddict['SelectionType'] = "MCA"
        else:
            ddict['SelectionType'] = arrayType
        if QTVERSION < '4.0.0':
            self.emit(qt.PYSIGNAL("otherSignals"), (ddict,))
        else:
            self.emit(qt.SIGNAL("otherSignals"), ddict)

    def __getParamWidget(self, widtype):
        if QTVERSION < '4.0.0':
            wid= self.paramWidget.visibleWidget()
            if self.paramWidget.id(wid)!=self.paramIndex[widtype]:
                self.paramWidget.raiseWidget(self.paramIndex[widtype])
                wid= self.paramWidget.visibleWidget()
        else:
            wid= self.paramWidget.currentWidget()
            if self.paramWidget.indexOf(wid) != self.paramIndex[widtype]:
                self.paramWidget.setCurrentIndex(self.paramIndex[widtype])
                wid = self.paramWidget.currentWidget()
        return wid

    def __replaceClicked(self):
        if DEBUG:
            print("replace clicked")
        selkeys= self.__getSelectedKeys()
        if len(selkeys):
            #self.eh.event(self.repEvent, selkeys)
            if DEBUG:
                print("Replace event")
            sel = {}
            sel['SourceType'] = SOURCE_TYPE
            sellistsignal = []
            for selection in selkeys:
                selsignal = {}
                selsignal['SourceType'] = self.data.sourceType
                selsignal['SourceName'] = self.data.sourceName
                selsignal['selection'] = None
                selsignal['Key'] = selection['Key']
                if 'SourceName' not in sel:
                    sel['SourceName'] = selection['SourceName']
                arrayname = selection['Key']
                if 'Key' not in sel:
                    sel['Key'] = selection['Key']
                if arrayname not in sel:
                    sel[arrayname] = {'rows':[],'cols':[]}
                if selection['plot'] == 'cols':
                     selsignal["selection"] = {"cols":{}}
                     selsignal["selection"]["cols"] = {}
                     selsignal["selection"]["cols"]["x"] = [selection['x']]
                     selsignal["selection"]["cols"]["y"] = [selection['y']]
                     if type(self.data.sourceName) == type(''):
                        sname = [self.data.sourceName]
                     else:
                        sname = self.data.sourceName
                     selsignal["legend"] = sname[0] +\
                                        " "+selsignal['Key']+\
                                        ".c.%d" % int(selection['y'])
                     sel[arrayname]['cols'].append({'x':selection['x'],'y':selection['y']})
                elif selection['plot'] == 'rows':
                     sel[arrayname]['rows'].append({'x':selection['x'],'y':selection['y']})                              
                     selsignal["selection"] = {"rows":{}}
                     selsignal["selection"]["rows"] = {}
                     selsignal["selection"]["rows"]["x"] = [selection['x']]
                     selsignal["selection"]["rows"]["y"] = [selection['y']]
                     if type(self.data.sourceName) == type(''):
                        sname = [self.data.sourceName]
                     else:
                        sname = self.data.sourceName
                     selsignal["legend"] = sname[0] +\
                                        " "+selsignal['Key']+\
                                        ".r.%d" % int(selection['y'])
                elif selection['plot'] == 'XIA':
                     sel[arrayname]['rows'].append({'x':selection['x'],
                                            'y':selection['y']})                              
                     #selsignal["Key"] += ".r.%d" % int(selection['y'])
                     selsignal["selection"] = {"rows":{}, "XIA":True}
                     selsignal["selection"]["rows"] = {}
                     selsignal["selection"]["rows"]["x"] = [selection['x']]
                     selsignal["selection"]["rows"]["y"] = [selection['y']]
                     if type(self.data.sourceName) == type(''):
                        sname = [self.data.sourceName]
                     else:
                        sname = self.data.sourceName
                     selsignal["legend"] = sname[0] +\
                                        " "+selsignal['Key']+\
                                        " #%02d" % int(selection['y'])

                elif selection['plot'] == 'scan':
                     if SCAN_MODE:
                         #sel[arrayname]['cols'].append({'x':selection['selection']['x'][0],
                         #                               'y':selection['selection']['y'][0]})
                         selsignal["selection"] = selection['selection']
                         selsignal['legend'] = self.data.sourceName + " " + \
                                              selsignal['Key']
                         selsignal['scanselection'] = True
                         #print "cheeting"
                         #selsignal['scanselection'] = False
                     else:
                         #do it as a col
                         sel[arrayname]['cols'].append({'x':selection['selection']['x'],
                                                        'y':selection['selection']['y']})
                         selsignal["selection"] = {'cols':{}}
                         selsignal["selection"]["cols"]["x"] = selection['selection']['x']
                         selsignal["selection"]["cols"]["y"] = selection['selection']['y']
                         selsignal['legend'] = self.data.sourceName + " " + \
                                               selsignal['Key']
                         selsignal['scanselection'] = True
                         #print "cheeting"
                         #selsignal['scanselection'] = False
                         
                elif selection['plot'] == 'image':
                    selsignal["selection"] = selection['selection']   
                    selsignal['legend'] = self.data.sourceName + " " + \
                                              selsignal['Key']
                    selsignal['scanselection']  = False
                    selsignal['imageselection'] = True
                sellistsignal.append(selsignal)
            self.setSelected([sel],reset=1)
            if QTVERSION < '4.0.0':
                self.emit(qt.PYSIGNAL("replaceSelection"), (sellistsignal,))
            else:
                self.emit(qt.SIGNAL("replaceSelection"), sellistsignal)

    def __addClicked(self):
        if DEBUG:
            print("select clicked")
        selkeys= self.__getSelectedKeys()
        if DEBUG:
            print("selected keys = ",selkeys )
        if len(selkeys):
            #self.eh.event(self.addEvent, selkeys)
            if DEBUG:
                print("Select event")
            sel = {}
            sel['SourceType'] = SOURCE_TYPE 
            sellistsignal = []           
            for selection in selkeys:
                selsignal = {}
                selsignal['SourceType'] = self.data.sourceType
                selsignal['SourceName'] = self.data.sourceName
                selsignal['selection'] = None
                selsignal['Key'] = selection['Key']
                if 'SourceName' not in sel:
                    sel['SourceName'] = selection['SourceName']
                arrayname = selection['Key']
                if 'Key' not in sel:
                    sel['Key'] = selection['Key']
                if arrayname not in sel:
                    sel[arrayname] = {'rows':[],'cols':[]}
                if selection['plot'] == 'XIA':
                     sel[arrayname]['rows'].append({'x':selection['x'],
                                            'y':selection['y']})                              
                     #selsignal["Key"] += ".r.%d" % int(selection['y'])
                     selsignal["selection"] = {"rows":{}, "XIA":True}
                     selsignal["selection"]["rows"] = {}
                     selsignal["selection"]["rows"]["x"] = [selection['x']]
                     selsignal["selection"]["rows"]["y"] = [selection['y']]
                     if type(self.data.sourceName) == type(''):
                        sname = [self.data.sourceName]
                     else:
                        sname = self.data.sourceName
                     selsignal["legend"] = sname[0] +\
                                        " "+selsignal['Key']+\
                                        " #%02d" % int(selection['y'])
                elif selection['plot'] == 'cols':
                     sel[arrayname]['cols'].append({'x':selection['x'],
                                                    'y':selection['y']})
                     #selsignal["Key"] += ".c.%d" % int(selection['y'])
                     selsignal["selection"] = {"cols":{}}
                     selsignal["selection"]["cols"] = {}
                     selsignal["selection"]["cols"]["x"] = [selection['x']]
                     selsignal["selection"]["cols"]["y"] = [selection['y']]
                     if type(self.data.sourceName) == type(''):
                        sname = [self.data.sourceName]
                     else:
                        sname = self.data.sourceName
                     selsignal["legend"] = sname[0] +\
                                        " "+selsignal['Key']+\
                                        ".c.%d" % int(selection['y'])
                elif selection['plot'] == 'rows':
                     sel[arrayname]['rows'].append({'x':selection['x'],
                                            'y':selection['y']})                              
                     #selsignal["Key"] += ".r.%d" % int(selection['y'])
                     selsignal["selection"] = {"rows":{}}
                     selsignal["selection"]["rows"] = {}
                     selsignal["selection"]["rows"]["x"] = [selection['x']]
                     selsignal["selection"]["rows"]["y"] = [selection['y']]
                     if type(self.data.sourceName) == type(''):
                        sname = [self.data.sourceName]
                     else:
                        sname = self.data.sourceName
                     selsignal["legend"] = sname[0] +\
                                        " "+selsignal['Key']+\
                                        ".r.%d" % int(selection['y'])
                elif selection['plot'] == 'scan':
                     if SCAN_MODE:
                         #sel[arrayname]['cols'].append({'x':selection['selection']['x'][0],
                         #                               'y':selection['selection']['y'][0]})
                         selsignal["selection"] = selection['selection']
                         selsignal['legend'] = self.data.sourceName + " " + \
                                              selsignal['Key']
                         selsignal['scanselection'] = True
                         #print "cheeting"
                         #selsignal['scanselection'] = False
                     else:
                         #do it as a col
                         sel[arrayname]['cols'].append({'x':selection['selection']['x'],
                                                        'y':selection['selection']['y']})
                         selsignal["selection"] = {'cols':{}}
                         selsignal["selection"]["cols"]["x"] = selection['selection']['x']
                         selsignal["selection"]["cols"]["y"] = selection['selection']['y']
                         selsignal['legend'] = self.data.sourceName + " " + \
                                               selsignal['Key']
                         selsignal['scanselection'] = True
                         #print "cheeting"
                         #selsignal['scanselection'] = False
                elif selection['plot'] == 'image':
                    selsignal["selection"] = selection['selection']   
                    selsignal['legend'] = self.data.sourceName + " " + \
                                              selsignal['Key']
                    selsignal['scanselection']  = False
                    selsignal['imageselection'] = True

                sellistsignal.append(selsignal)
            if self.selection is None: 
                self.setSelected([sel],reset=1)
            else:
                self.setSelected([sel],reset=0)
            if QTVERSION < '4.0.0':
                self.emit(qt.PYSIGNAL("addSelection"), (sellistsignal,))
            else:
                self.emit(qt.SIGNAL("addSelection"),  sellistsignal)
            
    def __getSelectedKeys(self):
        selkeys= []
        if QTVERSION < '4.0.0':
            parwid= self.paramWidget.visibleWidget()
        else:
            parwid= self.paramWidget.currentWidget()
        if self.currentArray is not None:
            for sel in parwid.getSelection():
                sel["SourceName"]= self.currentSpec
                sel['SourceType'] = SOURCE_TYPE            
                sel["Key"]= self.currentArray
                selkeys.append(sel)
        return selkeys

    def __removeClicked(self):
        if DEBUG:
            print("remove clicked")
        selkeys= self.__getSelectedKeys()
        if len(selkeys):
            #self.eh.event(self.delEvent, selkeys)
            if DEBUG:
                print("Remove Event")
                print("self.selection before = ",self.selection)
            returnedselection=[]
            sellistsignal = []           
            for selection in selkeys:
                selsignal = {}
                selsignal['SourceType'] = self.data.sourceType
                selsignal['SourceName'] = self.data.sourceName
                selsignal['selection'] = None
                selsignal['Key'] = selection['Key']
                sel = {}
                sel['SourceName'] = selection['SourceName']
                sel['SourceType'] = SOURCE_TYPE            
                sel['Key'] = selection['Key']
                arrayname = selection['Key']
                sel[arrayname] = {'rows':[],'cols':[]}
                if selection['plot'] == 'cols':
                     sel[arrayname]['cols'].append({'x':selection['x'],'y':selection['y']})
                     selsignal["selection"] = {"cols":{}}
                     selsignal["selection"]["cols"] = {}
                     selsignal["selection"]["cols"]["x"] = [selection['x']]
                     selsignal["selection"]["cols"]["y"] = [selection['y']]
                     if type(self.data.sourceName) == type(''):
                        sname = [self.data.sourceName]
                     else:
                        sname = self.data.sourceName
                     selsignal["legend"] = sname[0] +\
                                        " "+selsignal['Key']+\
                                        ".c.%d" % int(selection['y'])
                elif selection['plot'] == 'rows':
                     sel[arrayname]['rows'].append({'x':selection['x'],'y':selection['y']})
                     selsignal["selection"] = {"rows":{}}
                     selsignal["selection"]["rows"] = {}
                     selsignal["selection"]["rows"]["x"] = [selection['x']]
                     selsignal["selection"]["rows"]["y"] = [selection['y']]
                     if type(self.data.sourceName) == type(''):
                        sname = [self.data.sourceName]
                     else:
                        sname = self.data.sourceName
                     selsignal["legend"] = sname[0] +\
                                        " "+selsignal['Key']+\
                                        ".r.%d" % int(selection['y'])
                elif selection['plot'] == 'XIA':
                     sel[arrayname]['rows'].append({'x':selection['x'],
                                            'y':selection['y']})                              
                     #selsignal["Key"] += ".r.%d" % int(selection['y'])
                     selsignal["selection"] = {"rows":{}, "XIA":True}
                     selsignal["selection"]["rows"] = {}
                     selsignal["selection"]["rows"]["x"] = [selection['x']]
                     selsignal["selection"]["rows"]["y"] = [selection['y']]
                     if type(self.data.sourceName) == type(''):
                        sname = [self.data.sourceName]
                     else:
                        sname = self.data.sourceName
                     selsignal["legend"] = sname[0] +\
                                        " "+selsignal['Key']+\
                                        " #%02d" % int(selection['y'])
                elif selection['plot'] == 'scan':
                     #sel[arrayname]['cols'].append({'x':selection['selection']['x'][0],
                     #                               'y':selection['selection']['y'][0]})
                     selsignal["selection"] = selection['selection']
                     selsignal['legend'] = self.data.sourceName + " " + \
                                          selsignal['Key']
                     selsignal['scanselection'] = True

                elif selection['plot'] == 'image':
                    selsignal["selection"] = selection['selection']   
                    selsignal['legend'] = self.data.sourceName + " " + \
                                              selsignal['Key']
                    selsignal['scanselection']  = False
                    selsignal['imageselection'] = True

                sellistsignal.append(selsignal)
                returnedselection.append(sel)
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
                                for couple in  sel[arrayname]['cols']:
                                    if couple in  self.selection[sel['SourceName']][arrayname]['cols']:
                                        index= self.selection[sel['SourceName']][arrayname]['cols'].index(couple)
                                        del self.selection[sel['SourceName']][arrayname]['cols'][index]
                                seln = {}
                                seln['SourceName'] = sel['SourceName'] 
                                seln['SourceType'] = SOURCE_TYPE            
                                seln['Key']        = sel['Key']
                                seln[seln['Key']]  = self.selection[seln['SourceName']][seln['Key']]
                                self.setSelected([seln],reset=0)
            if QTVERSION < '4.0.0':
                self.emit(qt.PYSIGNAL("removeSelection"), (sellistsignal,))
            else:
                self.emit(qt.SIGNAL("removeSelection"), sellistsignal)
            
    def removeSelection(self,selection):
        if type(selection) != type([]):
            selection=[selection]
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
                                for couple in  sel[arrayname]['cols']:
                                    if couple in  self.selection[sel['SourceName']][arrayname]['cols']:
                                        index= self.selection[sel['SourceName']][arrayname]['cols'].index(couple)
                                        del self.selection[sel['SourceName']][arrayname]['cols'][index]
                                seln = {}
                                seln['SourceName'] = sel['SourceName'] 
                                seln['SourceType'] = SOURCE_TYPE            
                                seln['Key']        = sel['Key']
                                seln[seln['Key']]  = self.selection[seln['SourceName']][seln['Key']]
                                self.setSelected([seln],reset=0)
        if QTVERSION < '4.0.0':
            self.emit(qt.PYSIGNAL("removeSelection"), (selection,))
        else:
            self.emit(qt.SIGNAL("removeSelection"), (selection))
                             
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
            #selkey is the array name what to do if multiple array names?
            if type(sel["Key"]) == type([]):
                selkey = sel["Key"][0]
            else:
                selkey = sel["Key"]
            if specname not in self.selection:
                self.selection[specname]= {}
            if selkey not in self.selection[specname]:
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
                sel['SourceType']   = 'SPS'
                sel['Key']          = arraykey
                sel[arraykey]        = self.selection[sourcekey][arraykey]
                selection.append(sel)
        return selection

        
    def __refreshSelection(self):
        return
        if DEBUG:
            print("__refreshSelection(self) called")
            print(self.selection)
        if self.selection is not None:
            sel = self.selection.get(self.data.SourceName, {})
            selkeys = []
            for key in sel.keys():
                if (sel[key]['mca'] != []) or (sel[key]['scan']['Ycnt'] !=  []):
                    selkeys.append(key)
            if DEBUG:
                print("selected scans =",selkeys)
                print("but self.selection = ",self.selection)
                print("and self.selection.get(self.data.SourceName, {}) =",sel)
            self.scanList.markScanSelected(selkeys)
            scandict = sel.get(self.currentScan, {})
            if 'mca' in scandict:
                self.mcaTable.markMcaSelected(scandict['mca'])
            else:
                self.mcaTable.markMcaSelected([])
            if 'scan' in scandict:
                self.cntTable.markCntSelected(scandict['scan'])
            else:
                self.cntTable.markCntSelected({})
                
    def isSelectionUpdated(self,sellist):
        outsel = []
        if type(sellist) != type([]):
            sellist = [sellist]
        for ddict in  sellist:
            #for dict in selection:
                if 'SourceName' in ddict:
                    spec = ddict['SourceName']
                    if 'Key' in ddict:
                        shm  = ddict['Key']
                        if shm in ddict:
                            check = 0
                            rows = []
                            cols = []
                            if 'cols' in ddict[shm]:
                                cols = ddict[shm]['cols']
                                if len(cols):
                                    check =  1
                            if 'rows' in ddict[shm]:
                                rows = ddict[shm]['rows']
                                if len(rows):
                                    check =  1
                            if check and sps.specrunning(spec):
                                if sps.isupdated(spec,shm):
                                    outsel.append({'SourceName':spec,
                                                   'Key':shm,
                                                   shm:{'rows':rows,
                                                        'cols':cols},
                                                    'SourceType':'SPS'})
        return outsel                                    
                                        
                                        
        

def test():
    import sys
    if 0:
        from PyMca import SPSLayer
    else:
        from PyMca import QSpsDataSource

    a= qt.QApplication(sys.argv)
    a.connect(a, qt.SIGNAL("lastWindowClosed()"),a,qt.SLOT("quit()"))
    def repSelection(sel):
        print("replaceSelection", sel)
    def removeSelection(sel):
        print("removeSelection", sel)
    def addSelection(sel):
        print("addSelection", sel)

    w= QSpsWidget()
    if QTVERSION < '4.0.0':
        qt.QObject.connect(w,qt.PYSIGNAL("addSelection"),addSelection)
        qt.QObject.connect(w,qt.PYSIGNAL("removeSelection"),removeSelection)
        qt.QObject.connect(w,qt.PYSIGNAL("replaceSelection"),repSelection)
    else:
        qt.QObject.connect(w,qt.SIGNAL("addSelection"),addSelection)
        qt.QObject.connect(w,qt.SIGNAL("removeSelection"),removeSelection)
        qt.QObject.connect(w,qt.SIGNAL("replaceSelection"),repSelection)
    #d = QSpsDataSource.QSpsDataSource()
    #w.setData(d)
    """
    w.eh.register("addSelection", addSelection)
    w.eh.register("repSelection", repSelection)
    """
    w.show()
    if QTVERSION < '4.0.0':
        a.exec_loop()
    else:
        a.exec_()

if __name__=="__main__":
    test()
