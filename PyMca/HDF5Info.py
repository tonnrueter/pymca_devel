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
try:
    from PyMca import PyMcaQt as qt
    safe_str = qt.safe_str
except ImportError:
    # for people using this widget without PyMca installed
    import PyQt4.Qt as qt
    safe_str = str
import copy
import posixpath

class HDFInfoCustomEvent(qt.QEvent):
    def __init__(self, ddict):
        if ddict is None:
            ddict = {}
        self.dict = ddict
        qt.QEvent.__init__(self, qt.QEvent.User)

class VerticalSpacer(qt.QWidget):
    def __init__(self, *args):
        qt.QWidget.__init__(self, *args)
        self.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed,
                                          qt.QSizePolicy.Expanding))

class SimpleInfoGroupBox(qt.QGroupBox):
    def __init__(self, parent, title=None, keys=None):
        qt.QGroupBox.__init__(self, parent)
        self.mainLayout = qt.QGridLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(2)
        if title is not None:
            self.setTitle(title)
        if keys is None:
           keys = [] 
        self.keyList = keys
        self.keyDict = {}
        if len(self.keyList):
            self._build()

    def _build(self):
        i = 0
        for key in self.keyList:
            label = qt.QLabel(self)
            label.setText(key)
            line = qt.QLineEdit(self)
            line.setReadOnly(True)
            self.mainLayout.addWidget(label, i, 0)
            self.mainLayout.addWidget(line, i, 1)
            self.keyDict[key] = (label, line)
            i += 1

    def setInfoDict(self, ddict):
        if not len(self.keyList):
            self.keyList = ddict.keys()
            self._build()            
        self._fillInfo(ddict)

    def _fillInfo(self, ddict0):
        ddict = self._getMappedDict(ddict0)
        actualKeys = list(ddict.keys())
        dictKeys = []
        for key in actualKeys:
            l = key.lower()
            if l not in dictKeys:
                dictKeys.append(l)
        for key in self.keyList:
            l = key.lower()
            if l in dictKeys:
                self._fillKey(key,
                              ddict[actualKeys[dictKeys.index(l)]])

    def _getMappedDict(self, ddict):
        #Default implementation returns a copy of the input dictionnary
        return copy.deepcopy(ddict)

    def _fillKey(self, key, value):
        #This can be overwritten
        if type(value) == type(""):
            self.keyDict[key][1].setText(value)
        else:
            self.keyDict[key][1].setText(safe_str(value))

class NameGroupBox(SimpleInfoGroupBox):
    def __init__(self, parent, title=None, keys=[]):
        SimpleInfoGroupBox.__init__(self, parent, title=title, keys=["Name", "Path", "Type"])

    def setInfoDict(self, ddict):
        key = "Value"
        if key in ddict.keys():
            if key not in self.keyList:
                self.keyList.append(key)
                label = qt.QLabel(self)
                label.setText(key)
                line = qt.QLineEdit(self)
                line.setReadOnly(True)
                i = self.keyList.index(key)
                self.mainLayout.addWidget(label, i, 0)
                self.mainLayout.addWidget(line, i, 1)
                self.keyDict[key] = (label, line)
        if 'Path' in ddict:
            if ddict['Path'] == "/":
                if 'Name' in ddict:
                    self.keyDict['Name'][0].setText("File")
        SimpleInfoGroupBox.setInfoDict(self, ddict)

class DimensionGroupBox(SimpleInfoGroupBox):
    def __init__(self, parent, title=None, keys=None):
        keys = ["No. of Dimension(s)",
                "Dimension Size(s)",
                "Data Type"]
        SimpleInfoGroupBox.__init__(self, parent, title=title, keys=keys)

    def _getMappedDict(self, ddict):
        return copy.deepcopy(ddict)

class MembersGroupBox(qt.QGroupBox):
    def __init__(self, parent):
        qt.QGroupBox.__init__(self, parent)
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(2)
        self.setTitle("Group Members")
        self.label = qt.QLabel(self)
        self.label.setText("Number of members: 0")
        self.table = qt.QTableWidget(self)
        #labels = ["Name", "Type", "Shape", "Value"]
        labels = ["Name", "Value", "Type"]
        nlabels = len(labels)
        self.table.setColumnCount(nlabels)
        rheight = self.table.horizontalHeader().sizeHint().height()
        self.table.setMinimumHeight(12*rheight)
        self.table.setMaximumHeight(20*rheight)
        for i in range(nlabels):
            item = self.table.horizontalHeaderItem(i)
            if item is None:
                item = qt.QTableWidgetItem(labels[i],
                                           qt.QTableWidgetItem.Type)
            item.setText(labels[i])
            self.table.setHorizontalHeaderItem(i, item)
        self._tableLabels = labels
        self.mainLayout.addWidget(self.label)
        self.mainLayout.addWidget(self.table)

    def setInfoDict(self, ddict):
        keylist = ddict.keys()
        if "members" not in keylist:
            self.hide()
            return
        keylist = ddict['members']
        self.label.setText("Number of members: %d" % len(keylist))
        nrows = len(keylist)
        if not nrows:
            self.table.setRowCount(nrows)
            self.hide()
            return
        self.table.setRowCount(nrows)
        if ddict['Path'] != '/':
            #this could distroy ordering ...
            keylist.sort()
        row = 0
        for key in keylist:
            item = self.table.item(row, 0)
            if item is None:
                item = qt.QTableWidgetItem(key, qt.QTableWidgetItem.Type)
                item.setFlags(qt.Qt.ItemIsSelectable|
                                  qt.Qt.ItemIsEnabled)
                self.table.setItem(row, 0, item)
            else:
                item.setText(key)
            for label in self._tableLabels[1:]:
                if not label in ddict[key]:
                    continue
                col = self._tableLabels.index(label)
                info = ddict[key][label]
                item = self.table.item(row, col)    
                if item is None:
                    item = qt.QTableWidgetItem(info, qt.QTableWidgetItem.Type)
                    item.setFlags(qt.Qt.ItemIsSelectable|
                                      qt.Qt.ItemIsEnabled)
                    self.table.setItem(row, col, item)
                else:
                    item.setText(info)                    
            row += 1
        for i in range(self.table.columnCount()):
            self.table.resizeColumnToContents(i)
        
class HDF5GeneralInfoWidget(qt.QWidget):
    def __init__(self, parent=None, ddict=None):
        qt.QWidget.__init__(self, parent)
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(2)
        self.nameWidget = NameGroupBox(self)
        self.membersWidget = MembersGroupBox(self)
        self.dimensionWidget = DimensionGroupBox(self)
        self.mainLayout.addWidget(self.nameWidget)
        self.mainLayout.addWidget(self.membersWidget)
        self.mainLayout.addWidget(VerticalSpacer(self))
        self.mainLayout.addWidget(self.dimensionWidget)
        self._notifyCloseEventToWidget = None
        
        if ddict is not None:
            self.setInfoDict(ddict)

    def setInfoDict(self, ddict):
        if 'general' in ddict:
            self._setInfoDict(ddict['general'])
        else:
            self._setInfoDict(ddict)

    def _setInfoDict(self, ddict):
        self.nameWidget.setInfoDict(ddict)
        self.membersWidget.setInfoDict(ddict)
        self.dimensionWidget.setInfoDict(ddict)
        if 'members' in ddict:
            if len(ddict['members']):
                #it is a datagroup
                self.dimensionWidget.hide()  
        self.dimensionWidget.hide()


class HDF5AttributesInfoWidget(qt.QWidget):
    def __init__(self, parent):
        qt.QGroupBox.__init__(self, parent)
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setMargin(0)
        self.mainLayout.setSpacing(2)
        self.label = qt.QLabel(self)
        self.label.setText("Number of members: 0")
        self.table = qt.QTableWidget(self)
        labels = ["Name", "Value", "Type", "Size"]
        self.table.setColumnCount(len(labels))
        for i in range(len(labels)):
            item = self.table.horizontalHeaderItem(i)
            if item is None:
                item = qt.QTableWidgetItem(labels[i],
                                           qt.QTableWidgetItem.Type)
            item.setText(labels[i])
            self.table.setHorizontalHeaderItem(i, item)
        self._tableLabels = labels
        self.mainLayout.addWidget(self.label)
        self.mainLayout.addWidget(self.table)

    def setInfoDict(self, ddict):
        if 'attributes' in ddict:
            self._setInfoDict(ddict['attributes'])
        else:
            self._setInfoDict(ddict)


    def _setInfoDict(self, ddict):
        keylist = ddict['names']
        self.label.setText("Number of attributes: %d" % len(keylist))
        nrows = len(keylist)
        if not nrows:
            self.table.setRowCount(nrows)
            self.hide()
            return
        self.table.setRowCount(nrows)        
        keylist.sort()
        row = 0
        for key in keylist:
            for label in self._tableLabels:
                if not label in ddict[key]:
                    continue
                else:
                    text = ddict[key][label]
                col = self._tableLabels.index(label)
                item = self.table.item(row, col)
                if item is None:
                    item = qt.QTableWidgetItem(text, qt.QTableWidgetItem.Type)
                    item.setFlags(qt.Qt.ItemIsSelectable|
                                      qt.Qt.ItemIsEnabled)
                    self.table.setItem(row, col, item)
                else:
                    item.setText(text)
            row += 1

        #for i in range(self.table.columnCount()):
        self.table.resizeColumnToContents(1)
        self.table.resizeColumnToContents(3)
        
class HDF5InfoWidget(qt.QTabWidget):
    def __init__(self, parent=None, info=None):
        qt.QTabWidget.__init__(self, parent)
        self._notifyCloseEventToWidget = []
        self._build()
        if info is not None:
            self.setInfoDict(info)

    def sizeHint(self):
        return qt.QSize(2 * qt.QTabWidget.sizeHint(self).width(),
                        qt.QTabWidget.sizeHint(self).height())
                        
    def _build(self):
        self.generalInfoWidget = HDF5GeneralInfoWidget(self)
        self.attributesInfoWidget = HDF5AttributesInfoWidget(self)
        self.addTab(self.generalInfoWidget, 'General')
        self.addTab(self.attributesInfoWidget, 'Attributes')

    def setInfoDict(self, ddict):
        self.generalInfoWidget.setInfoDict(ddict)
        self.attributesInfoWidget.setInfoDict(ddict)

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
                try:
                    qt.QApplication.postEvent(widget, newEvent)
                except:
                    print("Error notifying close event to widget", widget)
            self._notifyCloseEventToWidget = []
        return qt.QWidget.closeEvent(self, event)

def getInfo(hdf5File, node):
    """
    hdf5File is and HDF5 file-like insance
    node is the posix path to the node
    """
    data = hdf5File[node]
    ddict = {}
    ddict['general'] = {}
    ddict['attributes'] = {}

    externalFile = False
    if hasattr(hdf5File, "file"):
        if hasattr(hdf5File.file, "filename"):
            if data.file.filename != hdf5File.file.filename:
                externalFile = True

    if externalFile:
        path = node
        if path != "/":
            name = posixpath.basename(node)
        else:
            name = hdf5File.file.filename
        if data.name != node:
            path += " external link to %s" % data.name
            path += " in %s" % safe_str(data.file.filename)
            name += " external link to %s" % data.name
            name += " in %s" % safe_str(data.file.filename)
        else:
            name = data.name
        ddict['general']['Path'] = path
        ddict['general']['Name'] = name
    else:
        ddict['general']['Path'] = data.name
        if ddict['general']['Path'] != "/":
            ddict['general']['Name'] = posixpath.basename(data.name)
        else:
            ddict['general']['Name'] = data.file.filename
    ddict['general']['Type'] = safe_str(data)
    if hasattr(data, 'dtype'):
        if ("%s" % data.dtype).startswith("|S") or\
           ("%s" % data.dtype).startswith("|O"):
            if hasattr(data, 'shape'):
                shape = data.shape
                if not len(shape):
                    ddict['general']['Value'] = "%s" % data.value
                elif shape[0] == 1:
                    ddict['general']['Value'] = "%s" % data.value[0]
                else:
                    print("Warning: Node %s not fully understood" % node)
                    ddict['general']['Value'] = "%s" % data.value
        elif hasattr(data, 'shape'):
            shape = data.shape
            if len(shape) == 1:
                if shape[0] == 1:
                    ddict['general']['Value'] = "%s" % data.value[0]
            elif len(shape) == 0:
                ddict['general']['Value'] = "%s" % data.value
    if hasattr(data, "keys"):
        ddict['general']['members'] = list(data.keys())
    elif hasattr(data, "listnames"):
        ddict['general']['members'] = list(data.listnames())
    else:
        ddict['general']['members'] = []
    for member in list(ddict['general']['members']):
        ddict['general'][member] = {}
        ddict['general'][member]['Name'] = safe_str(member)
        if ddict['general']['Path'] == "/":
            ddict['general'][member]['Type'] = safe_str(hdf5File[node+"/"+member])
            continue
        memberObject = hdf5File[node][member]            
        if hasattr(memberObject, 'shape'):
            ddict['general'][member]['Type'] = safe_str(hdf5File[node+"/"+member])
            dtype = memberObject.dtype
            if hasattr(memberObject, 'shape'):
                shape = memberObject.shape
                if ("%s" % dtype).startswith("|S") or\
                   ("%s" % dtype).startswith("|O"):
                    if not len(shape):
                        ddict['general'][member]['Shape'] = ""
                        ddict['general'][member]['Value'] = "%s" % memberObject.value
                    else:
                        ddict['general'][member]['Shape'] = shape[0]
                        ddict['general'][member]['Value'] = "%s" % memberObject.value[0]
                    continue
                if not len(shape):
                    ddict['general'][member]['Shape'] = ""
                    ddict['general'][member]['Value'] = "%s" % memberObject.value
                    continue
                ddict['general'][member]['Shape'] = "%d" % shape[0]
                for i in range(1, len(shape)):
                    ddict['general'][member]['Shape'] += " x %d" % shape[i]
                if len(shape) == 1:
                    if shape[0] == 1:
                        ddict['general'][member]['Value'] = "%s" % memberObject.value[0]
                    else:
                        ddict['general'][member]['Value'] = "%s, ..., %s" % (memberObject.value[0],
                                                                             memberObject.value[-1])
                elif len(shape) == 2:
                    ddict['general'][member]['Value'] = "%s, ..., %s" % (memberObject.value[0],
                                                                         memberObject.value[-1])                                                                         
                else:
                    pass
        else:
            ddict['general'][member]['Type'] = safe_str(hdf5File[node+"/"+member])

    if hasattr(data.attrs, "keys"):
        ddict['attributes']['names'] = data.attrs.keys()
    elif hasattr(data.attrs, "listnames"):
        ddict['attributes']['names'] = data.attrs.listnames()
    else:
        ddict['attributes']['names'] = []
    if sys.version >= '3.0.0':
        ddict['attributes']['names'] = list(ddict['attributes']['names'])
    ddict['attributes']['names'].sort()
    for key in ddict['attributes']['names']:
        ddict['attributes'][key] = {}
        Name = key
        Value = data.attrs[key]
        Type =  safe_str(type(Value))
        if type(Value) == type(""):
            Size = "%d" % len(Value)
        elif type(Value) in [type(1), type(0.0)]:
            Value = safe_str(Value)
            Size = "1"
        else:
            Value = safe_str(Value)
            Size = "Unknown"
        ddict['attributes'][key]['Name']  = Name
        ddict['attributes'][key]['Value'] = Value
        ddict['attributes'][key]['Type']  = Type
        ddict['attributes'][key]['Size']  = Size
    return ddict
        
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage:")
        print("python HDF5Info.py hdf5File node")
        sys.exit(0)
    import h5py
    h=h5py.File(sys.argv[1])
    node = sys.argv[2]
    info = getInfo(h, node)
    app = qt.QApplication([])
    w = HDF5InfoWidget()
    w.setInfoDict(info)
    w.show()
    sys.exit(app.exec_())
