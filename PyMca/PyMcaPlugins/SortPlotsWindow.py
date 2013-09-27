import numpy
from PyMca import PyMcaQt as qt
from PyMca.PyMca_Icons import IconDict

from PyMca import CloseEventNotifyingWidget
from PyMca import ScanWindow as sw

DEBUG = True
#class SortPlotsScanWindow(qt.QWidget):
class SortPlotsScanWindow(sw.ScanWindow):
    
    plotModifiedSignal = qt.pyqtSignal()
    
    def __init__(self,  origin,
                         name='',  
                         selection=[],  
                         parent=None, 
                         remove_spikes = False, 
                         noise_filter = False,
                         plot_xas = False, 
                         normalize = False):
        sw.ScanWindow.__init__(self, 
                                             parent, 
                                             name='Scan Window'+' '+name, 
                                             specfit=None)
        self.plotWindow = origin
        self.title = name
        
        self.scanWindowInfoWidget.hide()
        
        buttonAdd = qt.QPushButton('Add',  self)
        buttonReplace = qt.QPushButton('Replace',  self)
        buttonAddAll = qt.QPushButton('Add all',  self)
        buttonReplaceAll = qt.QPushButton('Replace all',  self)
        
        buttonLayout = qt.QHBoxLayout(None)
        buttonLayout.setContentsMargins(0, 0, 0, 0)
        buttonLayout.setSpacing(5)
        buttonLayout.addWidget(qt.HorizontalSpacer(self))
        buttonLayout.addWidget(buttonAdd)
        buttonLayout.addWidget(buttonAddAll)
        buttonLayout.addWidget(buttonReplace)
        buttonLayout.addWidget(buttonReplaceAll)
        
        self.mainLayout.addLayout(buttonLayout)
        
        buttonAdd.clicked.connect(self.add)
        buttonReplace.clicked.connect(self.replace)
        buttonAddAll.clicked.connect(self.addAll)
        buttonReplaceAll.clicked.connect(self.replaceAll)
        
        # Computation starts HERE
        pwActiveCurve = self.plotWindow.graph.getActiveCurve(justlegend=True)
        pwX1min,  pwX1max = self.plotWindow.graph.getX1AxisLimits()
        for legend in selection:
            self.plotWindow.setActiveCurve( legend )
            tmp = self.plotWindow.getActiveCurve()
            if tmp:
                (xVal,  yVal,  leg,  info) = tmp
            else:
                continue
            if (xVal[0] < pwX1min) or (xVal[-1] > pwX1max):
                # Clip image to zoomed in Values
                mask = numpy.nonzero((xVal > pwX1min) &
                                                     (xVal < pwX1max))[0]
                xVal = numpy.take(xVal,  mask)
                yVal = numpy.take(yVal,  mask)
            if remove_spikes:
                self.spikeRemoval(yVal)
            if 'selectionlegend' in info:
                newLegend = info['selectionlegend']
            elif 'operation' in info:
                newLegend = (str(operation) + ' ' + self.title)
            else:
                newLegend = leg
            self.addCurve(xVal,  yVal,  
                                      newLegend, 
                                      info)
        if pwActiveCurve:
            self.plotWindow.setActiveCurve(pwActiveCurve)
        
    def spikeRemoval(self,  inp,
                     threshold=0.001,
                     length=5):
        '''
        inp         -> narray type
        threshold   -> float
        length      -> odd integer
        '''
        n = length // 2
        N = len(inp)-1
        mask    = numpy.arange(-n, n+1, 
                               dtype=int)
        maskN   = numpy.array([N]*length)
        # Set up filter array
        for i in range(N+1):
            cut = mask + i
            chk0 = (cut<0)
            chk1 = (cut>N)
            if chk0.any():
                cut *= -chk0
            if chk1.any():
                cut *= -chk1
                cut += (maskN * chk1)
            med = numpy.median(inp.take(cut))
            # Spike removal
            if abs(inp[i]-med)>threshold:
                inp[i] = med

    def normalize(self):
        pass

    def noiseFilter(self):
        pass
    
    def copyCurves(self):
        pass

    def add(self):
        if DEBUG:
            print 'add():'
        activeCurve = self.getActiveCurve()
        if activeCurve is None:
            return
        (xVal,  yVal,  legend,  info) = activeCurve
        if 'selectionlegend' in info:
            newLegend = info['selectionlegend']
        elif 'operation' in info:
            newLegend = (str(operation) + ' ' + self.title)
        else:
            newLegend = (legend + ' ' + self.title)
        self.plotWindow.addCurve(xVal,  yVal,  
                                                 newLegend,
                                                 info)
        self.plotModifiedSignal.emit()

    def addAll(self):
        if DEBUG:
            print 'addAll():'
        for (xVal,  yVal,  legend,  info) in self.getAllCurves():
            if 'selectionlegend' in info:
                newLegend = info['selectionlegend']
            elif 'operation' in info:
                newLegend = (str(operation) + ' ' + self.title)
            else:
                newLegend = (legend + ' ' + self.title)
            self.plotWindow.addCurve(xVal,  yVal, 
                                                     newLegend, 
                                                     info)
        self.plotModifiedSignal.emit()

    def replace(self):
        if DEBUG:
            print 'replace():'
        activeCurve = self.getActiveCurve()
        if activeCurve is None:
            return
        (xVal,  yVal,  legend,  info) = activeCurve
        if 'selectionlegend' in info:
            newLegend = info['selectionlegend']
        elif 'operation' in info:
            newLegend = (str(operation) + ' ' + self.title)
        else:
            newLegend = (legend + self.title)
        self.plotWindow.addCurve(xVal,  yVal, 
                                                newLegend,  
                                                info,  
                                                replace=True)
        self.plotModifiedSignal.emit()

    def replaceAll(self):
        if DEBUG:
            print 'replaceAll()'
        allCurves = self.getAllCurves()
        for (i, (xVal,  yVal,  legend,  info)) in enumerate(allCurves):
            if 'selectionlegend' in info:
                newLegend = info['selectionlegend']
            elif 'operation' in info:
                newLegend = (str(operation) + ' ' + self.title)
            else:
                newLegend = (legend + ' ' + self.title)
            if i == 0:
                self.plotWindow.addCurve(xVal,  yVal,  
                                                         newLegend, 
                                                         info,  
                                                         replace=True)
            else:
                self.plotWindow.addCurve(xVal,  yVal, 
                                                         newLegend, 
                                                         info)
        self.plotModifiedSignal.emit()

    def closeEvent(self,  event):
        self.close()

class SortPlotsMenu(qt.QMenu):
    def __init__(self,  parent,  functionList):
        '''
        List functions has to have the form (functionName, function)
        
        Default is ('', function)
        '''
        qt.QMenu.__init__(self,  parent)
        for (name, function) in functionList:
            if name != '':
                fName = name
            else:
                fName = function.func_name
            act = qt.QAction( fName,  self )
            act.triggered.connect( function )
            self.addAction( act )

class SortPlotsTreeWidget(qt.QTreeWidget):
    def __init__(self,  parent):
        qt.QTreeWidget.__init__(self,  parent)

    def contextMenuEvent(self,  event):
        if event.reason() == event.Mouse:
            pos = event.globalPos()
            item = self.itemAt(event.pos())
        else:
            pos = None
            sel = self.selectedItems()
            if sel:
                item = sel[0]
            else:
                item = self.currentItem()
                if item is None:
                    self.invisibleRootItem().child(0)
            if item is not None:
                itemrect = self.visualItemRect(item)
                portrect = self.viewport().rect()
                itemrect.setLeft(portrect.left())
                itemrect.setWidth(portrect.width())
                pos = self.mapToGlobal(itemrect.bottomLeft())
        if pos is not None:
            menu = qt.QMenu('Perform..',  self)
            menu.addActions(self.parentWidget().actionList)
            menu.popup(pos)
        event.accept()
        
    def invertSelection(self):
        root = self.invisibleRootItem()
        for i in range(root.childCount()):
            if root.child(i).isSelected():
                root.child(i).setSelected( False )
            else:
                root.child(i).setSelected( True )

    def selectedItems(self):
        sel = super(SortPlotsTreeWidget,  self).selectedItems()
        ret = []
        for item in sel:
            # Convert from QTreeWidgetItem to str
            # Only use selected legend
            ret += [ str(item.text(0)) ]
        if DEBUG:
            print 'selectedItems: %d Item(s) selected'%len(sel)
        return ret

    def build(self,  items,  headerLabels):
        '''
        (Re-) Builds the tree display
        
        headerLabels must be of type QStringList
        items must be of type [QStringList]
        '''
        sel = self.selectedItems() # Remember selection
        self.clear()
        self.setHeaderLabels(headerLabels)
        for item in items:
            treeItem = qt.QTreeWidgetItem(self,  item)
            if treeItem.text(0) in sel:
                treeItem.setSelected(True)

class SortPlotsInstructionWidget(qt.QWidget):
    def __init__(self,  parent):
        qt.QWidget.__init__(self,  parent)
        buttonPerfom = qt.QPushButton("Perform action",  self)
        buttonProcessA =  qt.QPushButton("Process as A",  self)
        buttonProcessB =  qt.QPushButton("Process as B",  self)
        self.actionCBox = qt.QComboBox(self)
        
        buttonLayout = qt.QHBoxLayout(None)
        buttonLayout.addWidget(buttonProcessA)
        buttonLayout.addWidget(buttonProcessB)
        buttonLayout.addWidget(qt.HorizontalSpacer(self))
        buttonLayout.addWidget(self.actionCBox) 
        buttonLayout.addWidget(buttonPerfom)
        
        buttonPerfom.clicked.connect( self.performAction )
        buttonProcessA.clicked.connect( parent.processAsA )
        buttonProcessB.clicked.connect( parent.processAsB )
        
        buttonLayout.setContentsMargins(1, 1, 1, 1)
        buttonLayout.setSpacing(2)
        
        self.setLayout(buttonLayout)

    def performAction(self):
        index = self.actionCBox.currentIndex()
        self.parentWidget().actionList[index].trigger()

class SortPlotsWidget(CloseEventNotifyingWidget.CloseEventNotifyingWidget):
    
    actionListModifiedSignal = qt.pyqtSignal()
    
    def __init__(self,  parent,
                        legends,
                        motorValues,
                        plotWindow = None,
                        actions = [],
                        nSelectors = 2, 
                        instructions = False):
        """
        legends            List contains Plotnames
        motorValues     List contains names and values of the motors
        """
        CloseEventNotifyingWidget.CloseEventNotifyingWidget.__init__(self,  parent)
        self.setWindowTitle("Sort Plots Window")
        
        self.legendList = legends
        self.motorsList = motorValues
        self.motorNamesList = [''] + self.getAllMotorNames()
        self.motorNamesList.sort()
        self.numCurves = len(legends)
        self.actionList = actions
        self.cBoxList = []
        self.ScanWindowA = None
        self.ScanWindowB = None
        self.closeWidget = CloseEventNotifyingWidget.CloseEventNotifyingWidget()
        self.plotWindow = plotWindow

        updatePixmap = qt.QPixmap(IconDict["reload"])
        buttonUpdate = qt.QPushButton(qt.QIcon(updatePixmap), 
                                                         '', 
                                                         self)
        
        cBoxLabel = qt.QLabel( qt.QString('Select motor:'),  self)
        for i in range(nSelectors):
            cBox = qt.QComboBox(self)
            cBox.addItems( self.motorNamesList )
            cBox.activated['QString'].connect(self.updateTree)
            self.cBoxList += [cBox]
        
        self.list = SortPlotsTreeWidget(self)
        labels = ['Legend'] + nSelectors*['']
        ncols  = len(labels)
        self.list.setColumnCount(ncols)
        self.list.setHeaderLabels(labels)
        self.list.setSortingEnabled(True)
        self.list.setSelectionMode(qt.QAbstractItemView.ExtendedSelection)
        
        mainLayout = qt.QGridLayout(self)
        cBoxLayout = qt.QHBoxLayout(None)
        mainLayout.setContentsMargins(1, 1, 1, 1)
        mainLayout.setSpacing(2)
        mainLayout.addLayout(cBoxLayout,  0, 0)
#        mainLayout.addLayout(buttonLayout,  2, 0)
        self.setLayout(mainLayout)

        cBoxLayout.addWidget(cBoxLabel)
        for cBox in self.cBoxList:
            cBoxLayout.addWidget(cBox)
        cBoxLayout.addWidget(qt.HorizontalSpacer(self))
        cBoxLayout.addWidget(buttonUpdate)
        mainLayout.addWidget(self.list,  1,  0)
        if instructions:
            self.instWidget = SortPlotsInstructionWidget(self)
            mainLayout.addWidget( self.instWidget )
        else:
            self.instWidget = None
        self.resize(500,  300)
        
        buttonUpdate.clicked.connect( self.updatePlots )

        self.updateTree()
        self.updateActionList(
              [('Invert selection',  self.list.invertSelection), 
                ('Process as A',  self.processAsA), 
                ('Process as B',  self.processAsB), 
                ('Remove curve(s)',  self.removeCurve_)])

    def getAllMotorNames(self):
        nameSet = set()
        for dic in self.motorsList:
            for key in dic.keys(): nameSet.add( key )
        return list( nameSet )

    def removeCurve_(self):
        sel = self.list.selectedItems()
        self.plotWindow.removeCurves(sel)
        self.updatePlots()

    def updatePlots(self,  newLegends = None,  
                                            newMotorValues = None):
        self._setLists()
        self.motorNamesList = [''] + self.getAllMotorNames()
        self.motorNamesList.sort()
        for cBox in self.cBoxList:
            index = cBox.currentIndex()
            cBox.clear()
            cBox.addItems( self.motorNamesList )
            cBox.setCurrentIndex(index)
        self.updateTree()

    def updateTree(self):
        mList = [ str( cBox.currentText() ) for cBox in self.cBoxList ]
        labels = ['Legend'] + mList
        items = []
        for i in range(len(self.legendList)):
            tmp = qt.QStringList()
            legend = self.legendList[i]
            values = self.motorsList[i]
            tmp.append( legend )
            for m in mList:
                if m == '':
                    tmp.append( '' )
                else:
                    tmp.append( str(values.get( m, '---' )) )
            items.append(tmp)
        self.list.build(items,  labels)
    
    def updateActionList(self,  functionList):
        if self.instWidget:
            self.instWidget.actionCBox.clear()
        for (name, function) in functionList:
            fName = name if name is not '' else function.func_name
            act = qt.QAction( fName,  self )
            act.triggered.connect( function )
            self.actionList += [ act ]
            if self.instWidget:
                self.instWidget.actionCBox.addItem( fName )  
    
    def processAsA(self):
        if self.ScanWindowA:
            if self.ScanWindowA.isVisible():
                self.ScanWindowA.close()
        self.ScanWindowA = self.spawnNewScanWindow('A')
    
    def processAsB(self):
        if self.ScanWindowB:
            if self.ScanWindowB.isVisible():
                self.ScanWindowB.close()
        self.ScanWindowB = self.spawnNewScanWindow('B')

    def spawnNewScanWindow(self,  name=''):
        sel = self.list.selectedItems()
        swin  = SortPlotsScanWindow(self.plotWindow, 
                                    name, 
                                    sel,
                                    None, 
                                    remove_spikes = True)
        swin.plotModifiedSignal.connect(self.updatePlots)
        swin.show()
        swin.raise_()
        self.notifyCloseEventToWidget(swin)
        return swin

    def _convertInfoDictionary(self,  infosList):
        ret = []
        for info in infosList :
            motorNames = info.get('MotorNames',  None)
            if motorNames is not None:
                if type(motorNames) == str:
                    namesList = motorNames.split()
                elif type(motorNames) == list:
                    namesList = motorNames
                else:
                    namesList = []
            else:
                namesList = []
            motorValues = info.get('MotorValues',  None)
            if motorNames is not None:
                if type(motorValues) == str:
                    valuesList = motorValues.split()
                elif type(motorValues) == list:
                    valuesList = motorValues
                else:
                    valuesList = []
            else:
                valuesList = []
            if len(namesList) == len(valuesList):
                ret.append( dict( zip( namesList,  valuesList ) ) )
            else:
                print("Number of motors and values does not match!")
        return ret

    def _setLists(self):
        curves = self.plotWindow.getAllCurves()
        nCurves = len(curves)
        if DEBUG:
            print ("Received %d curve(s).." % nCurves)
        self.legendList = [leg for (xvals, yvals,  leg,  info) in curves] 
        infoList = [info for (xvals, yvals,  leg,  info) in curves] 
        self.motorsList = self._convertInfoDictionary( infoList )

    def closeEvent(self,  event):
        for swin in [self.ScanWindowA,  self.ScanWindowB]:
            if swin:
                swin.close()
        self.close()
 
def main():
    import sys,  numpy
    app = qt.QApplication(sys.argv)
    swin = sw.ScanWindow()
    legends = ['Curve0', 'Curve1', 'Curve2']
    motors = [{'Motor12': 1, 'Motor11': 8.692713996985609, 'Motor10': 21.98364185388587, 'Motor 8': 0.19806882661182112, 'Motor 9': 0.4844754557916431, 'Motor 4': 0.3502522172639875, 'Motor 5': 0.6639252709334457, 'Motor 6': 0.8130332644206067, 'Motor 7': 0.22114941021809853, 'Motor 0': 0.5931882588655031, 'Motor 1': 0.6780103928805297, 'Motor 2': 0.26738924783290086, 'Motor 3': 0.6778906178576761}, {'Motor18': 0.4707468826876532, 'Motor17': 0.6958160702991127, 'Motor16': 0.8257808117546283, 'Motor15': 0.2587637453100148, 'Motor14': 0.7392644674355958, 'Motor13': 0.09084289261899736, 'Motor12': 2, 'Motor11': 0.21344565983311958, 'Motor10': 0.823400550314221, 'Motor 8': 0.020278096856981342, 'Motor 9': 0.5687440213219551, 'Motor 4': 0.8537811553701731, 'Motor 5': 0.6967303868907243, 'Motor 6': 0.2691963139564302, 'Motor 7': 0.7932933343951395, 'Motor 0': 0.7692165677566825, 'Motor 1': 0.9590927095265979, 'Motor 2': 0.010926468369733544, 'Motor 3': 0.5382649725528551}, {'Motor12': 2, 'Motor11': 0.44400576643956124, 'Motor10': 0.613870067851634, 'Motor 8': 0.901968648110583, 'Motor 9': 0.3197687710845185, 'Motor 4': 0.5714322786278168, 'Motor 5': 0.2786758361634877, 'Motor 6': 0.15443677487828655, 'Motor 7': 0.41623199933237764, 'Motor 0': 0.294201017230741, 'Motor 1': 0.813913587747513, 'Motor 2': 0.5775729031053222, 'Motor 3': 0.8690451825680668}, {'Motor13': 0.6491598094029021, 'Motor12': 10, 'Motor11': 0.006312468992195397, 'Motor10': 0.06727805971206435, 'Motor 8': 0.0929878987747117, 'Motor 9': 0.014325738753558803, 'Motor 4': 0.8185362197656616, 'Motor 5': 0.6643614796103005, 'Motor 6': 0.6479279384366304, 'Motor 7': 0.3485172683358245, 'Motor 0': 0.9858738343685299, 'Motor 1': 0.9330130170323839, 'Motor 2': 0.7550180320112966, 'Motor 3': 0.8814284215685484}, {'Motor19': 0.39846564175862953, 'Motor18': 0.2745751180457152, 'Motor17': 0.42793840508599434, 'Motor16': 0.5335910248322966, 'Motor15': 0.14010423968992758, 'Motor14': 0.27948624022431734, 'Motor13': 0.1737756266389101, 'Motor12': 0.6425110521350722, 'Motor11': 0.9040646490476784, 'Motor10': 0.22997142790156133, 'Motor 8': 0.3520106476992403, 'Motor 9': 0.37023110928070235, 'Motor 4': 0.8110924828319052, 'Motor 5': 0.854155188450653, 'Motor 6': 0.12438157550841666, 'Motor 7': 0.3303770832430888, 'Motor 0': 0.4583273673870403, 'Motor 1': 0.40863603059350373, 'Motor 2': 0.7396799985670546, 'Motor 3': 0.5532134465740317, 'Motor22': 0.7154261407207922, 'Motor20': 0.6735594219326284, 'Motor21': 0.24068704947080943}, {'Motor18': 0.7501922139242619, 'Motor17': 0.067572631661458, 'Motor16': 0.23941863624378346, 'Motor15': 0.543195970137226, 'Motor14': 0.5045110454536483, 'Motor13': 0.47129338234441986, 'Motor12': 0.7039345533241258, 'Motor11': 0.5496976809598649, 'Motor10': 0.028685484457880994, 'Motor 8': 0.3736138811685542, 'Motor 9': 0.6200990287805606, 'Motor 4': 0.30138047598948403, 'Motor 5': 0.15683187764664286, 'Motor 6': 0.061169736595949264, 'Motor 7': 0.35931932492621954, 'Motor 0': 0.7241839150429988, 'Motor 1': 0.7985803970529565, 'Motor 2': 0.5239059568843569, 'Motor 3': 0.7404964999807312}, {'Motor10': 0.90828582481094, 'Motor 8': 0.8424405354748069, 'Motor 9': 0.021278797555318363, 'Motor 4': 0.8593234401902958, 'Motor 5': 0.2638651881043157, 'Motor 6': 0.281687767263718, 'Motor 7': 0.48283570902507555, 'Motor 0': 0.659487116102895, 'Motor 1': 24.591253182578376, 'Motor 2': 3.032078904732739, 'Motor 3': 0.17860013910027928}, {'Motor 8': 0.7246181445974952, 'Motor 9': 0.5375876404160089, 'Motor 4': 0.7608877399780997, 'Motor 5': 0.6164359666836775, 'Motor 6': 0.3910546574315933, 'Motor 7': 0.5287834048239588, 'Motor 0': 0.9700467881758079, 'Motor 1': 0.9064128957850547, 'Motor 2': 0.4434306640093745, 'Motor 3': 0.2783396189782661}, {'Motor19': 0.4741833534896892, 'Motor18': 0.1884371839846597, 'Motor17': 0.660882814263354, 'Motor16': 0.25871486157318313, 'Motor15': 0.6181192138005907, 'Motor14': 0.11534451504645371, 'Motor13': 0.3356756251510249, 'Motor12': 0.8578128852052718, 'Motor11': 0.002943123668270098, 'Motor10': 0.08980970319869397, 'Motor 8': 0.40586648583549123, 'Motor 9': 0.7700310455423328, 'Motor 4': 0.8389920867382025, 'Motor 5': 0.2560110245056251, 'Motor 6': 0.671297941874289, 'Motor 7': 0.7041220063735543, 'Motor 0': 0.4865107750866541, 'Motor 1': 0.8623573559114868, 'Motor 2': 0.8378911209243649, 'Motor 3': 0.056056301247044416, 'Motor24': 0.8535082807686701, 'Motor22': 0.4362354327544248, 'Motor23': 0.17386904782647783, 'Motor20': 0.11001296204329247, 'Motor21': 0.5653716280128318}, {'Motor13': 0.5900826517637087, 'Motor12': 0.2876746207456713, 'Motor11': 0.1829075413610104, 'Motor10': 0.9677552520998641, 'Motor 8': 0.47506344789108046, 'Motor 9': 0.32097198197020305, 'Motor 4': 0.5708449042766175, 'Motor 5': 0.06093583375842648, 'Motor 6': 0.10172375432338043, 'Motor 7': 0.989917381621416, 'Motor 0': 0.8047039621208083, 'Motor 1': 0.9477209087673744, 'Motor 2': 0.46582818765280054, 'Motor 3': 0.0511893987634543},  {'Motor12': 0.6504890336783156, 'Motor11': 0.44400576643956124, 'Motor10': 0.613870067851634, 'Motor 8': 0.901968648110583, 'Motor 9': 0.3197687710845185, 'Motor 4': 0.5714322786278168, 'Motor 5': 0.2786758361634877, 'Motor 6': 0.15443677487828655, 'Motor 7': 0.41623199933237764, 'Motor 0': 0.294201017230741, 'Motor 1': 0.813913587747513, 'Motor 2': 0.5775729031053222, 'Motor 3': 0.8690451825680668}, {'Motor13': 0.6491598094029021, 'Motor12': 0.2975843286841311, 'Motor11': 0.006312468992195397, 'Motor10': 0.06727805971206435, 'Motor 8': 0.0929878987747117, 'Motor 9': 0.014325738753558803, 'Motor 4': 0.8185362197656616, 'Motor 5': 0.6643614796103005, 'Motor 6': 0.6479279384366304, 'Motor 7': 0.3485172683358245, 'Motor 0': 0.9858738343685299, 'Motor 1': 0.9330130170323839, 'Motor 2': 0.7550180320112966, 'Motor 3': 0.8814284215685484},  {'Motor12': 0.6504890336783156, 'Motor11': 0.44400576643956124, 'Motor10': 0.613870067851634, 'Motor 8': 0.901968648110583, 'Motor 9': 0.3197687710845185, 'Motor 4': 0.5714322786278168, 'Motor 5': 0.2786758361634877, 'Motor 6': 0.15443677487828655, 'Motor 7': 0.41623199933237764, 'Motor 0': 0.294201017230741, 'Motor 1': 0.813913587747513, 'Motor 2': 0.5775729031053222, 'Motor 3': 0.8690451825680668}, {'Motor13': 0.6491598094029021, 'Motor12': 0.2975843286841311, 'Motor11': 0.006312468992195397, 'Motor10': 0.06727805971206435, 'Motor 8': 0.0929878987747117, 'Motor 9': 0.014325738753558803, 'Motor 4': 0.8185362197656616, 'Motor 5': 0.6643614796103005, 'Motor 6': 0.6479279384366304, 'Motor 7': 0.3485172683358245, 'Motor 0': 0.9858738343685299, 'Motor 1': 0.9330130170323839, 'Motor 2': 0.7550180320112966, 'Motor 3': 0.8814284215685484}]
    x = numpy.arange(1000.)
    y0 =  10 * x + 10000. * numpy.exp(-0.5*(x-500)*(x-500)/400)
    y1 =  10 * x + 10000. * numpy.exp(-0.5*(x-600)*(x-600)/400)
    y2 =  10 * x + 10000. * numpy.exp(-0.5*(x-400)*(x-400)/400)
    swin.addCurve(x, y2, legend="Curve2", replot=False, replace=False)
    swin.addCurve(x, y0, legend="Curve0", replot=False, replace=False)
    swin.addCurve(x, y1, legend="Curve1", replot=False, replace=False)

    w = SortPlotsWidget(None, swin.getAllCurves(just_legend=True),  motors,  swin,  instructions = True)
    w.show()
    app.exec_()
    
if __name__ == '__main__':
    main()
