from PyMca import PyMcaQt as qt
from PyMca.PyMca_Icons import IconDict

class SortPlotsMenu(qt.QMenu):
    def __init__(self,  parent,  functionList):
        '''
        List functions has to have the form (functionName, function)
        
        Default is ('', function)
        '''
        qt.QMenu.__init__(self,  parent)
        for (name, function) in functionList:
            fName = name if name is not '' else function.func_name
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
#            menu = SortPlotsMenu(self, [('Invert selection',  self.invertSelection)] )
            menu = qt.QMenu('Perform..',  self)
            menu.addActions(self.parentWidget().actionList)
            menu.popup(pos)
        event.accept()
        
    def invertSelection(self ):
        root = self.invisibleRootItem()
        for i in range(root.childCount()):
            if root.child(i).isSelected():
                root.child(i).setSelected( False )
            else:
                root.child(i).setSelected( True )

class SortPlotsWidget(qt.QWidget):
    def __init__(self,  parent,  legends,  motorValues,  actions = [],  nSelectors = 2):
        """
        legends            List contains Plotnames
        motorValues     List contains names and values of the motors
        """
        qt.QWidget.__init__(self,  parent)
        self.setWindowTitle("Sort Plots Window")
        
        self.legendList = legends
        self.motorsList = motorValues
        self.motorNamesList = [''] + self.getAllMotorNames()
        self.motorNamesList.sort()
        self.numCurves = len(legends)
        self.actionList = actions
        self.cBoxList = []

        # Buttons
        buttonPerfom = qt.QPushButton("Perform action",  self)
        self.buttonUpdate = qt.QPushButton(qt.QIcon(qt.QPixmap(IconDict["reload"])),  '',  self)
        
        # Comboboxes
        cBoxLabel = qt.QLabel( qt.QString('Select motor:'),  self)
        for i in range(nSelectors):
            cBox = qt.QComboBox(self)
            cBox.addItems( self.motorNamesList )
            cBox.activated['QString'].connect(self.updateTree)
            self.cBoxList += [cBox]
        self.actionCBox = qt.QComboBox(self)
        
        # TreeWidget
        self.list = SortPlotsTreeWidget(self)
        labels = ['Legend'] + nSelectors*['']
        ncols  = len(labels)
        self.list.setColumnCount(ncols)
        self.list.setHeaderLabels(labels)
        self.list.setSortingEnabled(True)
        self.list.setSelectionMode(qt.QAbstractItemView.ExtendedSelection)
        
        # Layout
        mainLayout = qt.QGridLayout(self)
        cBoxLayout = qt.QHBoxLayout(None)
        buttonLayout = qt.QHBoxLayout(None)
        mainLayout.setContentsMargins(1, 1, 1, 1)
        mainLayout.setSpacing(2)
        mainLayout.addLayout(cBoxLayout,  0, 0)
        mainLayout.addLayout(buttonLayout,  2, 0)
        self.setLayout(mainLayout)

        # Add widgets to layouts
        cBoxLayout.addWidget(cBoxLabel)
        for cBox in self.cBoxList:
            cBoxLayout.addWidget(cBox)
        cBoxLayout.addWidget(qt.HorizontalSpacer(self))
        cBoxLayout.addWidget(self.buttonUpdate)
        mainLayout.addWidget(self.list,  1,  0)
        buttonLayout.addWidget(self.actionCBox)
        buttonLayout.addWidget(qt.HorizontalSpacer(self))
        buttonLayout.addWidget(buttonPerfom)
        self.resize(500,  300)
        # Make connections
        buttonPerfom.clicked.connect( self.performAction )
        # Update TreeView, ActionList
        self.updateTree()
        self.updateActionList([('Invert selection',  self.list.invertSelection)])

    def getAllMotorNames(self):
        nameSet = set()
        for dic in self.motorsList:
            for key in dic.keys(): nameSet.add( key )
        return list( nameSet )

    def performAction(self):
#        items = self.list.selectedItems()
#        for item in items:
#            print item.text(0)
        index = self.actionCBox.currentIndex()
        self.actionList[index].trigger()
        
    def updateActionList(self,  functionList):
#        self.actionList = []
        self.actionCBox.clear()
        for (name, function) in functionList:
            fName = name if name is not '' else function.func_name
            act = qt.QAction( fName,  self )
            act.triggered.connect( function )
            self.actionList += [ act ]
            self.actionCBox.addItem( fName )

    def updatePlots(self,  newLegends,  newMotorValues):
        self.legendList = newLegends
        self.motorsList = newMotorValues
        self.motorNamesList = [''] + self.getAllMotorNames()
        self.motorNamesList.sort()
        for cBox in self.cBoxList:
            index = cBox.currentIndex()
            cBox.clear()
            cBox.addItems( self.motorNamesList )
            cBox.setCurrentIndex(index)
        self.updateTree()

    def updateTree(self):
        self.list.clear()
        mList = [ str( cBox.currentText() ) for cBox in self.cBoxList ]
        labels = ['Legend'] + mList
        self.list.setHeaderLabels(labels)
        for i in range(len(self.legendList)):
            legend = self.legendList[i]
            values = self.motorsList[i]
            tmp = qt.QStringList( legend )
            for m in mList:
                if m == '':
                    tmp.append( '' )
                else:
                    tmp.append( str(values.get( m, '---' )) )
            item = qt.QTreeWidgetItem(self.list,  tmp  )
        for i in range(self.list.columnCount()):
            self.list.resizeColumnToContents(i)
 
def main():
    import sys
    legends = ['Curve0', 'Curve1', 'Curve2']
    motors = [{'Motor12': 1, 'Motor11': 8.692713996985609, 'Motor10': 21.98364185388587, 'Motor 8': 0.19806882661182112, 'Motor 9': 0.4844754557916431, 'Motor 4': 0.3502522172639875, 'Motor 5': 0.6639252709334457, 'Motor 6': 0.8130332644206067, 'Motor 7': 0.22114941021809853, 'Motor 0': 0.5931882588655031, 'Motor 1': 0.6780103928805297, 'Motor 2': 0.26738924783290086, 'Motor 3': 0.6778906178576761}, {'Motor18': 0.4707468826876532, 'Motor17': 0.6958160702991127, 'Motor16': 0.8257808117546283, 'Motor15': 0.2587637453100148, 'Motor14': 0.7392644674355958, 'Motor13': 0.09084289261899736, 'Motor12': 2, 'Motor11': 0.21344565983311958, 'Motor10': 0.823400550314221, 'Motor 8': 0.020278096856981342, 'Motor 9': 0.5687440213219551, 'Motor 4': 0.8537811553701731, 'Motor 5': 0.6967303868907243, 'Motor 6': 0.2691963139564302, 'Motor 7': 0.7932933343951395, 'Motor 0': 0.7692165677566825, 'Motor 1': 0.9590927095265979, 'Motor 2': 0.010926468369733544, 'Motor 3': 0.5382649725528551}, {'Motor12': 2, 'Motor11': 0.44400576643956124, 'Motor10': 0.613870067851634, 'Motor 8': 0.901968648110583, 'Motor 9': 0.3197687710845185, 'Motor 4': 0.5714322786278168, 'Motor 5': 0.2786758361634877, 'Motor 6': 0.15443677487828655, 'Motor 7': 0.41623199933237764, 'Motor 0': 0.294201017230741, 'Motor 1': 0.813913587747513, 'Motor 2': 0.5775729031053222, 'Motor 3': 0.8690451825680668}, {'Motor13': 0.6491598094029021, 'Motor12': 10, 'Motor11': 0.006312468992195397, 'Motor10': 0.06727805971206435, 'Motor 8': 0.0929878987747117, 'Motor 9': 0.014325738753558803, 'Motor 4': 0.8185362197656616, 'Motor 5': 0.6643614796103005, 'Motor 6': 0.6479279384366304, 'Motor 7': 0.3485172683358245, 'Motor 0': 0.9858738343685299, 'Motor 1': 0.9330130170323839, 'Motor 2': 0.7550180320112966, 'Motor 3': 0.8814284215685484}, {'Motor19': 0.39846564175862953, 'Motor18': 0.2745751180457152, 'Motor17': 0.42793840508599434, 'Motor16': 0.5335910248322966, 'Motor15': 0.14010423968992758, 'Motor14': 0.27948624022431734, 'Motor13': 0.1737756266389101, 'Motor12': 0.6425110521350722, 'Motor11': 0.9040646490476784, 'Motor10': 0.22997142790156133, 'Motor 8': 0.3520106476992403, 'Motor 9': 0.37023110928070235, 'Motor 4': 0.8110924828319052, 'Motor 5': 0.854155188450653, 'Motor 6': 0.12438157550841666, 'Motor 7': 0.3303770832430888, 'Motor 0': 0.4583273673870403, 'Motor 1': 0.40863603059350373, 'Motor 2': 0.7396799985670546, 'Motor 3': 0.5532134465740317, 'Motor22': 0.7154261407207922, 'Motor20': 0.6735594219326284, 'Motor21': 0.24068704947080943}, {'Motor18': 0.7501922139242619, 'Motor17': 0.067572631661458, 'Motor16': 0.23941863624378346, 'Motor15': 0.543195970137226, 'Motor14': 0.5045110454536483, 'Motor13': 0.47129338234441986, 'Motor12': 0.7039345533241258, 'Motor11': 0.5496976809598649, 'Motor10': 0.028685484457880994, 'Motor 8': 0.3736138811685542, 'Motor 9': 0.6200990287805606, 'Motor 4': 0.30138047598948403, 'Motor 5': 0.15683187764664286, 'Motor 6': 0.061169736595949264, 'Motor 7': 0.35931932492621954, 'Motor 0': 0.7241839150429988, 'Motor 1': 0.7985803970529565, 'Motor 2': 0.5239059568843569, 'Motor 3': 0.7404964999807312}, {'Motor10': 0.90828582481094, 'Motor 8': 0.8424405354748069, 'Motor 9': 0.021278797555318363, 'Motor 4': 0.8593234401902958, 'Motor 5': 0.2638651881043157, 'Motor 6': 0.281687767263718, 'Motor 7': 0.48283570902507555, 'Motor 0': 0.659487116102895, 'Motor 1': 24.591253182578376, 'Motor 2': 3.032078904732739, 'Motor 3': 0.17860013910027928}, {'Motor 8': 0.7246181445974952, 'Motor 9': 0.5375876404160089, 'Motor 4': 0.7608877399780997, 'Motor 5': 0.6164359666836775, 'Motor 6': 0.3910546574315933, 'Motor 7': 0.5287834048239588, 'Motor 0': 0.9700467881758079, 'Motor 1': 0.9064128957850547, 'Motor 2': 0.4434306640093745, 'Motor 3': 0.2783396189782661}, {'Motor19': 0.4741833534896892, 'Motor18': 0.1884371839846597, 'Motor17': 0.660882814263354, 'Motor16': 0.25871486157318313, 'Motor15': 0.6181192138005907, 'Motor14': 0.11534451504645371, 'Motor13': 0.3356756251510249, 'Motor12': 0.8578128852052718, 'Motor11': 0.002943123668270098, 'Motor10': 0.08980970319869397, 'Motor 8': 0.40586648583549123, 'Motor 9': 0.7700310455423328, 'Motor 4': 0.8389920867382025, 'Motor 5': 0.2560110245056251, 'Motor 6': 0.671297941874289, 'Motor 7': 0.7041220063735543, 'Motor 0': 0.4865107750866541, 'Motor 1': 0.8623573559114868, 'Motor 2': 0.8378911209243649, 'Motor 3': 0.056056301247044416, 'Motor24': 0.8535082807686701, 'Motor22': 0.4362354327544248, 'Motor23': 0.17386904782647783, 'Motor20': 0.11001296204329247, 'Motor21': 0.5653716280128318}, {'Motor13': 0.5900826517637087, 'Motor12': 0.2876746207456713, 'Motor11': 0.1829075413610104, 'Motor10': 0.9677552520998641, 'Motor 8': 0.47506344789108046, 'Motor 9': 0.32097198197020305, 'Motor 4': 0.5708449042766175, 'Motor 5': 0.06093583375842648, 'Motor 6': 0.10172375432338043, 'Motor 7': 0.989917381621416, 'Motor 0': 0.8047039621208083, 'Motor 1': 0.9477209087673744, 'Motor 2': 0.46582818765280054, 'Motor 3': 0.0511893987634543},  {'Motor12': 0.6504890336783156, 'Motor11': 0.44400576643956124, 'Motor10': 0.613870067851634, 'Motor 8': 0.901968648110583, 'Motor 9': 0.3197687710845185, 'Motor 4': 0.5714322786278168, 'Motor 5': 0.2786758361634877, 'Motor 6': 0.15443677487828655, 'Motor 7': 0.41623199933237764, 'Motor 0': 0.294201017230741, 'Motor 1': 0.813913587747513, 'Motor 2': 0.5775729031053222, 'Motor 3': 0.8690451825680668}, {'Motor13': 0.6491598094029021, 'Motor12': 0.2975843286841311, 'Motor11': 0.006312468992195397, 'Motor10': 0.06727805971206435, 'Motor 8': 0.0929878987747117, 'Motor 9': 0.014325738753558803, 'Motor 4': 0.8185362197656616, 'Motor 5': 0.6643614796103005, 'Motor 6': 0.6479279384366304, 'Motor 7': 0.3485172683358245, 'Motor 0': 0.9858738343685299, 'Motor 1': 0.9330130170323839, 'Motor 2': 0.7550180320112966, 'Motor 3': 0.8814284215685484},  {'Motor12': 0.6504890336783156, 'Motor11': 0.44400576643956124, 'Motor10': 0.613870067851634, 'Motor 8': 0.901968648110583, 'Motor 9': 0.3197687710845185, 'Motor 4': 0.5714322786278168, 'Motor 5': 0.2786758361634877, 'Motor 6': 0.15443677487828655, 'Motor 7': 0.41623199933237764, 'Motor 0': 0.294201017230741, 'Motor 1': 0.813913587747513, 'Motor 2': 0.5775729031053222, 'Motor 3': 0.8690451825680668}, {'Motor13': 0.6491598094029021, 'Motor12': 0.2975843286841311, 'Motor11': 0.006312468992195397, 'Motor10': 0.06727805971206435, 'Motor 8': 0.0929878987747117, 'Motor 9': 0.014325738753558803, 'Motor 4': 0.8185362197656616, 'Motor 5': 0.6643614796103005, 'Motor 6': 0.6479279384366304, 'Motor 7': 0.3485172683358245, 'Motor 0': 0.9858738343685299, 'Motor 1': 0.9330130170323839, 'Motor 2': 0.7550180320112966, 'Motor 3': 0.8814284215685484}]
    app = qt.QApplication(sys.argv)
#    w = SortPlotsWidget(None, legends,  motors)
    w = SortPlotsWidget(None, legends,  motors)
    w.show()
    app.exec_()
    
if __name__ == '__main__':
    main()
