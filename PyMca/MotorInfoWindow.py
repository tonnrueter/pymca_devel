import weakref,  time
from PyMca import PyMcaQt as qt

DEBUG = 1

class MotorInfoComboBox(qt.QComboBox):
    
    loadColumnSignal = qt.pyqtSignal(object)
    
    def __init__(self,  parent,  mlist,  nCol):
        qt.QComboBox.__init__(self,  parent)
        self.setSizeAdjustPolicy( qt.QComboBox.AdjustToContents )
        self.motorNamesList = ["       "] + mlist
        self.nColumn = nCol
        self.addItems( [ qt.QString(elem) for elem in self.motorNamesList ] )
        self.activated.connect(self.emitLoadColumnSignal)

    def emitLoadColumnSignal(self):
        ddict = {}
        ddict['column'] = self.nColumn
        ddict['motor'] = str(self.currentText())
        ddict['event'] = "activated"
        self.loadColumnSignal.emit(ddict)

    def currentMotor(self):
        return str( self.currentText() )
        
    def updateMotorNamesList(self,  newMotorNamesList):
        currentMotorName = self.currentMotor()
        self.clear()
        newMotorNamesList = ['       '] + newMotorNamesList
        self.motorNamesList = newMotorNamesList 
        if currentMotorName in self.motorNamesList:
            newIndex = newMotorNamesList.index(currentMotorName)
        else:
            newIndex = 0
        self.addItems( [ qt.QString(elem) for elem in self.motorNamesList ] )
        self.setCurrentIndex(newIndex)

class MotorInfoDialog(qt.QWidget):
    def __init__(self,  parent,  legends,  motorValues):
        """
        legends            List contains Plotnames
        motorValues     List contains names and values of the motors
        """
        qt.QWidget.__init__(self,  parent)
        self.setWindowTitle("Motor Info Plugin")
        if len(legends) <> len(motorValues):
            print('Consistency error: legends and motorValues do not have same length!')
        self.numCurves = len(legends)
        # Buttons
        self.buttonAddColumn = qt.QPushButton("Add",  self)
        self.buttonDeleteColumn = qt.QPushButton("Del",  self)
        self.buttonClose = qt.QPushButton("Close",  self)
        self.buttonUpdate = qt.QPushButton("Update",  self)
        # Table
        self.table = MotorInfoTable(self, self.numCurves + 1, 5, legends,  motorValues)
        # Layout
        self.mainLayout = qt.QGridLayout(self)
        self.mainLayout.setContentsMargins(1, 1, 1, 1)
        self.mainLayout.setSpacing(2)
        self.buttonLayout = qt.QHBoxLayout(None)
        self.buttonLayout .setSpacing(1)
        # Add widgets to layour
        self.mainLayout.addWidget(self.table,  0,  0)
        self.mainLayout.addLayout(self.buttonLayout,  1,  0)
        self.buttonLayout.addWidget(self.buttonUpdate)
        self.buttonLayout.addWidget(self.buttonAddColumn)
        self.buttonLayout.addWidget(self.buttonDeleteColumn)
        self.buttonLayout.addWidget(qt.HorizontalSpacer(self))
        self.buttonLayout.addWidget(self.buttonClose)
        self.resize(700,  500)
        # Make connections
        self.buttonClose.clicked.connect(self.close)
        self.buttonAddColumn.clicked.connect(self.table.addColumn)
        self.buttonDeleteColumn.clicked.connect(self.table.delColumn)
        # buttonUpdate is connected in MotorInfoPlugin

    def keyPressEvent(self,  event):
        if (event.key() == qt.Qt.Key_Escape):
            self.close()

class MotorInfoTable(qt.QTableWidget):
    def __init__(self,  parent,  numRows,  numColumns,  legList, motList):
        qt.QTableWidget.__init__(self,  1, 1,  parent)
        self.currentComboBox = 1
        self.legendsList = legList
        self.motorsList  = motList
        self.motorNamesList = self.getAllMotorNames()
        self.motorNamesList.sort()
        self.infoDict = dict( zip( self.legendsList,  self.motorsList ) )
        self.horizHeader = self.horizontalHeader()
        self.horizHeader.setResizeMode(qt.QHeaderView.ResizeToContents)
        self.horizHeader.hide()
        self.verticalHeader().hide()
        self.setItem(0, 0, self.textItem("Plot Name"))
        for i in range(1,  numColumns):
            self.addColumn()
        for curveLegend in self.legendsList:
            self.insertRow(self.rowCount())
            self.setItem(self.rowCount()-1, 0,  self.textItem( curveLegend ))

    def addColumn(self):
        currentColumn = self.columnCount()
        self.insertColumn(currentColumn)
        newBox = MotorInfoComboBox(self,  self.motorNamesList ,  currentColumn)
        newBox.loadColumnSignal.connect(self.loadColumn)
        self.setCellWidget(0, currentColumn, newBox)

    def delColumn(self):
        if self.columnCount() > 1:
            self.removeColumn(self.columnCount()-1)

    def addRow(self):
        currentRow = self.rowCount()
        self.insertRow(currentRow)
        legend = self.legendsList[currentRow-1]
        self.setItem(currentRow, 0,  self.textItem( legend ))
        for i in range(1,  self.columnCount()):
            cBox = self.cellWidget(0,  i)
            motorName = cBox.currentMotor()
            motors = self.motorsList[currentRow-1]
            value = motors.get(motorName,  None)
            if motorName == "       ":
                self.setItem(currentRow, i,  self.textItem( '' ))
            elif value is not None:
                self.setItem(currentRow, i,  self.textItem( str(value) ))
            else:
                self.setItem(currentRow, i,  self.textItem( '---' ))

    def delRow(self,  n):
        if n > 0:
            self.removeRow( n )

    def updateTable(self,  legList, motList):
        if self.legendsList == legList and self.motorsList == motList:
            pass
        else:
            for i in range(0, self.rowCount()):
                self.delRow( self.rowCount()-1 )
            self.legendsList = legList
            self.motorsList = motList
            self.infoDict = dict( zip( self.legendsList,  self.motorsList ) )
            self.motorNamesList = self.getAllMotorNames()
            self.motorNamesList.sort()
            for i in range(1, self.columnCount()):
                cBox = self.cellWidget(0, i)
                cBox.updateMotorNamesList(self.motorNamesList)
            for i in range(len(legList)):
                self.addRow()

    def loadColumn(self,  ddict):
        for key in ddict.keys():
            if str(key) == str("motor"):
                motorName = ddict[key]
            elif str(key) == str("column"):
                column = ddict[key]
        if motorName <> "       ":
            for i in range(1,  self.rowCount()):
                legend = str( self.item(i, 0).text() )
                curveInfo = self.infoDict.get(legend,  None)
                if curveInfo is not None:
                    motorValue = curveInfo.get(motorName,  '---')
                else:
                    motorValue = '---'
                self.setItem(i,  column,  self.textItem( str(motorValue) ))
        else:
            for i in range(1,  self.rowCount()):
                self.setItem(i,  column,  self.textItem( "" ))

    def getAllMotorNames(self):
        nameSet = set()
        for dic in self.motorsList:
            for key in dic.keys(): nameSet.add( key )
        return list( nameSet )

    def textItem(self,  text):
        item = qt.QTableWidgetItem(text)
        item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
        return item

def generateDummies(n,  m):
    dummyDict = {}
    legendList = [("Curve" + str(i)) for i in range(n)]
    motorsList = []
    for elem in legendList:
        randInt = random.randint(1,1+m)
        tmp = dict( [(('Motor%2d'%i), random.random()) for i in range(0, randInt)] )
        motorsList  += [tmp]
    return legendList,  motorsList

def updateTest(w):
#    legends = ['Curve1', 'Curve2', 'Curve3', 'Curve4', 'Curve5', 'Curve6', 'Curve7', 'Curve8']
    legends = ['Curve1']
#    motors = [{'Motor18': 0.4707468826876532, 'Motor17': 0.6958160702991127, 'Motor16': 0.8257808117546283, 'Motor15': 0.2587637453100148, 'Motor14': 0.7392644674355958, 'Motor13': 0.09084289261899736, 'Motor12': 0.5190253643331453, 'Motor11': 0.21344565983311958, 'Motor10': 0.823400550314221, 'Motor 8': 0.020278096856981342, 'Motor 9': 0.5687440213219551, 'Motor 4': 0.8537811553701731, 'Motor 5': 0.6967303868907243, 'Motor 6': 0.2691963139564302, 'Motor 7': 0.7932933343951395, 'Motor 0': 0.7692165677566825, 'Motor 1': 0.9590927095265979, 'Motor 2': 0.010926468369733544, 'Motor 3': 0.5382649725528551}, {'Motor12': 0.6504890336783156, 'Motor11': 0.44400576643956124, 'Motor10': 0.613870067851634, 'Motor 8': 0.901968648110583, 'Motor 9': 0.3197687710845185, 'Motor 4': 0.5714322786278168, 'Motor 5': 0.2786758361634877, 'Motor 6': 0.15443677487828655, 'Motor 7': 0.41623199933237764, 'Motor 0': 0.294201017230741, 'Motor 1': 0.813913587747513, 'Motor 2': 0.5775729031053222, 'Motor 3': 0.8690451825680668}, {'Motor13': 0.6491598094029021, 'Motor12': 0.2975843286841311, 'Motor11': 0.006312468992195397, 'Motor10': 0.06727805971206435, 'Motor 8': 0.0929878987747117, 'Motor 9': 0.014325738753558803, 'Motor 4': 0.8185362197656616, 'Motor 5': 0.6643614796103005, 'Motor 6': 0.6479279384366304, 'Motor 7': 0.3485172683358245, 'Motor 0': 0.9858738343685299, 'Motor 1': 0.9330130170323839, 'Motor 2': 0.7550180320112966, 'Motor 3': 0.8814284215685484}, {'Motor19': 0.39846564175862953, 'Motor18': 0.2745751180457152, 'Motor17': 0.42793840508599434, 'Motor16': 0.5335910248322966, 'Motor15': 0.14010423968992758, 'Motor14': 0.27948624022431734, 'Motor13': 0.1737756266389101, 'Motor12': 0.6425110521350722, 'Motor11': 0.9040646490476784, 'Motor10': 0.22997142790156133, 'Motor 8': 0.3520106476992403, 'Motor 9': 0.37023110928070235, 'Motor 4': 0.8110924828319052, 'Motor 5': 0.854155188450653, 'Motor 6': 0.12438157550841666, 'Motor 7': 0.3303770832430888, 'Motor 0': 0.4583273673870403, 'Motor 1': 0.40863603059350373, 'Motor 2': 0.7396799985670546, 'Motor 3': 0.5532134465740317, 'Motor22': 0.7154261407207922, 'Motor20': 0.6735594219326284, 'Motor21': 0.24068704947080943}, {'Motor18': 0.7501922139242619, 'Motor17': 0.067572631661458, 'Motor16': 0.23941863624378346, 'Motor15': 0.543195970137226, 'Motor14': 0.5045110454536483, 'Motor13': 0.47129338234441986, 'Motor12': 0.7039345533241258, 'Motor11': 0.5496976809598649, 'Motor10': 0.028685484457880994, 'Motor 8': 0.3736138811685542, 'Motor 9': 0.6200990287805606, 'Motor 4': 0.30138047598948403, 'Motor 5': 0.15683187764664286, 'Motor 6': 0.061169736595949264, 'Motor 7': 0.35931932492621954, 'Motor 0': 0.7241839150429988, 'Motor 1': 0.7985803970529565, 'Motor 2': 0.5239059568843569, 'Motor 3': 0.7404964999807312}, {'Motor10': 0.90828582481094, 'Motor 8': 0.8424405354748069, 'Motor 9': 0.021278797555318363, 'Motor 4': 0.8593234401902958, 'Motor 5': 0.2638651881043157, 'Motor 6': 0.281687767263718, 'Motor 7': 0.48283570902507555, 'Motor 0': 0.659487116102895, 'Motor 1': 0.24591253182578376, 'Motor 2': 0.3032078904732739, 'Motor 3': 0.17860013910027928}, {'Motor 8': 0.7246181445974952, 'Motor 9': 0.5375876404160089, 'Motor 4': 0.7608877399780997, 'Motor 5': 0.6164359666836775, 'Motor 6': 0.3910546574315933, 'Motor 7': 0.5287834048239588, 'Motor 0': 0.9700467881758079, 'Motor 1': 0.9064128957850547, 'Motor 2': 0.4434306640093745, 'Motor 3': 0.2783396189782661}, {'Motor19': 0.4741833534896892, 'Motor18': 0.1884371839846597, 'Motor17': 0.660882814263354, 'Motor16': 0.25871486157318313, 'Motor15': 0.6181192138005907, 'Motor14': 0.11534451504645371, 'Motor13': 0.3356756251510249, 'Motor12': 0.8578128852052718, 'Motor11': 0.002943123668270098, 'Motor10': 0.08980970319869397, 'Motor 8': 0.40586648583549123, 'Motor 9': 0.7700310455423328, 'Motor 4': 0.8389920867382025, 'Motor 5': 0.2560110245056251, 'Motor 6': 0.671297941874289, 'Motor 7': 0.7041220063735543, 'Motor 0': 0.4865107750866541, 'Motor 1': 0.8623573559114868, 'Motor 2': 0.8378911209243649, 'Motor 3': 0.056056301247044416, 'Motor24': 0.8535082807686701, 'Motor22': 0.4362354327544248, 'Motor23': 0.17386904782647783, 'Motor20': 0.11001296204329247, 'Motor21': 0.5653716280128318}]
    motors = [{'Motor New': 0.12345, 'Motor16': 0.8257808117546283, 'Motor15': 0.2587637453100148, 'Motor14': 0.7392644674355958, 'Motor12': 0.5190253643331453, 'Motor11': 0.21344565983311958, 'Motor10': 0.823400550314221, 'Motor 8': 0.020278096856981342, 'Motor 9': 0.5687440213219551, 'Motor 4': 0.8537811553701731, 'Motor 5': 0.6967303868907243, 'Motor 6': 0.2691963139564302, 'Motor 7': 0.7932933343951395, 'Motor 0': 0.7692165677566825, 'Motor 1': 0.9590927095265979, 'Motor 2': 0.010926468369733544, 'Motor 3': 0.5382649725528551}]
    w.table.updateTable(legends,  motors)

def main():
    import sys,  random
#    legends,  motors = generateDummies(10, 20)
    legends = ['Curve0', 'Curve1', 'Curve2', 'Curve3', 'Curve4', 'Curve5', 'Curve6', 'Curve7', 'Curve8', 'Curve9']
    motors = [{'Motor12': 0.5283546103038855, 'Motor11': 0.8692713996985609, 'Motor10': 0.2198364185388587, 'Motor 8': 0.19806882661182112, 'Motor 9': 0.4844754557916431, 'Motor 4': 0.3502522172639875, 'Motor 5': 0.6639252709334457, 'Motor 6': 0.8130332644206067, 'Motor 7': 0.22114941021809853, 'Motor 0': 0.5931882588655031, 'Motor 1': 0.6780103928805297, 'Motor 2': 0.26738924783290086, 'Motor 3': 0.6778906178576761}, {'Motor18': 0.4707468826876532, 'Motor17': 0.6958160702991127, 'Motor16': 0.8257808117546283, 'Motor15': 0.2587637453100148, 'Motor14': 0.7392644674355958, 'Motor13': 0.09084289261899736, 'Motor12': 0.5190253643331453, 'Motor11': 0.21344565983311958, 'Motor10': 0.823400550314221, 'Motor 8': 0.020278096856981342, 'Motor 9': 0.5687440213219551, 'Motor 4': 0.8537811553701731, 'Motor 5': 0.6967303868907243, 'Motor 6': 0.2691963139564302, 'Motor 7': 0.7932933343951395, 'Motor 0': 0.7692165677566825, 'Motor 1': 0.9590927095265979, 'Motor 2': 0.010926468369733544, 'Motor 3': 0.5382649725528551}, {'Motor12': 0.6504890336783156, 'Motor11': 0.44400576643956124, 'Motor10': 0.613870067851634, 'Motor 8': 0.901968648110583, 'Motor 9': 0.3197687710845185, 'Motor 4': 0.5714322786278168, 'Motor 5': 0.2786758361634877, 'Motor 6': 0.15443677487828655, 'Motor 7': 0.41623199933237764, 'Motor 0': 0.294201017230741, 'Motor 1': 0.813913587747513, 'Motor 2': 0.5775729031053222, 'Motor 3': 0.8690451825680668}, {'Motor13': 0.6491598094029021, 'Motor12': 0.2975843286841311, 'Motor11': 0.006312468992195397, 'Motor10': 0.06727805971206435, 'Motor 8': 0.0929878987747117, 'Motor 9': 0.014325738753558803, 'Motor 4': 0.8185362197656616, 'Motor 5': 0.6643614796103005, 'Motor 6': 0.6479279384366304, 'Motor 7': 0.3485172683358245, 'Motor 0': 0.9858738343685299, 'Motor 1': 0.9330130170323839, 'Motor 2': 0.7550180320112966, 'Motor 3': 0.8814284215685484}, {'Motor19': 0.39846564175862953, 'Motor18': 0.2745751180457152, 'Motor17': 0.42793840508599434, 'Motor16': 0.5335910248322966, 'Motor15': 0.14010423968992758, 'Motor14': 0.27948624022431734, 'Motor13': 0.1737756266389101, 'Motor12': 0.6425110521350722, 'Motor11': 0.9040646490476784, 'Motor10': 0.22997142790156133, 'Motor 8': 0.3520106476992403, 'Motor 9': 0.37023110928070235, 'Motor 4': 0.8110924828319052, 'Motor 5': 0.854155188450653, 'Motor 6': 0.12438157550841666, 'Motor 7': 0.3303770832430888, 'Motor 0': 0.4583273673870403, 'Motor 1': 0.40863603059350373, 'Motor 2': 0.7396799985670546, 'Motor 3': 0.5532134465740317, 'Motor22': 0.7154261407207922, 'Motor20': 0.6735594219326284, 'Motor21': 0.24068704947080943}, {'Motor18': 0.7501922139242619, 'Motor17': 0.067572631661458, 'Motor16': 0.23941863624378346, 'Motor15': 0.543195970137226, 'Motor14': 0.5045110454536483, 'Motor13': 0.47129338234441986, 'Motor12': 0.7039345533241258, 'Motor11': 0.5496976809598649, 'Motor10': 0.028685484457880994, 'Motor 8': 0.3736138811685542, 'Motor 9': 0.6200990287805606, 'Motor 4': 0.30138047598948403, 'Motor 5': 0.15683187764664286, 'Motor 6': 0.061169736595949264, 'Motor 7': 0.35931932492621954, 'Motor 0': 0.7241839150429988, 'Motor 1': 0.7985803970529565, 'Motor 2': 0.5239059568843569, 'Motor 3': 0.7404964999807312}, {'Motor10': 0.90828582481094, 'Motor 8': 0.8424405354748069, 'Motor 9': 0.021278797555318363, 'Motor 4': 0.8593234401902958, 'Motor 5': 0.2638651881043157, 'Motor 6': 0.281687767263718, 'Motor 7': 0.48283570902507555, 'Motor 0': 0.659487116102895, 'Motor 1': 0.24591253182578376, 'Motor 2': 0.3032078904732739, 'Motor 3': 0.17860013910027928}, {'Motor 8': 0.7246181445974952, 'Motor 9': 0.5375876404160089, 'Motor 4': 0.7608877399780997, 'Motor 5': 0.6164359666836775, 'Motor 6': 0.3910546574315933, 'Motor 7': 0.5287834048239588, 'Motor 0': 0.9700467881758079, 'Motor 1': 0.9064128957850547, 'Motor 2': 0.4434306640093745, 'Motor 3': 0.2783396189782661}, {'Motor19': 0.4741833534896892, 'Motor18': 0.1884371839846597, 'Motor17': 0.660882814263354, 'Motor16': 0.25871486157318313, 'Motor15': 0.6181192138005907, 'Motor14': 0.11534451504645371, 'Motor13': 0.3356756251510249, 'Motor12': 0.8578128852052718, 'Motor11': 0.002943123668270098, 'Motor10': 0.08980970319869397, 'Motor 8': 0.40586648583549123, 'Motor 9': 0.7700310455423328, 'Motor 4': 0.8389920867382025, 'Motor 5': 0.2560110245056251, 'Motor 6': 0.671297941874289, 'Motor 7': 0.7041220063735543, 'Motor 0': 0.4865107750866541, 'Motor 1': 0.8623573559114868, 'Motor 2': 0.8378911209243649, 'Motor 3': 0.056056301247044416, 'Motor24': 0.8535082807686701, 'Motor22': 0.4362354327544248, 'Motor23': 0.17386904782647783, 'Motor20': 0.11001296204329247, 'Motor21': 0.5653716280128318}, {'Motor13': 0.5900826517637087, 'Motor12': 0.2876746207456713, 'Motor11': 0.1829075413610104, 'Motor10': 0.9677552520998641, 'Motor 8': 0.47506344789108046, 'Motor 9': 0.32097198197020305, 'Motor 4': 0.5708449042766175, 'Motor 5': 0.06093583375842648, 'Motor 6': 0.10172375432338043, 'Motor 7': 0.989917381621416, 'Motor 0': 0.8047039621208083, 'Motor 1': 0.9477209087673744, 'Motor 2': 0.46582818765280054, 'Motor 3': 0.0511893987634543}]
    app = qt.QApplication(sys.argv)
    w = MotorInfoDialog(None, legends,  motors)
    w.show()
    app.exec_()
    
if __name__ == '__main__':
    main()
