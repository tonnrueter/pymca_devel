import sys,  random
from PyQt4 import QtGui

from PyMca import PyMcaQt as qt
from PyMca.PyMca_Icons import IconDict
#try:
#    from PyMca import MotorInfo
#except:
#    print("WARNING: MotorInfo performing local import")
#    from . import MotorInfo

DEBUG = 0

class MotorInfoComboBox(qt.QComboBox):
    def __init__(self,  parent,  mlist,  nCol):
        qt.QComboBox.__init__(self,  parent)
        self.table = parent
        self.motorNamesList = [""] + mlist
        self.nColumn = nCol
        self.infoDict = dict( zip( parent.legendsList,  parent.motorsList ) )
        self.addItems( [ qt.QString(elem) for elem in self.motorNamesList ] )

    def _getColumn(self):
        return self.nColumn

    def _loadColumn(self,  key):
        if DEBUG:
            print("-- Load column %d with key '%s'"  % (self.nColumn,  str(key)))
        """
        - Get curve legend from first Column -> rowLegend
        - Access dictionary entry for rowLegend and get motorName
        """
        if key <> 0:
            motorName = self.motorNamesList[key]
            for i in range(1,  self.table.rowCount()):
                legend = str( self.table.item(i, 0).text() ) # Problem: QTableWidgetItem -> QString -> String ?!
                curveInfo = self.infoDict.get(legend,  None)
                if curveInfo is not None:
                    motorValue = curveInfo.get(motorName,  '---')
                else:
                    motorValue = '---'
                if DEBUG:
                    print("---- '%s' has motorValue '%s'"  % (legend,  motorValue))
                self.table.setItem(i,  self.nColumn,  self.table._textItem( str(motorValue) ))
        else:
            for i in range(1,  self.table.rowCount()):
                self.table.setItem(i,  self.nColumn,  qt.QTableWidgetItem( str('') ))

    
class MotorInfoWidget(qt.QDialog): # Changed from QDialog
    def __init__(self,  parent,  legends,  motorValues):
        # legends and motorValues are both Lists
        qt.QWidget.__init__(self,  parent)
        self.setWindowTitle("Motor Info Plugin")
        if len(legends) <> len(motorValues):
            print('!!! Consistency error !!!')
        self.numCurves = len(legends)
        self.numMotors = 4
        
        # Define Buttons
        self.buttonAddColumn = qt.QPushButton("Add",  self)
#        self.buttonAddColumn.setAutoDefault(False)
        self.buttonDeleteColumn = qt.QPushButton("Del",  self)
        self.buttonClose = qt.QPushButton("Close")
        
        # Define Table
        self.table = MotorInfoTable(self, self.numCurves + 1, self.numMotors + 1, legends,  motorValues)

        # Layout Section
        self.layout = qt.QGridLayout(self)
        self.layout.setContentsMargins(1, 1, 1, 1)
        self.layout.setSpacing(2)
        self.buttonLayout = qt.QHBoxLayout(None)
        self.buttonLayout .setSpacing(1)

        # Add Widgets to Layout
        self.layout.addWidget(self.table,  0,  0)
        self.layout.addLayout(self.buttonLayout,  1,  0)
        self.buttonLayout.addWidget(self.buttonAddColumn)
        self.buttonLayout.addWidget(self.buttonDeleteColumn)
        self.buttonLayout.addWidget(qt.HorizontalSpacer(self))
        self.buttonLayout.addWidget(self.buttonClose)
        
        self.resize(700,  500)
        
        # Make connections
        self.buttonClose.clicked.connect(self.close)
        self.buttonAddColumn.clicked.connect(self.table._addColumn)
        self.buttonDeleteColumn.clicked.connect(self.table._delColumn)
        
    def keyPressEvent(self,  event):
        if (event.key() == qt.Qt.Key_Escape):
            self.close()

class MotorInfoTable(qt.QTableWidget):
    def __init__(self,  parent,  numRows,  numColumns,  legList, motList):
        qt.QTableWidget.__init__(self,  1, 1,  parent)
        self.legendsList = legList
        self.motorsList  = motList
        self.motorNamesList = self._getAllMotorNames()
        self.motorNamesList.sort()
        if DEBUG:
            print("Created new table")
        # Build table
        # Remove real headers
        self.horizHeader = self.horizontalHeader()
        self.horizHeader.setResizeMode(qt.QHeaderView.ResizeToContents)
        self.horizHeader.sectionClicked.connect(self.sortByColumn)
        self.horizHeader.hide()
        self.verticalHeader().hide()
        self.setItem(0, 0, self._textItem("Plot Name"))
        # Add columns
        for i in range(1,  numColumns):
            self._addColumn()
        # Add rows
        for curveLegend in self.legendsList:
            self.insertRow(self.rowCount())
            self.setItem(self.rowCount()-1, 0,  self._textItem( curveLegend ))

    def _addColumn(self):
        currentColumn = self.columnCount()
        self.insertColumn(currentColumn)
        newBox = MotorInfoComboBox(self,  self.motorNamesList ,  currentColumn)
        newBox.activated.connect(newBox._loadColumn) # New Style Connect
        self.setCellWidget(0, currentColumn, newBox)

    def _delColumn(self):
        if self.columnCount() > 1:
            self.removeColumn(self.columnCount()-1)

    def _getAllMotorNames(self):
        nameSet = set()
        for dic in self.motorsList:
            for key in dic.keys(): nameSet.add( key )
        return list( nameSet )
        
    def _textItem(self,  text):
        # Return QTableWidgetItem
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

def main():
#    legends,  motors = generateDummies(10, 20)
    legends = ['Curve0', 'Curve1', 'Curve2', 'Curve3', 'Curve4', 'Curve5', 'Curve6', 'Curve7', 'Curve8', 'Curve9']
    motors = [{'Motor12': 0.5283546103038855, 'Motor11': 0.8692713996985609, 'Motor10': 0.2198364185388587, 'Motor 8': 0.19806882661182112, 'Motor 9': 0.4844754557916431, 'Motor 4': 0.3502522172639875, 'Motor 5': 0.6639252709334457, 'Motor 6': 0.8130332644206067, 'Motor 7': 0.22114941021809853, 'Motor 0': 0.5931882588655031, 'Motor 1': 0.6780103928805297, 'Motor 2': 0.26738924783290086, 'Motor 3': 0.6778906178576761}, {'Motor18': 0.4707468826876532, 'Motor17': 0.6958160702991127, 'Motor16': 0.8257808117546283, 'Motor15': 0.2587637453100148, 'Motor14': 0.7392644674355958, 'Motor13': 0.09084289261899736, 'Motor12': 0.5190253643331453, 'Motor11': 0.21344565983311958, 'Motor10': 0.823400550314221, 'Motor 8': 0.020278096856981342, 'Motor 9': 0.5687440213219551, 'Motor 4': 0.8537811553701731, 'Motor 5': 0.6967303868907243, 'Motor 6': 0.2691963139564302, 'Motor 7': 0.7932933343951395, 'Motor 0': 0.7692165677566825, 'Motor 1': 0.9590927095265979, 'Motor 2': 0.010926468369733544, 'Motor 3': 0.5382649725528551}, {'Motor12': 0.6504890336783156, 'Motor11': 0.44400576643956124, 'Motor10': 0.613870067851634, 'Motor 8': 0.901968648110583, 'Motor 9': 0.3197687710845185, 'Motor 4': 0.5714322786278168, 'Motor 5': 0.2786758361634877, 'Motor 6': 0.15443677487828655, 'Motor 7': 0.41623199933237764, 'Motor 0': 0.294201017230741, 'Motor 1': 0.813913587747513, 'Motor 2': 0.5775729031053222, 'Motor 3': 0.8690451825680668}, {'Motor13': 0.6491598094029021, 'Motor12': 0.2975843286841311, 'Motor11': 0.006312468992195397, 'Motor10': 0.06727805971206435, 'Motor 8': 0.0929878987747117, 'Motor 9': 0.014325738753558803, 'Motor 4': 0.8185362197656616, 'Motor 5': 0.6643614796103005, 'Motor 6': 0.6479279384366304, 'Motor 7': 0.3485172683358245, 'Motor 0': 0.9858738343685299, 'Motor 1': 0.9330130170323839, 'Motor 2': 0.7550180320112966, 'Motor 3': 0.8814284215685484}, {'Motor19': 0.39846564175862953, 'Motor18': 0.2745751180457152, 'Motor17': 0.42793840508599434, 'Motor16': 0.5335910248322966, 'Motor15': 0.14010423968992758, 'Motor14': 0.27948624022431734, 'Motor13': 0.1737756266389101, 'Motor12': 0.6425110521350722, 'Motor11': 0.9040646490476784, 'Motor10': 0.22997142790156133, 'Motor 8': 0.3520106476992403, 'Motor 9': 0.37023110928070235, 'Motor 4': 0.8110924828319052, 'Motor 5': 0.854155188450653, 'Motor 6': 0.12438157550841666, 'Motor 7': 0.3303770832430888, 'Motor 0': 0.4583273673870403, 'Motor 1': 0.40863603059350373, 'Motor 2': 0.7396799985670546, 'Motor 3': 0.5532134465740317, 'Motor22': 0.7154261407207922, 'Motor20': 0.6735594219326284, 'Motor21': 0.24068704947080943}, {'Motor18': 0.7501922139242619, 'Motor17': 0.067572631661458, 'Motor16': 0.23941863624378346, 'Motor15': 0.543195970137226, 'Motor14': 0.5045110454536483, 'Motor13': 0.47129338234441986, 'Motor12': 0.7039345533241258, 'Motor11': 0.5496976809598649, 'Motor10': 0.028685484457880994, 'Motor 8': 0.3736138811685542, 'Motor 9': 0.6200990287805606, 'Motor 4': 0.30138047598948403, 'Motor 5': 0.15683187764664286, 'Motor 6': 0.061169736595949264, 'Motor 7': 0.35931932492621954, 'Motor 0': 0.7241839150429988, 'Motor 1': 0.7985803970529565, 'Motor 2': 0.5239059568843569, 'Motor 3': 0.7404964999807312}, {'Motor10': 0.90828582481094, 'Motor 8': 0.8424405354748069, 'Motor 9': 0.021278797555318363, 'Motor 4': 0.8593234401902958, 'Motor 5': 0.2638651881043157, 'Motor 6': 0.281687767263718, 'Motor 7': 0.48283570902507555, 'Motor 0': 0.659487116102895, 'Motor 1': 0.24591253182578376, 'Motor 2': 0.3032078904732739, 'Motor 3': 0.17860013910027928}, {'Motor 8': 0.7246181445974952, 'Motor 9': 0.5375876404160089, 'Motor 4': 0.7608877399780997, 'Motor 5': 0.6164359666836775, 'Motor 6': 0.3910546574315933, 'Motor 7': 0.5287834048239588, 'Motor 0': 0.9700467881758079, 'Motor 1': 0.9064128957850547, 'Motor 2': 0.4434306640093745, 'Motor 3': 0.2783396189782661}, {'Motor19': 0.4741833534896892, 'Motor18': 0.1884371839846597, 'Motor17': 0.660882814263354, 'Motor16': 0.25871486157318313, 'Motor15': 0.6181192138005907, 'Motor14': 0.11534451504645371, 'Motor13': 0.3356756251510249, 'Motor12': 0.8578128852052718, 'Motor11': 0.002943123668270098, 'Motor10': 0.08980970319869397, 'Motor 8': 0.40586648583549123, 'Motor 9': 0.7700310455423328, 'Motor 4': 0.8389920867382025, 'Motor 5': 0.2560110245056251, 'Motor 6': 0.671297941874289, 'Motor 7': 0.7041220063735543, 'Motor 0': 0.4865107750866541, 'Motor 1': 0.8623573559114868, 'Motor 2': 0.8378911209243649, 'Motor 3': 0.056056301247044416, 'Motor24': 0.8535082807686701, 'Motor22': 0.4362354327544248, 'Motor23': 0.17386904782647783, 'Motor20': 0.11001296204329247, 'Motor21': 0.5653716280128318}, {'Motor13': 0.5900826517637087, 'Motor12': 0.2876746207456713, 'Motor11': 0.1829075413610104, 'Motor10': 0.9677552520998641, 'Motor 8': 0.47506344789108046, 'Motor 9': 0.32097198197020305, 'Motor 4': 0.5708449042766175, 'Motor 5': 0.06093583375842648, 'Motor 6': 0.10172375432338043, 'Motor 7': 0.989917381621416, 'Motor 0': 0.8047039621208083, 'Motor 1': 0.9477209087673744, 'Motor 2': 0.46582818765280054, 'Motor 3': 0.0511893987634543}]
    app = qt.QApplication(sys.argv)
    w = MotorInfoWidget(None, legends,  motors)
    w.show()
    app.exec_()
    
if __name__ == '__main__':
    main()
