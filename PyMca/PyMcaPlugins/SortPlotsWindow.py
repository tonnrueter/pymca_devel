import numpy
from PyMca import PyMcaQt as qt
from PyMca.PyMca_Icons import IconDict

from PyMca import CloseEventNotifyingWidget
from PyMca import ScanWindow as sw
from PyMca.PyMcaSciPy.signal import median
from PyMca.SGModule import getSavitzkyGolay

DEBUG = True
class SPExtension(qt.QWidget):
    def __init__(self,  parent,  psel = None,  msel = None):
        qt.QWidget.__init__(self, parent)

        if psel:
            self.pselection = psel
        else:
            self.pselection = qt.QStringListModel(self)

        if msel:
            self.mselection = msel
        else:
            self.mselection = qt.QStringListModel(self)
        
        self.pview = qt.QListView(self)
        self.mview = qt.QListView(self)
        
        self.pview.setModel(self.pselection)
        self.mview.setModel(self.mselection)
        
        mainLayout = qt.QHBoxLayout(self)
        mainLayout.addWidget(self.pview)
        mainLayout.addWidget(self.mview)
        self.setLayout(mainLayout)

class CheckBoxWindow(qt.QWidget):
    def __init__(self,  parent):
        qt.QWidget.__init__(self, parent)
       
#        groupBox = qt.QGroupBox('XMCD Options',  self)
        
        self.cbRemoveSpikes = qt.QCheckBox('Remove &spikes', self)
        self.cbNoiseFilter = qt.QCheckBox('Apply noise &filter', self)
        self.cbXAS = qt.QCheckBox('Calculate &XAS', self)
        self.cbNormalize = qt.QCheckBox('&Normalize scans', self)
        
#        widthLayout = qt.QHBoxLayout()
#        self.removeSpikesThreshold = qt.QLineEdit('0.001',  self)
#        self.removeSpikesWidth = qt.QLineEdit('5',  self)
#        thresholdLayout = qt.QHBoxLayout()
#        thresholdLayout.addWidget(qt.QLabel('Threshold'))
#        thresholdLayout.addWidget(qt.HorizontalSpacer())
#        thresholdLayout.addWidget(self.removeSpikesThreshold)
#        widhtLayout = qt.QHBoxLayout()
#        widthLayout.addWidget(qt.QLabel('Width'))
#        widthLayout.addWidget(qt.HorizontalSpacer())
#        widthLayout.addWidget(self.removeSpikesWidth)
        
        mainLayout = qt.QVBoxLayout(self)
        mainLayout.addWidget(self.cbRemoveSpikes)
#        mainLayout.addLayout(widthLayout)
#        mainLayout.addLayout(thresholdLayout)        
        mainLayout.addWidget(self.cbNoiseFilter)
        mainLayout.addWidget(self.cbXAS)
        mainLayout.addWidget(self.cbNormalize)
        self.setLayout(mainLayout)
        
#        groupBox.setLayout(mainLayout)
    
    def check(self):
        return (self.cbRemoveSpikes.isChecked(),
                self.cbNoiseFilter.isChecked(),
                self.cbXAS.isChecked(),
                self.cbNormalize.isChecked())

    def removeSpikeValues(self):
        pass
    

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
        
        if selection != []:
            self.XMCD(selection, 
                      remove_spikes,
                      noise_filter,
                      plot_xas,
                      normalize)

    def XMCD(self,
             selection = [],
             remove_spikes = False,
             noise_filter = False,
             plot_xas = False,
             normalize = False):
        '''
        Function to perform XMCD Anlysis before plotting
        '''
        xvalList   = []
        yvalList   = []
        legendList = []
        infoList   = []
        newInfo    = {}
        
        pwActiveCurve = self.plotWindow.getActiveCurve()
        if pwActiveCurve:
            print len(pwActiveCurve)
            (pwActiveXVal,
             pwActiveYVal,
             pwActiveLeg,
             pwActiveInfo) = pwActiveCurve
        else:
            (pwActiveXVal,
             pwActiveYVal,
             pwActiveLeg,
             pwActiveInfo) = (None, None, None, None)
        pwX1min, pwX1max = self.plotWindow.graph.getX1AxisLimits()
        for legend in selection:
            self.plotWindow.setActiveCurve( legend )
            tmp = self.plotWindow.getActiveCurve()
            if tmp:
                (xVal,  yVal,  leg,  info) = tmp[0:4]
            else:
                # TODO: Errorhandling, curve not found
                print "Retrieved none type curve"
                continue
            if (xVal[0] < pwX1min) or (xVal[-1] > pwX1max):
                # Clip image to zoomed in Values
                mask = numpy.nonzero((xVal > pwX1min) &
                                     (xVal < pwX1max))[0]
                xVal = numpy.take(xVal,  mask)
                yVal = numpy.take(yVal,  mask)
            if remove_spikes:
                yVal = self.spikeRemoval(yVal)
            if noise_filter:
                # Savitzky Golay smoothing
                # instead of Wiener filter
                yVal = getSavitzkyGolay(yVal)
            if normalize:
                yVal = self.normalize(xVal, yVal)
            print( len(xVal), len(yVal) )
            xvalList.append(xVal)
            yvalList.append(yVal)
            legendList.append(leg)
            infoList.append(info)
            # Set info-Dictionary values
            newInfo['xlabel'] = info.get('xlabel','X')
            newInfo['ylabel'] = info.get('ylabel','Y')
            
        # Average curves
#        x, y = self.simpleMath.average(xvalList, yvalList)
        x, y = self.specAverage(xvalList, 
                                yvalList,
                                pwActiveXVal)
            
        # Naming policy
        newLegend = (self.title + '_avg')
            
        self.addCurve(x,
                      y, 
                      newLegend,
                      newInfo)
        # Reset active curve in Plotwindow 
        if pwActiveCurve:
            self.plotWindow.setActiveCurve(pwActiveLeg)
    
    def spikeRemoval(self,  inp,
                     threshold=0.001,
                     length=5):
        '''
        inp :  Numpy array 
            Contains the y-Values of the spectrum
        threshold : Float
            Threshold value for replacement
        length : Int
            Window length of medium filter. Must be
            an odd integer.
        
        Returns
        -------
        out : Numpy array
            Replaced values in the spectra whose difference
            to medium spectrum is above threshold level.
        '''
        if DEBUG:
            print 'Perform spike removal.. ', 
            print 'Threshold: %.3f, ' % threshold, 
            print 'Window length: %d' % length
        filtered = median.medfilt1d(inp, length)
        out = numpy.where(abs(filtered-inp)>threshold,
                          filtered, 
                          inp)
        return out
                           
    def normalize(self,  x, y):
        '''
        x : Numpy array 
            Contains x-Values
        y : Numpy array 
            Contains y Values
        
        Returns
        -------
        ynorm : Numpy array
            Spectrum normalized by its integral
        '''
        # Check for non-zero values?
        ynorm = y - numpy.min(y)
        ymax  = numpy.trapz(ynorm,  x)
        ynorm /= ymax
        return ynorm
        
    def specAverage(self, xarr, yarr, xrange = None):
        '''
        xarr : list
            List containing x-Values in 1-D numpy arrays
        yarr : list
            List containing y-Values in 1-D numpy arrays
        xrange : Numpy array
            x-Values used for interpolation
        
        From the spectra given in xarr & yarr, the method
        determines the overlap in the x-range. For spectra
        with unequal x-ranges, the method interpolates all
        spectra on the values given in xrange and averages
        them.
        
        Returns
        -------
        xnew, ynew : Numpy arrays
            Average spectrum
        '''
        if len(xarr) != len(yarr) and DEBUG:
            print 'specAverage: arrays of unequal length!'
        if DEBUG:
            print 'specAverage: dimension of arrays received:'
            for (i, (x, y)) in enumerate(zip(xarr, yarr)):
                print '\tarr[%d]: %d' % (i, len(y))

        same = True
        x0 = xarr[0]
        for x in xarr:
            if len(x0) == len(x):
                if numpy.all(x0 == x):
                    pass
                else:
                    same = False
                    break
            else:
                same = False
                break    

        if DEBUG:
            print 'specAverage: same = ', same

        xsort = []
        ysort = []
        for (x,y) in zip(xarr, yarr):
            if numpy.all(numpy.diff(x) > 0.):
                # All values sorted
                xsort.append(x)
                ysort.append(y)
            else:
                # Sort values
                mask = numpy.argsort(x)
                xsort.append(x.take(mask))
                ysort.append(y.take(mask))
      
        xmin0 = xsort[0][0]
        xmax0 = xsort[0][-1]
        if not same:
            # Determine global xmin0 & xmax0
            for x in xsort:
                xmin = x.min()
                xmax = x.max()
                if xmin > xmin0:
                    xmin0 = xmin
                    if DEBUG:
                        print 'New xmin0: ', xmin0
                if xmax < xmax0:
                    xmax0 = xmax
                    if DEBUG:
                        print 'New xmax0: ', xmax0
            if xmax <= xmin:
                print 'No overlap between spectra!'
                return numpy.array([]), numpy.array([])
        
        # Clip xrange to maximal overlap in spectra
        if xrange is None:
            xrange = xsort[0]
        mask = numpy.nonzero((xrange>=xmin0) & 
                             (xrange<=xmax0))[0]
        xnew = numpy.take(xrange, mask)
        ynew = numpy.zeros(len(xnew))
        
        # Perform average
        for (x, y) in zip(xsort, ysort):
            if same:
                ynew += y  
            else:
                yinter = numpy.interp(xnew, x, y)
                ynew   += numpy.asarray(yinter)
        num = len(yarr) # TODO: Cast as numpy.dtype?
        ynew /= num
        if DEBUG:
            print 'specAverage: dimension of average: ', len(xnew)
        return xnew, ynew

    
    def noiseFilter(self, y):
        size  = asarray([3] * len(y.shape))
        mean  = numpy.correlate(y, ones(size), 'same') / product(size, axis=0)
        var   = (numpy.correlate(y**2, ones(size), 'same') / product(size, axis=0)\
                - mean**2)
        noise = numpy.mean(numpy.ravel(var), axis=0)
        filt  = y - mean
        filt *= (1 - noise/var)
        filt += mean
        out   = numpy.where(var < noise, mean, res)
        
        return out
        
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
        self.plotWindow.addCurve(xVal,
                                 yVal,  
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
            self.plotWindow.addCurve(xVal,
                                     yVal, 
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
        self.plotWindow.addCurve(xVal,
                                 yVal, 
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
                self.plotWindow.addCurve(xVal,
                                         yVal,
                                         newLegend, 
                                         info,  
                                         replace=True)
            else:
                self.plotWindow.addCurve(xVal,
                                         yVal, 
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
        if self.columnCount() > 1:
            self.sortItems(1, qt.Qt.AscendingOrder)

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
                        nSelectors = 1, 
                        defaultMotor = None,
                        instructions = False,
                        selView = False,
                        cbWindow = True,
                        singleScanWindow = False):
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
        if singleScanWindow:
            self.ScanWindow = SortPlotsScanWindow(origin=plotWindow, parent=self)
        else:
            self.ScanWindow = None
        self.ScanWindowA = None
        self.ScanWindowB = None
        self.closeWidget = CloseEventNotifyingWidget.CloseEventNotifyingWidget()
        self.plotWindow = plotWindow
        self.pselection = qt.QStringListModel(self)
        self.mselection = qt.QStringListModel(self)
        self.defaultMotor = defaultMotor

        updatePixmap = qt.QPixmap(IconDict["reload"])
        buttonUpdate = qt.QPushButton(qt.QIcon(updatePixmap), 
                                      '', 
                                      self)
        
        cBoxLabel = qt.QLabel( qt.QString('Select motor:'),  self)
        for i in range(nSelectors):
            cBox = qt.QComboBox(self)
            cBox.addItems( self.motorNamesList )
            if self.defaultMotor in self.motorNamesList:
                idx = self.motorNamesList.index(self.defaultMotor)
                cBox.setCurrentIndex(idx)
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

        cBoxLayout.addWidget(cBoxLabel)
        for cBox in self.cBoxList:
            cBoxLayout.addWidget(cBox)
        cBoxLayout.addWidget(qt.HorizontalSpacer(self))
        cBoxLayout.addWidget(buttonUpdate)
        mainLayout.addWidget(self.list, 1, 0)
        if instructions:
            self.instWidget = SortPlotsInstructionWidget(self)
            mainLayout.addWidget(self.instWidget)
        else:
            self.instWidget = None    
        if cbWindow:
            print ('Creating CheckBoxWindow')
            self.checkBoxWidget = CheckBoxWindow(self)
            mainLayout.addWidget(self.checkBoxWidget)
        else:
            self.checkBoxWidget = None
        if selView:
            print 'Adding selectionView'
            self.selectionView = SPExtension(self, self.pselection,  self.mselection)
            mainLayout.addWidget(self.selectionView)
        else:
            self.selectionView = None
        
#        self.resize(500,  500)

        buttonUpdate.clicked.connect(self.updatePlots)
        
        if singleScanWindow:
            # Add ScanWindow to the Layout
            nRowsLayout = mainLayout.rowCount()
            mainLayout.addWidget( self.ScanWindow, 0, 1, nRowsLayout, 1 )

        self.setLayout(mainLayout)

        self.updateTree()
        self.updateActionList(
              [('Invert selection', self.list.invertSelection), 
                ('Process as A',    self.processAsA), 
                ('Process as B',    self.processAsB), 
                ('Remove curve(s)', self.removeCurve_)])

    def getAllMotorNames(self):
        names = []
        for dic in self.motorsList:
            for key in dic.keys():
                if key not in names:
                    names.append( key )
        names.sort()
        return names

    def removeCurve_(self):
        sel = self.list.selectedItems()
        self.plotWindow.removeCurves(sel)
        self.updatePlots()

    def updatePlots(self,
                    newLegends = None,  
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
            act = qt.QAction(fName, self)
            act.triggered.connect(function)
            self.actionList += [act]
            if self.instWidget:
                self.instWidget.actionCBox.addItem(fName)
    
    def processAs(self,  pol):
        if not self.selectionView:
            return
        sel = self.list.selectedItems()
        if pol == 'p':
            self.pselection.setStringList(sel)
        elif pol == 'm':
            self.mselection.setStringList(sel)

    def processAsA(self):
        if self.ScanWindowA:
            if self.ScanWindowA.isVisible():
                self.ScanWindowA.close()
        self.ScanWindowA = self.spawnNewScanWindow('A')

        if self.selectionView:
            sel = self.list.selectedItems()
            self.pselection.setStringList(sel)
        
        print 'ProcessAsA finished'
    
    def processAsB(self):
        if self.ScanWindowB:
            if self.ScanWindowB.isVisible():
                self.ScanWindowB.close()
        self.ScanWindowB = self.spawnNewScanWindow('B')
        print 'ProcessAsB finished'

    def spawnNewScanWindow(self,  name=''):
        removeSpikes = False
        noiseFilter  = False
        plotXAS      = True
        normalize    = False
        if self.checkBoxWidget:
            (removeSpikes, 
             noiseFilter, 
             plotXAS, 
             normalize) = self.checkBoxWidget.check()
        if DEBUG:
            print 'removeSpikes = ',  removeSpikes
            print 'noiseFilter = ',  noiseFilter
            print 'plotXAS = ', plotXAS
            print 'normalize = ', normalize
        sel = self.list.selectedItems()
        swin  = SortPlotsScanWindow(self.plotWindow, 
                                    name, 
                                    sel,
                                    None,
                                    remove_spikes = removeSpikes,
                                    noise_filter = noiseFilter,
                                    plot_xas = plotXAS,
                                    normalize = normalize)
        swin.plotModifiedSignal.connect(self.updatePlots)
        if name == 'A':
            self.pselection.setStringList(sel)
        elif name == 'B':
            self.mselection.setStringList(sel)
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
    y0 =  10 * x + 10000. * numpy.exp(-0.5*(x-500)*(x-500)/400) + 1500 * numpy.random.random(1000.)
    y1 =  10 * x + 10000. * numpy.exp(-0.5*(x-600)*(x-600)/400) + 1500 * numpy.random.random(1000.)
#    y1 =  10 *x + 1000 
#    y2 =  10 * x + 10000. * numpy.exp(-0.5*(x-400)*(x-400)/400)
    y2 =  10 * x + 10000. * numpy.exp(-0.5*(x-400)*(x-400)/400) + 1500 * numpy.random.random(1000.)
    y2[320:322] = 50000.
    swin.addCurve(x, y2, legend="Curve2", replot=False, replace=False)
    swin.addCurve(x, y0, legend="Curve0", replot=False, replace=False)
    swin.addCurve(x, y1, legend="Curve1", replot=False, replace=False)

    w = SortPlotsWidget(None, swin.getAllCurves(just_legend=True),  motors,  swin,  selView = False, defaultMotor = 'Motor11')
#    w = CheckBoxWindow(None)
    w.show()
    app.exec_()
    
if __name__ == '__main__':
    main()
