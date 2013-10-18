import numpy, copy
from PyMca import PyMcaDirs
from os.path import splitext
from os import linesep as newline
from PyMca import PyMcaQt as qt
from PyMca.PyMca_Icons import IconDict

from PyMca import CloseEventNotifyingWidget
from PyMca import ScanWindow as sw

DEBUG = True
if DEBUG:
    numpy.set_printoptions(threshold=50)

class SortPlotsScanWindow(sw.ScanWindow):

    plotModifiedSignal = qt.pyqtSignal()

    def __init__(self,
                 origin,
                 parent=None):
        sw.ScanWindow.__init__(self, 
                               parent, 
                               name='XMCD Analysis', 
                               specfit=None)
        self.plotWindow = origin

        self.scanWindowInfoWidget.hide()

        # Buttons to push spectra to main Window
        buttonShowXAS = qt.QPushButton('Hide XAS', self)
        buttonShowXMCD = qt.QPushButton('Hide XMCD', self)
        buttonAdd = qt.QPushButton('Add',  self)
        buttonReplace = qt.QPushButton('Replace',  self)
        buttonAddAll = qt.QPushButton('Add all',  self)
        buttonReplaceAll = qt.QPushButton('Replace all',  self)
        buttonLayout = qt.QHBoxLayout(None)
        buttonLayout.setContentsMargins(0, 0, 0, 0)
        buttonLayout.setSpacing(5)
        # Show XAS & XMCD Buttons
        buttonLayout.addWidget(buttonShowXAS)
        buttonLayout.addWidget(buttonShowXMCD)
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

        # Copy spectra from origin
        self.selectionDict = {'m':[], 'p':[]}
        self.curvesDict = {}
        self.xmcdCheck = (None, None, None, None)

        self.xRange = None
        # Keep track of Averages, XMCD and XAS curves
        self.avgM = None
        self.avgP = None
        self.xmcd = None
        self.xas  = None

    def interpXRange(self, equidistant=True, xRangeList=None):
        '''
        Input
        -----
        fromCurves : Bool
            Uses curves present in self.curvesDict
            if set to true. If set to false, an
            ndarray with equistant interpolation
            points is returned
        xRangeList : List
            List of ndarray from whose overlap
            is determined
        
        Checks dataObjectsDictionary for curves
        present in the ScanWindow and tries to
        find the overlap in the x-range of these
        curves.

        If equidistant is True:
        The x-ranges of the curves containing
        the minmal resp. the maximal x-value are
        used to determine the minimal number of 
        points (n) in their overlap in order to
        avoid oversampling.
        
        If equidistant is False:
        If an active curve in the plotWindow is
        set, its xRange is used as interpolation 
        points. If not active curve is set, the
        first curve in the scan sequence is taken
        
        Returns
        -------
        out : numpy array
            (Evenly spaced) x-range between xmin and
            xmax containing n points
        '''
        if not xRangeList:
            xRangeList = [v.x[0] for v in self.curvesDict.values()]
        if not len(xRangeList):
            # Nothing to do..
            print 'Nothing to do'
            return None
            
        num = 0
        xmin, xmax = self.plotWindow.getGraphXLimits()
        for x in xRangeList:
            if x.min() > xmin:
                xmin = x.min()
                if DEBUG:
                    print 'New minimum: ', xmin
            if x.max() < xmax:
                xmax = x.max()
                if DEBUG:
                    print 'New maximum: ', xmax
        
        if xmin >= xmax:
            # TODO: Somekind of Error..
            pass
        if equidistant:
            for x in xRangeList:
                curr = numpy.nonzero((x >= xmin) & 
                                     (x <= xmax))[0].size
                num = curr if curr>num else num
            # Exclude first and last point
            out = numpy.linspace(xmin, xmax, num, endpoint=False)[1:]
        else:
            active = self.plotWindow.graph.getActiveCurve(just_legend=True)
            if active:
                curve = self.plotWindow.dataObjectsDict[active]
                if DEBUG:
                    print 'interpXRange -- Active curve is \'%s\''%active
            else:
                first = sorted(self.curvesDict.keys())[0]
                if DEBUG:
                    print 'interpXRange -- No active curve,',
                    print 'proceed with \'%s\' as first'%first
                curve = self.plotWindow.dataObjectsDict[first]
            x = curve.x[0]
            mask = numpy.nonzero((x >= xmin) &
                                 (x <= xmax))[0]
            out = numpy.sort(numpy.take(x, mask))
        if DEBUG:
            print 'interpXRange -- Resulting xrange:'
            print '\tmin = ', out.min()
            print '\tmax = ', out.max()
            print '\tnum = ', len(out)
        return out

    def processSelection(self, msel, psel):
        all = self.getAllCurves(just_legend=True)
        self.removeCurves(all)
        self.avgP, self.avgM = None, None
        
        self.selectionDict['m'] = msel[:]
        self.selectionDict['p'] = psel[:]
        if DEBUG:
            print 'm: ', self.selectionDict['m']
            print 'p: ', self.selectionDict['p']
        self.curvesDict = self.copyCurves(msel + psel)
        
        if (len(self.curvesDict) == 0) or\
           ((len(self.selectionDict['m']) == 0) and\
           (len(self.selectionDict['p']) == 0)):
            # Nothing to do
            return

        xRange = self.interpXRange(equidistant=True) # TODO: Make equidistant optional
        activeLegend = self.plotWindow.graph.getActiveCurve(justlegend=True)
        if (not activeLegend) or (activeLegend not in self.curvesDict.keys()):
            # Use first curve in the series as xrange
            activeLegend = sorted(self.curvesDict.keys())[0]
        active = self.curvesDict[activeLegend]
        xlabel, ylabel = self.extractLabels(active.info)
        
        # Calculate averages and add them to the plot
        for id in ['p','m']:
            sel = self.selectionDict[id]
            if not len(sel):
                continue
            xvalList = []
            yvalList = []
            for legend in sel:
                tmp = self.curvesDict[legend]
                xvalList.append(tmp.x[0])
                yvalList.append(tmp.y[0])
            avg_x, avg_y = self.specAverage(xvalList,
                                            yvalList,
                                            xRange)
            avgName = 'avg_' + id
#            print 'id  :', id
#            print 'xrng:', xRange
#            print 'xVal:', xvalList
#            print 'yVal:', yvalList
#            print 'avgx:', avg_x
#            print 'avgy:', avg_y
            self.newCurve(avg_x,
                          avg_y,
                          avgName,
                          xlabel,
                          ylabel)
            if id == 'p':
                self.avgP = self.dataObjectsList[-1]
            if id == 'm':
                self.avgM = self.dataObjectsList[-1]
            
        if (self.avgM and self.avgP):
            # TODO: if self.showXMCD() ...
            self.performXMCD()
            # TODO: if self.showXAS() ...
            self.performXAS()

    def copyCurves(self, selection):
        '''
        selection : List
            Contains names of curves to be processed
    
        Returns
        -------
        out : Dictionary
            Contains legends as keys and dataObjects
            as values.
        '''
        if not len(selection):
            return {}
        out = {}
        for legend in selection:
            tmp = self.plotWindow.dataObjectsDict.get(legend, None)
            if tmp:
                out[legend] = copy.deepcopy(tmp)
            else:
                # TODO: Errorhandling, curve not found
                if DEBUG:
                    print "copyCurves -- Retrieved none type curve"
                continue
        if DEBUG:
            for c in out.values():
                x = c.x[0]
                print '\tmin:',x.min(),'\tmax:',x.max()
        return out

    def specAverage(self, xarr, yarr, xrange=None):
        '''
        xarr : list
            List containing x-Values in 1-D numpy arrays
        yarr : list
            List containing y-Values in 1-D numpy arrays
        xrange : Numpy array
            x-Values used for interpolation. Must overlap
            with all arrays in xarr

        From the spectra given in xarr & yarr, the method
        determines the overlap in the x-range. For spectra
        with unequal x-ranges, the method interpolates all
        spectra on the values given in xrange and averages
        them.

        Returns
        -------
        xnew, ynew : Numpy arrays or None
            Average spectrum. In case of invalid input,
            (None, None) tuple is returned.
        '''
        if (len(xarr) != len(yarr)) or\
           (len(xarr) == 0) or (len(yarr) == 0):
            if DEBUG:
                print 'specAverage -- invalid input!'
            return None, None 

        same = True
        if xrange == None:
            x0 = xarr[0]
        else:
            x0 = xrange
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
            print 'specAverage -- same = ', same

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

        if xrange != None:
            xmin0 = xrange.min()
            xmax0 = xrange.max()
        else:
            xmin0 = xsort[0][0]
            xmax0 = xsort[0][-1]
        if (not same) or (xrange == None):
            # Determine global xmin0 & xmax0
            for x in xsort:
                xmin = x.min()
                xmax = x.max()
                if xmin > xmin0:
                    xmin0 = xmin
                    if DEBUG:
                        print 'specAverage -- New xmin0: ', xmin0
                if xmax < xmax0:
                    xmax0 = xmax
                    if DEBUG:
                        print 'specAverage -- New xmax0: ', xmax0
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
            print 'specAverage -- xrange: '
            print '\tmin = ', xnew.min()
            print '\tmax = ', xnew.max()
            print '\tnum = ', len(xnew)
        return xnew, ynew

    def extractLabels(self, info):
        xlabel = 'X'
        ylabel = 'Y'
        sel = info.get('selection', None)
        labelNames = info.get('LabelNames',[])
        if sel:
            xsel = sel.get('x',[])
            ysel = sel.get('y',[])
            if len(xsel) > 0:
                x = xsel[0]
            else:
                y = -1
            if len(ysel) > 0:
                y = ysel[0]
            else:
                y = -1
            if len(labelNames) == 2:
                [xlabel, ylabel] = labelNames
            elif (len(labelNames) > max(x,y)):
                if y > 0:
                    ylabel = labelNames[y]
                if x > 0:
                    xlabel = labelNames[x]
        return xlabel, ylabel


    def performXAS(self):
        keys = self.dataObjectsDict.keys()
        if (self.avgM in keys) and (self.avgP in keys):
            m = self.dataObjectsDict[self.avgM]
            p = self.dataObjectsDict[self.avgP]
        else:
            if DEBUG:
                print 'performXAS -- Data not found: '
                print '\tavg_m = ', self.avgM
                print '\tavg_p = ', self.avgP
            return
        if numpy.all( m.x[0] == p.x[0] ):
            avg = .5*(p.y[0] + m.y[0])
        else:
            if DEBUG:
                print 'performXAS -- x ranges are not the same! ',
                print 'Force interpolation'
            avg = self.performAverage([m.x[0], p.x[0]],
                                      [m.y[0], p.y[0]],
                                       p.x[0])
        xmcdLegend = 'XAS'
        xlabel, ylabel = self.extractLabels(m.info)
        self.newCurve(m.x[0],
                      avg, 
                      xmcdLegend,
                      xlabel,
                      ylabel)
        self.xas = self.dataObjectsList[-1]

    def performXMCD(self):
        if DEBUG:
            print 'performXMCD -- xmcdCheck: ',
            print "avgP = '%s', "%self.avgP,
            print "avgM = '%s'"%self.avgM
        keys = self.dataObjectsDict.keys()
        if (self.avgM in keys) and (self.avgP in keys):
            m = self.dataObjectsDict[self.avgM]
            p = self.dataObjectsDict[self.avgP]
        else:
            if DEBUG:
                print 'performXMCD -- Data not found:'
            return
        if numpy.all( m.x[0] == p.x[0] ):
            diff = p.y[0] - m.y[0]
        else:
            if DEBUG:
                print 'performXMCD -- x ranges are not the same! ',
                print 'Force interpolation using p Average xrange'
            # Use performAverage d = 2 * avg(y1, -y2)
            # and force interpolation on p-xrange
            diff = 2. * self.performAverage([m.x[0], p.x[0]],
                                            [-m.y[0], p.y[0]],
                                            p.x[0])
        xmcdLegend = 'XMCD'
        xlabel, ylabel = self.extractLabels(m.info)
        # TODO: replace newCurve with addCurve -> use Plugin inteface
        self.newCurve(p.x[0],
                      diff, 
                      xmcdLegend,
                      xlabel,
                      ylabel)
        self.graph.mapToY2(' '.join([xmcdLegend, ylabel]))
        self._zoomReset()
        self.xmcd = self.dataObjectsList[-1]

    def _saveIconSignal(self):
        saveDir = PyMcaDirs.outputDir
        filter = 'CSV Files (*.csv);;Any File (*.*)'
        # TODO: QFileDialog.getSaveFileName returns
        # filename w/o extension!
        filename = qt.QFileDialog.getSaveFileName(self,
                                                'Save XMCD Analysis',
                                                saveDir,
                                                filter)
        filename = str(filename)
#        if not len(filename):
#            raise IOError('Invalid filename')
        if not len(splitext(filename)[1]):
            # Add '.csv' or some extension
            pass
        # Get all values shown. Possible ways:
        # - Loop through dataObjectsDict
        # - Access avgM, avgP, XMCD and XAS directly
        try:
            filehandle = open(filename, 'w')
        except IOError:
            msg = qt.QMessageBox(text="Unable to write to '%s'"%filename)
            return
        curves = self.dataObjectsDict.values()
        yVals = [curve.y[0] for curve in curves]
        xVals = [curves[0].x[0]]
        outArray = numpy.vstack([xVals, yVals]).T
        # Optional: Use numpy.savetxt
        for line in outArray:
            tmp = ','.join([str(num) for num in line])
            filehandle.write(tmp + newline)
        filehandle.close()
    
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

    # TODO: NAMING in add ... replaceAll
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

    def detrend(self):
        for (k,v) in self.dataObjectsDict.items():
            if k.startswith('XMCD'):
                xmcd = v
                xmcdLegend = k
                break
        xmcd = self.dataObjectsDict[xmcdLegend]
        x = xmcd.x[0]
        y = xmcd.y[0]
        a, b = numpy.polyfit(x,y,1)
        ynew = y - a*x - b
        y = ynew
        self.graph.checky2scale()
        

    
class SortPlotsMenu(qt.QMenu):
    def __init__(self,  parent, title=None):
        qt.QMenu.__init__(self,  parent)
        if title:
            self.setTitle(title)
        
    def setActionList(self, actionList, update=False):
        '''
        List functions has to have the form (functionName, function)

        Default is ('', function)
        '''
        if not update:
            self.clear()
        for (name, function) in actionList:
            if name == '$SEPERATOR':
                self.addSeparator()
                continue
            if name != '':
                fName = name
            else:
                fName = function.func_name
            act = qt.QAction(fName,  self)
            # Force triggered() instead of triggered(bool)
            # to ensure proper interaction with default parameters
            act.triggered[()].connect(function)
            self.addAction(act)

class SortPlotsTreeWidget(qt.QTreeWidget):

    selectionModifiedSignal = qt.pyqtSignal()

    def __init__(self,  parent, identifiers = ['p','m','d']):
        qt.QTreeWidget.__init__(self,  parent)
        self.identifiers = identifiers
        self.actionList  = []
        self.contextMenu = qt.QMenu('Perform',  self)

    def sizeHint(self):
        vscrollbar = self.verticalScrollBar()
        width = vscrollbar.width()
        for i in range(self.columnCount()):
            width += (2 + self.columnWidth(i))
        return qt.QSize( width, 200 )

    def setContextMenu(self, menu):
        self.contextMenu = menu

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
            self.contextMenu.popup(pos)
        event.accept()

    def invertSelection(self):
        root = self.invisibleRootItem()
        for i in range(root.childCount()):
            if root.child(i).isSelected():
                root.child(i).setSelected(False)
            else:
                root.child(i).setSelected(True)

    def getColumn(self, ncol, selectionOnly=False, convertType=str):
        '''
        Returns items in tree column ncol and converts them
        to convertType. If the conversion fails, the default
        type is a python string.
        
        If selectionOnly is set to True, only the selected
        the items of selected rows are returned.
        '''
        out = []
        convert = (convertType != str)
        if ncol > (self.columnCount()-1):
            if DEBUG:
                print 'getColum -- Selected column out of bounds'
            raise IndexError("Selected column '%d' out of bounds" % ncol)
            return out
        if selectionOnly:
            sel = self.selectedItems()
        else:
            root = self.invisibleRootItem()
            sel = [root.child(i) for i in range(root.childCount())]
        for item in sel:
            # TODO: To raise or not to raise?
            tmp = str(item.text(ncol))
            if convert:
                try:
                    tmp = convertType(tmp)
                except (TypeError, ValueError):
                    if DEBUG:
                        print 'getColum -- Conversion failed!'
                    if convertType == float:
                        tmp = float('NaN')
                    else:
                        raise TypeError        
            out += [tmp]
        return out

    def build(self,  items,  headerLabels):
        '''
        (Re-) Builds the tree display

        headerLabels must be of type QStringList
        items must be of type [QStringList]
        '''
        # Remember selection, then clear list
        sel = self.getColumn(1, True)
        self.clear()
        self.setHeaderLabels(headerLabels)
        for item in items:
            treeItem = qt.QTreeWidgetItem(self,  item)
            if treeItem.text(1) in sel:
                treeItem.setSelected(True)
        self.resizeColumnToContents(0)
        self.resizeColumnToContents(1)

    def setSelectionAs(self, id):
        '''
        Sets the items currently selected to 
        the identifier given in id.
        '''
        if id not in self.identifiers:
            # TODO: Error! Raise exception or something..
            id = ''
        sel = self.selectedItems()
        if id == self.identifiers[-1]:
            id = ''
        for item in sel:
            item.setText(0,id)
        self.selectionModifiedSignal.emit()

    def setSelectionToSequence(self, seq=None, selectionOnly=False):
        '''
        Sets the id column (col 0) to seq. If
        sequence is None, a dialog window is 
        shown.
        '''
        chk = True
        if selectionOnly:
            sel = self.selectedItems()
        else:
            root = self.invisibleRootItem()
            sel = [root.child(i) for i in range(root.childCount())]
        # Ensure alphabetically ordered list
        self.sortItems(1, qt.Qt.AscendingOrder)
        if not seq:
            seq, chk = qt.QInputDialog.\
                getText(None, 
                        'Sequence Dialog', 
                        'Valid identifiers are: ' + ', '.join(self.identifiers),
                        qt.QLineEdit.Normal, 
                        'Enter sequence')
        seq = str(seq)
        if not chk:
            return
        for id in seq:
            if id not in self.identifiers:
                invalidMsg = qt.QMessageBox(None)
                invalidMsg.setText('Invalid sequence. Try again.')
                invalidMsg.setStandardButtons(qt.QMessageBox.Ok)
                invalidMsg.exec_()
                return
        if len(sel) != len(seq):
            invalidMsg = qt.QMessageBox(None)
            invalidMsg.setText('Sequence length does not match item count.')
            invalidMsg.setStandardButtons(qt.QMessageBox.Ok)
            invalidMsg.exec_()
            return
        for (id, item) in zip(seq, sel):
            if id == self.identifiers[-1]:
                id = ''
            item.setText(0, id)
        self.selectionModifiedSignal.emit()

    def clearSelection(self, selectionOnly=True):
        '''
        Empties the id column for the selected rows.
        '''
        if selectionOnly:
            sel = self.selectedItems()
        else:
            root = self.invisibleRootItem()
            sel = [root.child(i) for i in range(root.childCount())]
        for item in sel:
            item.setText(0,'')
        self.selectionModifiedSignal.emit()

    def getSelection(self):
        '''
        Returns dictionary with where the keys
        are the identifiers and the values are
        (sorted) lists containing legends to 
        which the respective identifier is
        assigned to.
        '''
        out = dict((id, []) for id in self.identifiers)
        root = self.invisibleRootItem()
        for i in range(root.childCount()):
            item   = root.child(i)
            id     = str(item.text(0))
            legend = str(item.text(1))
            if len(id) == 0:
                id = self.identifiers[-1]
            out[id] += [legend]
        for value in out.values():
            value.sort()
        return out


class SortPlotsWidget(CloseEventNotifyingWidget.CloseEventNotifyingWidget):

    setSelectionSignal = qt.pyqtSignal(object, object)

    def __init__(self,  parent,
                        plotWindow,
                        beamline,
                        nSelectors = 2):
        """
        Input
        -----
        plotWindow : ScanWindow instance
            ScanWindow from which curves are passed for
            XMCD Analysis
        nSelectors : Int
            Number of Comboboxes shown in the widget.
            Every Combobox allows to select a different motor
        
        legendList : List 
            Contains curve legends. Format
            ['Legend0', ... , 'LegendX']
        motorsList : List 
            Contains dictionaries. Format:
            [{'MotorName0':MotorValue0, ... , 'MotorNameN':MotorValueN},
             ...,
             {'MotorName2':MotorValue2, ... , 'MotorNameM':MotorValueM}]
        """
        CloseEventNotifyingWidget.\
            CloseEventNotifyingWidget.__init__(self,  parent)
        self.plotWindow = plotWindow
        self.legendList = []
        self.motorsList = []
        # Set self.plotWindow before calling self._setLists!
        self._setLists()
        self.motorNamesList = [''] + self._getAllMotorNames()
        self.motorNamesList.sort()
        self.numCurves = len(self.legendList)
        self.cBoxList = []
        self.ScanWindow = SortPlotsScanWindow(origin=plotWindow, 
                                              parent=None)
                                              
        self.notifyCloseEventToWidget(self)
        self.notifyCloseEventToWidget(self.ScanWindow)
        self.selectionDict = {'d': [],
                              'p': [],
                              'm': []}
        self.beamline = beamline
        self.setSizePolicy(qt.QSizePolicy.MinimumExpanding, 
                           qt.QSizePolicy.Expanding)

        self.setWindowTitle("Sort Plots Window")
        updatePixmap = qt.QPixmap(IconDict["reload"])
        buttonUpdate = qt.QPushButton(qt.QIcon(updatePixmap), '', self)
        for i in range(nSelectors):
            cBox = qt.QComboBox(self)
            cBox.addItems(self.motorNamesList)
            cBox.activated['QString'].connect(self.updateTree)
            self.cBoxList += [cBox]

        self.list = SortPlotsTreeWidget(self)
        labels = ['Legend'] + nSelectors*['']
        ncols  = len(labels)
        self.list.setColumnCount(ncols)
        self.list.setHeaderLabels(labels)
        self.list.setSortingEnabled(True)
        self.list.setSelectionMode(qt.QAbstractItemView.ExtendedSelection)
        listContextMenu = SortPlotsMenu(None)
        listContextMenu.setActionList(
              [('XMCD analysis', self.triggerXMCD),
               ('$SEPERATOR', None),
               ('Set as p', self.setAsP),
               ('Set as m', self.setAsM),
               ('Enter sequence', self.list.setSelectionToSequence),
               ('Remove selection', self.list.clearSelection),
               ('$SEPERATOR', None),
               ('Invert selection', self.list.invertSelection), 
               ('Remove curve(s)', self.removeCurve_)])
        self.list.setContextMenu(listContextMenu)
        
        self.splitter = qt.QSplitter(qt.Qt.Horizontal, self)
        mainLayout = qt.QGridLayout(None)
        cBoxLayout = qt.QHBoxLayout(None)
        topLayout  = qt.QHBoxLayout(None)
        mainLayout.setContentsMargins(1, 1, 1, 1)
        mainLayout.setSpacing(2)
        mainLayout.addLayout(cBoxLayout, 2, 0)
        mainLayout.addWidget(self.list, 1, 0)
        mainLayout.addLayout(topLayout, 0, 0)
        leftWidget = qt.QWidget(self)
        leftWidget.setLayout(mainLayout)
        self.splitter.addWidget(leftWidget)
        self.splitter.addWidget(self.ScanWindow)
        self.splitter.setSizes([leftWidget.minimumSize(),self.list.size()])
        
        self.autoselectCBox = qt.QComboBox(self)
        self.autoselectCBox.addItems(
                        ['Select Experiment',
                         'ID08: Linear Dichorism',
                         'ID08: XMCD'])
        self.autoselectCBox.activated['QString'].connect(self.autoselect)
        
        topLayout.addWidget(buttonUpdate)
        topLayout.addWidget(qt.HorizontalSpacer(self))
        topLayout.addWidget(self.autoselectCBox)
        
        cBoxLayout.addWidget(qt.HorizontalSpacer(self))
        cBoxLayout.addWidget(
                qt.QLabel('Selected motor(s):',  self))
        for cBox in self.cBoxList:
            cBoxLayout.addWidget(cBox)   

        buttonUpdate.clicked.connect(self.updatePlots)
        self.list.selectionModifiedSignal.connect(self.updateSelectionDict)
        self.setSelectionSignal.connect(self.ScanWindow.processSelection)

#        self.setLayout(mainLayout)
        self.updateTree()
        self.list.sortByColumn(1, qt.Qt.AscendingOrder)
        self._setBeamlineSpecific(self.beamline)
    
    def autoselect(self, option):
        '''
        Beamline-specific experiments rely on specific
        motors. In order to implement new experiments:
        
        1.  Pick a name: The first characters must contain
            the beamline name (ex: 'ID99: Foobar Experiment')
        2.  Set related motor names: variables motor0 and
            motor1 must be given accordingly to the motors
            the experiment depends on.
        3.  Retrieve motor settings: Returned in numpy arrays
            values0 and values1. Non-float values are masked
            by NaN-values.
        4.  Set values: The curves are assigned into a
            p/m-selection depeding on the numpy-array values
            and a pivot element. values needs to be calculated
            from value0 and value1, vpivot needs to be set.
        5.  Add experiment to the items of self.autoselectCBox
        '''
        option = str(option)
        cBox0  = self.cBoxList[0]
        cBox1  = self.cBoxList[1]
        # Set motor names
        if option == 'ID08: Linear Dichorism':
            motor0 = 'PhaseD'
            motor1 = ''
        elif option == 'ID08: XMCD':
            motor0 = 'PhaseD'
            motor1 = 'oxPS'
        else:
            motor0 = ''
            motor1 = ''
        # Check if motor names are valid and retrieve motor settings
        try:
            self.setComboBoxToMotor(0, motor0)
            self.setComboBoxToMotor(1, motor1)
        except ValueError as e:
            if DEBUG:
                return 'autoselect --', e
            return
        values0 = numpy.array(self.list.getColumn(2, convertType=float))
        values1 = numpy.array(self.list.getColumn(3, convertType=float))
        if option == 'ID08: Linear Dichorism':
            values = values0
            vmin = values.min()
            vmax = values.max()
            vpivot = .5 * (vmax + vmin)
        elif option == 'ID08: XMCD':
            values = values0 * values1
            vpivot = 0.
        else:
            values = numpy.array([float('NaN')]*len(self.legendList))
            vpivot = 0.
        seq = ''
        for x in values:
            if str(x) == 'nan':
                seq += 'd'
            elif x>vpivot:
                seq += 'p'
            else:
                seq += 'm'
        self.list.setSelectionToSequence(seq)
        
    def setComboBoxToMotor(self, idx, motor):
        '''
        idx : int
            Index of the combobox to be selected.
            Must be between [0, nSelectors-1]
        motor : python str
            Must be in motorNamesList
            
        Sets the idx-th combobox in self.cBoxList
        to motor and updates the tree.
        '''
        max = len(self.cBoxList)
        if not (idx < max):
            raise IndexError('Only %d Comboboxes present'%len(max))
        cBox  = self.cBoxList[idx]
        motors = [str(cBox.itemText(i)) for i in range(cBox.count())]
        if motor in motors:
            index = motors.index(motor)
            cBox.setCurrentIndex(index)
            self.updateTree()
        else:
            raise ValueError('"%s" not in motors'%motor)

    def triggerXMCD(self):
        if self.ScanWindow:
            # TODO: Delete me
            self.ScanWindow.show()
            self.ScanWindow.raise_()
        msel = self.selectionDict['m']
        psel = self.selectionDict['p']
#        self.ScanWindow.performXMCD()
        self.ScanWindow.processSelection(msel, psel)

    def removeCurve_(self):
        sel = self.list.getColumn(1, True, str)
        if DEBUG:
            print 'removeCurve_ -- selection:\n\t',
            print '\n\t'.join(sel)
        self.plotWindow.removeCurves(sel)
        # TODO: also remove in ScanWindow
        self.updatePlots()

    def updateSelectionDict(self):
        self.selectionDict = self.list.getSelection()
        if DEBUG:
            print 'updateSelectionDict -- selectionDict: ',
            print self.selectionDict
        self.setSelectionSignal.emit(self.selectionDict['m'],
                                     self.selectionDict['p'])        

    def updatePlots(self,
                    newLegends = None,  
                    newMotorValues = None):
        self._setLists()
        self.motorNamesList = [''] + self._getAllMotorNames()
        self.motorNamesList.sort()
        for cBox in self.cBoxList:
            motorName = cBox.currentText()
            index = cBox.findText(motorName)
            # index = -1, if motorName not found
            index = 0 if index<0 else index
            cBox.clear()
            cBox.addItems(self.motorNamesList)
            cBox.setCurrentIndex(index)
        self.updateTree()
        # TODO: Determine x ranges in plotWindow
#        if plotWindow.isZoomed
        
        # in order to check if they have changed
        experiment = str(self.autoselectCBox.currentText())
        if experiment != 'Select Experiment':
            # Unable to open 
            self.autoselect(experiment)
        

    def updateTree(self):
        mList  = [ str(cBox.currentText()) for cBox in self.cBoxList ]
        labels = ['#','Legend'] + mList
        items  = []
        for i in range(len(self.legendList)):
            legend = self.legendList[i]
            values = self.motorsList[i]
            selection = ''
            for (id,v) in self.selectionDict.items():
                if legend in v:
                    if id != 'd':
                        selection = id
                    break
            tmp = qt.QStringList([selection, legend])
            for m in mList:
                if m == '':
                    tmp.append('')
                else:
                    tmp.append(str(values.get(m, '---')))
            items.append(tmp)
        self.list.build(items,  labels)

    def setAsM(self):
        self.list.setSelectionAs('m')

    def setAsP(self):
        self.list.setSelectionAs('p')

    def _getAllMotorNames(self):
        names = []
        for dic in self.motorsList:
            for key in dic.keys():
                if key not in names:
                    names.append(key)
        names.sort()
        return names

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
                ret.append(dict(zip(namesList,  valuesList)))
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
        self.motorsList = self._convertInfoDictionary(infoList)

    def _setBeamlineSpecific(self, beamline):
        '''
        beamline : python str
            Beamline identifier, all upper case
            (e.g. ID08, ID12, ..)
        '''
        options = [str(self.autoselectCBox.itemText(i))\
                   for i in range(self.autoselectCBox.count())]
        for (i, option) in enumerate(options):
            if option.startswith(beamline):
                self.autoselectCBox.setCurrentIndex(i)
                self.autoselectCBox.activated[qt.QString].emit(option)
                break

def main():
    import sys,  numpy
    app = qt.QApplication(sys.argv)
    swin = sw.ScanWindow()
    info0 = {'xlabel': 'foo', 'ylabel': 'arb', 'MotorNames': 'oxPS Motor11 Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3',  'MotorValues': '1 8.69271399699 21.9836418539 0.198068826612 0.484475455792 0.350252217264 0.663925270933 0.813033264421 0.221149410218 0.593188258866 0.678010392881 0.267389247833 0.677890617858'}
    info1 = {'MotorNames': 'PhaseD oxPS Motor16 Motor15 Motor14 Motor13 Motor12 Motor11 Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3', 'MotorValues': '0.470746882688 0.695816070299 0.825780811755 0.25876374531 0.739264467436 0.090842892619 2 0.213445659833 0.823400550314 0.020278096857 0.568744021322 0.85378115537 0.696730386891 0.269196313956 0.793293334395 0.769216567757 0.959092709527 0.0109264683697 0.538264972553'}
    info2 = {'MotorNames': 'PhaseD oxPS Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3',  'MotorValues': '2 0.44400576644 0.613870067852 0.901968648111 0.319768771085 0.571432278628 0.278675836163 0.154436774878 0.416231999332 0.294201017231 0.813913587748 0.577572903105 0.869045182568'}
    x = numpy.arange(100.,1100.)
    y0 =  10 * x + 10000. * numpy.exp(-0.5*(x-500)*(x-500)/400) + 1500 * numpy.random.random(1000.)
    y1 =  10 * x + 10000. * numpy.exp(-0.5*(x-600)*(x-600)/400) + 1500 * numpy.random.random(1000.)
    y2 =  10 * x + 10000. * numpy.exp(-0.5*(x-400)*(x-400)/400) + 1500 * numpy.random.random(1000.)
    y2[320:322] = 50000.
    swin.newCurve(x, y2, legend="Curve2", xlabel='ene_st2', ylabel='zratio2', info=info2, replot=False, replace=False)
    swin.newCurve(x, y0, legend="Curve0", xlabel='ene_st0', ylabel='zratio0', info=info0, replot=False, replace=False)
    swin.newCurve(x, y1, legend="Curve1", xlabel='ene_st1', ylabel='zratio1', info=info1, replot=False, replace=False)
    
    w = SortPlotsWidget(None, swin, 'ID08', nSelectors = 2)
    w.show()
    app.exec_()

if __name__ == '__main__':
    main()
