import numpy, copy
from PyMca import PyMcaDirs, PyMcaFileDialogs
from PyMca import ConfigDict
from os.path import splitext
from os import linesep as NEWLINE
from PyMca import PyMcaQt as qt
from PyMca.PyMca_Icons import IconDict

from PyMca import ScanWindow as sw

DEBUG = True
if DEBUG:
    numpy.set_printoptions(threshold=50)

class TreeWidgetItem(qt.QTreeWidgetItem):
    def __init__(self, parent, itemList):
        qt.QTreeWidgetItem.__init__(self, parent, itemList)
        
    def __lt__(self, other):
        col = self.treeWidget().sortColumn()
        val      = self.text(col)
        valOther = other.text(col)
        if val == '---':
                ret = True
        elif col > 0:
            try:
                ret  = (float(val) < float(valOther))
            except ValueError:
                ret  = qt.QTreeWidgetItem.__lt__(self, other)
        else:
            ret  = qt.QTreeWidgetItem.__lt__(self, other)
        return ret

class XMCDOptions(qt.QDialog):
    
    def __init__(self, parent, mList, full=True):
        qt.QDialog.__init__(self, parent)
        self.setWindowTitle('XLD/XMCD Options')
        self.setModal(True)
        self.motorList = mList
        self.saved = False

        # Buttons
        buttonOK = qt.QPushButton('OK')
        buttonCancel = qt.QPushButton('Cancel')
        if full:
            buttonSave = qt.QPushButton('Save')
        buttonLoad = qt.QPushButton('Load')

        # OptionLists and ButtonGroups
        # GroupBox can be generated from self.getGroupBox
        normOpts = ['No &normalization',
                    'Normalize &after average',
                    'Normalize &before average']
        xrangeOpts = ['&First curve in sequence',
                      'Active &curve',
                      '&Use equidistant x-range']
        # ButtonGroups
        normBG   = qt.QButtonGroup(self)
        xrangeBG = qt.QButtonGroup(self)
        # ComboBoxes
        normMeth = self.getComboBox(['(y-min(y))/trapz(max(y)-min(y),x)',
                                     'y/max(y)',
                                     '(y-min(y))/(max(y)-min(y))',
                                     '(y-min(y))/sum(max(y)-min(y))'])
        normMeth.setEnabled(False)
        motor0 = self.getComboBox(mList)
        motor1 = self.getComboBox(mList)
        self.optsDict = {
            'normalization' : normBG,
            'normalizationMethod' : normMeth,
            'xrange' : xrangeBG,
            'motor0': motor0,
            'motor1': motor1
        }
        # Subdivide into GroupBoxes
        normGroupBox = self.getGroupBox('Normalization', 
                                        normOpts, 
                                        normBG)
        xrangeGroupBox = self.getGroupBox('Interpolation x-range',
                                          xrangeOpts,
                                          xrangeBG)
        motorGroupBox = qt.QGroupBox('Motors')
        
        # Layouts
        mainLayout = qt.QVBoxLayout()
        buttonLayout = qt.QHBoxLayout()
        normLayout = qt.QHBoxLayout()
        motorLayout = qt.QGridLayout()
        if full: buttonLayout.addWidget(buttonSave)
        buttonLayout.addWidget(buttonLoad)
        buttonLayout.addWidget(qt.HorizontalSpacer())
        buttonLayout.addWidget(buttonOK)
        buttonLayout.addWidget(buttonCancel)
        normLayout.addWidget(qt.QLabel('Method:'))
        normLayout.addWidget(normMeth)
        motorLayout.addWidget(qt.QLabel('Motor 1:'),0,0)
        motorLayout.addWidget(motor0,0,1)
        motorLayout.addWidget(qt.QLabel('Motor 2:'),1,0)
        motorLayout.addWidget(motor1,1,1)
        motorGroupBox.setLayout(motorLayout)
        normGroupBox.layout().addLayout(normLayout)
        mainLayout.addWidget(normGroupBox)
        mainLayout.addWidget(xrangeGroupBox)
        mainLayout.addWidget(motorGroupBox)
        mainLayout.addLayout(buttonLayout)
        self.setLayout(mainLayout)
        
        # Connects
        if full:
            buttonOK.clicked.connect(self.accept)
        else:
            buttonOK.clicked[()].connect(self.saveOptionsAndClose)
        buttonCancel.clicked.connect(self.close)
        if full:
            buttonSave.clicked[()].connect(self.saveOptions)
        buttonLoad.clicked[()].connect(self.loadOptions)
        # Keep normalization method selector disabled 
        # when 'no normalization' selected
        normBG.button(0).toggled.connect(normMeth.setDisabled)

    def showEvent(self, event):
        self.saved = False
        qt.QDialog.showEvent(self, event)

    def updateMotorList(self, mList):
        for (key, obj) in self.optsDict.items():
            if key.startswith('motor') and isinstance(obj, qt.QComboBox):
                curr = obj.currentText()
                obj.clear()
                obj.addItems(mList)
                idx = obj.findText(curr)
                if idx < 0:
                    obj.setCurrentIndex(idx)
                else:
                    # Motor not found in Motorlist, set to default
                    obj.setCurrentIndex(idx)
    
    def getComboBox(self, itemList):
        tmp = qt.QComboBox()
        tmp.addItems(itemList)
        return tmp
        
    def getGroupBox(self, title, optionList, buttongroup=None):
        '''
        title : string
        optionList : List of strings
        buttongroup : qt.QButtonGroup
        
        Returns
        -------
        GroupBox of QRadioButtons build from a
        given optionList. If buttongroup is 
        specified, the buttons are organized in
        a QButtonGroup.
        '''
        first = True
        groupBox = qt.QGroupBox(title, None)
        gbLayout = qt.QVBoxLayout(None)
        gbLayout.addStretch(1)
        for (id, radioText) in enumerate(optionList):
            radio = qt.QRadioButton(radioText)
            gbLayout.addWidget(radio)
            if buttongroup:
                buttongroup.addButton(radio, id)
            if first:
                radio.setChecked(True)
                first = False
        groupBox.setLayout(gbLayout)
        return groupBox

    def normalizationMethod(self, ident):
        ret = None
        normDict = {
            'toMaximum'       : r'y/max(y)',
            'offsetAndMaximum': r'(y-min(y))/(max(y)-min(y))',
            'offsetAndCounts' : r'(y-min(y))/sum(max(y)-min(y))',
            'offsetAndArea'   : r'(y-min(y))/trapz(max(y)-min(y),x)'
        }
        for (name, eq) in normDict.iteritems():
            if ident == name:
                return eq
            if ident == eq:
                return name
        raise ValueError("'%s' not found.")

    def saveOptionsAndClose(self):
        if not self.saved:
            if not self.saveOptions():
                return
        self.accept()

    def saveOptions(self, filename=None):
        saveDir = PyMcaDirs.outputDir
        filter = ['PyMca (*.cfg)']
        if filename is None:
            try:
                filename = PyMcaFileDialogs.\
                            getFileList(parent=self,
                                filetypelist=filter,
                                message='Save XLD/XMCD Analysis Configuration',
                                mode='SAVE',
                                single=True)[0]
            except IndexError:
                # Returned list is empty
                return
            print 'saveOptions -- Filename: "%s"'%filename
        if len(filename) == 0:
            self.saved = False
            return False
        if not str(filename).endswith('.cfg'):
            filename += '.cfg'
        confDict = ConfigDict.ConfigDict()
        tmp = self.getOptions()
        for (key, value) in tmp.items():
            if key.startswith('Motor') and len(value) == 0:
                tmp[key] = 'None'
        confDict['XMCDOptions'] = tmp
        try:
            confDict.write(filename)
        except IOError:
            msg = qt.QMessageBox()
            msg.setWindowTitle('XLD/XMCD Options Error')
            msg.setText('Unable to write configuration to \'%s\''%filename)
            msg.exec_()
        self.saved = True
        return True

    def loadOptions(self):
        openDir = PyMcaDirs.outputDir
        filter = 'PyMca (*.cfg)'
        filename = qt.QFileDialog.\
                    getOpenFileName(self,
                                    'Load XLD/XMCD Analysis Configuration',
                                    openDir,
                                    filter)
        confDict = ConfigDict.ConfigDict()
        try:
            confDict.read(filename)
        except IOError:
            msg = qt.QMessageBox()
            msg.setTitle('XMCD Options Error')
            msg.setText('Unable to read configuration file \'%s\''%filename)
            return
        if 'XMCDOptions'not in confDict:
            return
        try:
            self.setOptions(confDict['XMCDOptions'])
        except ValueError as e:
            if DEBUG:
                print 'loadOptions -- int conversion failed:',
                print 'Invalid value for option \'%s\''%e
            else:
                msg = qt.QMessageBox()
                msg.setWindowTitle('XMCD Options Error')
                msg.setText('Configuration file \'%s\' corruted'%filename)
                msg.exec_()
                return
        except KeyError as e:
            if DEBUG:
                print 'loadOptions -- invalid identifier:',
                print 'option \'%s\' not found'%e
            else:
                msg = qt.QMessageBox()
                msg.setWindowTitle('XMCD Options Error')
                msg.setText('Configuration file \'%s\' corruted'%filename)
                msg.exec_()
                return
        self.saved = True

    def getOptions(self):
        ddict = {}
        for (option, obj) in self.optsDict.items():
            if isinstance(obj, qt.QButtonGroup):
                ddict[option] = obj.checkedId()
            elif isinstance(obj, qt.QComboBox):
                tmp = str(obj.currentText())
                if option == 'normalizationMethod':
                    tmp = self.normalizationMethod(tmp)
                if option.startswith('motor') and (not len(tmp)):
                    tmp = 'None'
                ddict[option] = tmp
            else:
                ddict[option] = 'None'
        return ddict

    def getMotors(self):
        motors = sorted([key for key in self.optsDict.keys() if key.startswith('motor')])
        return [str(self.optsDict[motor].currentText()) for motor in motors]

    def setOptions(self, ddict):
        for option in ddict.keys():
            obj = self.optsDict[option]
            if isinstance(obj, qt.QComboBox):
                name = ddict[option]
                if option == 'normalizationMethod':
                    name = self.normalizationMethod(name)
                if option.startswith('Motor') and name=='None':
                    name = ''
                id = obj.findText(qt.QString(name))
                obj.setCurrentIndex(id)
            elif isinstance(obj, qt.QButtonGroup):
                try:
                    id = int(ddict[option])
                except ValueError:
                    raise ValueError(option)
                button = self.optsDict[option].button(id)
                if type(button) == type(qt.QRadioButton()):
                        button.setChecked(True)

class XMCDScanWindow(sw.ScanWindow):

    plotModifiedSignal = qt.pyqtSignal()
    saveOptionsSignal = qt.pyqtSignal('QString')

    def __init__(self,
                 origin,
                 parent=None):
        sw.ScanWindow.__init__(self, 
                               parent, 
                               name='XLD/XMCD Analysis', 
                               specfit=None)
        self.plotWindow = origin
        self.scanWindowInfoWidget.hide()

        # Buttons to push spectra to main Window
        buttonWidget = qt.QWidget()
        buttonAdd = qt.QPushButton('Add', self)
        buttonAdd.setToolTip('Add active curve to scan window') 
        buttonReplace = qt.QPushButton('Replace', self)
        buttonReplace.setToolTip('Replace all curves in scan window with active curve') 
        buttonAddAll = qt.QPushButton('Add all', self)
        buttonAddAll.setToolTip('Add all curves to scan window') 
        buttonReplaceAll = qt.QPushButton('Replace all', self)
        buttonReplaceAll.setToolTip('Replace all curves in scan window with all curves') 
        buttonLayout = qt.QHBoxLayout(None)
        buttonLayout.setContentsMargins(0, 0, 0, 0)
        buttonLayout.setSpacing(5)
        # Show XAS & XMCD Buttons
        buttonLayout.addWidget(qt.HorizontalSpacer(self))
        buttonLayout.addWidget(buttonAdd)
        buttonLayout.addWidget(buttonAddAll)
        buttonLayout.addWidget(buttonReplace)
        buttonLayout.addWidget(buttonReplaceAll)
        buttonWidget.setLayout(buttonLayout)
#        self.mainLayout.addLayout(buttonLayout)
        self.mainLayout.addWidget(buttonWidget)
        
        buttonAdd.clicked.connect(self.add)
        buttonReplace.clicked.connect(self.replace)
        buttonAddAll.clicked.connect(self.addAll)
        buttonReplaceAll.clicked.connect(self.replaceAll)

        # Copy spectra from origin
        self.selectionDict = {'m':[], 'p':[]}
        self.curvesDict = {}
        self.optsDict = {
            'normAfterAvg'  : False,
            'normBeforeAvg' : False,
            'useActive'     : False,
            'equidistant'   : False,
            'normalizationMethod' : self.NormOffsetAndArea
        }
        self.xRange = None
        # Keep track of Averages, XMCD and XAS curves
        self.avgM = None
        self.avgP = None
        self.xmcd = None
        self.xas  = None

    def processOptions(self, options):
        tmp = { 'equidistant': False,
                'useActive': False,
                'normAfterAvg': False,
                'normBeforeAvg': False,
                'normalizationMethod': None
        }
        xRange = options['xrange']
        normalization = options['normalization']
        normMethod = options['normalizationMethod']
        # xRange Options. Default: Use first scan
        if xRange == 1:
            tmp['useActive']   = True
        elif xRange == 2:
            tmp['equidistant'] = True
        # Normalization Options. Default: No Normalization
        if normalization == 1:
            tmp['normAfterAvg']  = True
        elif normalization == 2:
            tmp['normBeforeAvg'] = True
        # Normalization Method. Default: offsetAndArea
        tmp['normalizationMethod'] = self.setNormalizationMethod(normMethod)
        # Trigger reclaculation
        self.optsDict = tmp
        msel = self.selectionDict['m']
        psel = self.selectionDict['p']
        self.processSelection(msel, psel)

    def setNormalizationMethod(self, fname):
        if fname == 'toMaximum':
            func = self.NormToMaximum
        elif fname == 'offsetAndMaximum':
            func = self.NormToOffsetAndMaximum
        elif fname == 'offsetAndCounts':
            func = self.NormOffsetAndCounts
        else:
            func = self.NormOffsetAndArea
        return func

    def NormToMaximum(self,x,y):
        ymax  = numpy.max(y)
        ynorm = y/ymax
        return ynorm

    def NormToOffsetAndMaximum(self,x,y):
        ynorm = y - numpy.min(y)
        ymax  = numpy.max(ynorm)
        ynorm /= ymax
        return ynorm

    def NormOffsetAndCounts(self, x, y):
        # Check for non-zero values?
        ynorm = y - numpy.min(y)
        ymax  = numpy.sum(ynorm)
        ynorm /= ymax
        return ynorm

    def NormOffsetAndArea(self,  x, y):
        # Check for non-zero values?
        ynorm = y - numpy.min(y)
        ymax  = numpy.trapz(ynorm,  x)
        ynorm /= ymax
        return ynorm

    def interpXRange(self,
                     xRange=None,
                     equidistant=False,
                     xRangeList=None):
        '''
        Input
        -----
        fromCurves : Bool
            Uses curves present in self.curvesDict
            if set to true. If set to false, an
            ndarray with equistant interpolation
            points is returned
        xRangeList : List
            List of ndarray from whose the overlap
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
            # Default xRangeList: curvesDict sorted for legends
            keys = sorted(self.curvesDict.keys())
            xRangeList = [self.curvesDict[k].x[0] for k in keys]
        if not len(xRangeList):
            if DEBUG:
                print 'interpXRange -- Nothing to do'
            return None
            
        num = 0
        xmin, xmax = self.plotWindow.getGraphXLimits()
        for x in xRangeList:
            if x.min() > xmin:
                xmin = x.min()
            if x.max() < xmax:
                xmax = x.max()
        if xmin >= xmax:
            raise ValueError('No overlap between curves')
            pass
        
        if equidistant:
            for x in xRangeList:
                curr = numpy.nonzero((x >= xmin) & 
                                     (x <= xmax))[0].size
                num = curr if curr>num else num
            # Exclude first and last point
            out = numpy.linspace(xmin, xmax, num, endpoint=False)[1:]
        else:
            if xRange is not None:
                x = xRange
            else:
                x = xRangeList[0]
            mask = numpy.nonzero((x > xmin) &
                                 (x < xmax))[0]
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
        self.curvesDict = self.copyCurves(msel + psel)
        
        if (len(self.curvesDict) == 0) or\
           ((len(self.selectionDict['m']) == 0) and\
           (len(self.selectionDict['p']) == 0)):
            # Nothing to do
            return

        # Make sure to use active curve when specified
        if self.optsDict['useActive']:
            # Get active curve
            active = self.plotWindow.getActiveCurve()
            if active:
                if DEBUG:
                    print 'processSelection -- xrange: use active'
                x, y, leg, info = active[0:4]
                xRange = self.interpXRange(xRange=x)
            else:
                return
        elif self.optsDict['equidistant']:
            if DEBUG:
                print 'processSelection -- xrange: use equidistant'
            xRange = self.interpXRange(equidistant=True)
        else:
            if DEBUG:
                print 'processSelection -- xrange: use first'
            xRange = self.interpXRange()
        activeLegend = self.plotWindow.graph.getActiveCurve(justlegend=True)
        if (not activeLegend) or (activeLegend not in self.curvesDict.keys()):
            # Use first curve in the series as xrange
            activeLegend = sorted(self.curvesDict.keys())[0]
        active = self.curvesDict[activeLegend]
        xlabel, ylabel = self.extractLabels(active.info)
        
        # Calculate averages and add them to the plot
        normalization = self.optsDict['normalizationMethod']
        normBefore = self.optsDict['normBeforeAvg']
        normAfter  = self.optsDict['normAfterAvg']
        for id in ['p','m']:
            sel = self.selectionDict[id]
            if not len(sel):
                continue
            xvalList = []
            yvalList = []
            for legend in sel:
                tmp = self.curvesDict[legend]
                if normBefore:
                    xVals = tmp.x[0]
                    yVals = normalization(xVals, tmp.y[0])
                else:
                    xVals = tmp.x[0]
                    yVals = tmp.y[0]
                xvalList.append(xVals)
                yvalList.append(yVals)
            avg_x, avg_y = self.specAverage(xvalList,
                                            yvalList,
                                            xRange)
            if normAfter:
                avg_y = normalization(avg_x, avg_y)
            avgName = 'avg_' + id
            info = {'xlabel': xlabel, 'ylabel': ylabel}
            self.addCurve(avg_x, avg_y, avgName, info)
            if id == 'p':
                self.avgP = self.dataObjectsList[-1]
            if id == 'm':
                self.avgM = self.dataObjectsList[-1]
            
        if (self.avgM and self.avgP):
            self.performXMCD()
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
        return out

    def specAverage(self, xarr, yarr, xRange=None):
        '''
        xarr : list
            List containing x-Values in 1-D numpy arrays
        yarr : list
            List containing y-Values in 1-D numpy arrays
        xRange : Numpy array
            x-Values used for interpolation. Must overlap
            with all arrays in xarr

        From the spectra given in xarr & yarr, the method
        determines the overlap in the x-range. For spectra
        with unequal x-ranges, the method interpolates all
        spectra on the values given in xRange and averages
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
                print 'specAverage -- invalid input!',
                print 'Array lengths do not match or are 0'
            return None, None 

        same = True
        if xRange == None:
            x0 = xarr[0]
        else:
            x0 = xRange
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

        if xRange != None:
            xmin0 = xRange.min()
            xmax0 = xRange.max()
        else:
            xmin0 = xsort[0][0]
            xmax0 = xsort[0][-1]
        if (not same) or (xRange == None):
            # Determine global xmin0 & xmax0
            for x in xsort:
                xmin = x.min()
                xmax = x.max()
                if xmin > xmin0:
                    xmin0 = xmin
                if xmax < xmax0:
                    xmax0 = xmax
            if xmax <= xmin:
                if DEBUG:
                    print 'specAverage -- ',
                    print 'No overlap between spectra!'
                return numpy.array([]), numpy.array([])

        # Clip xRange to maximal overlap in spectra
        if xRange is None:
            xRange = xsort[0]
        mask = numpy.nonzero((xRange>=xmin0) & 
                             (xRange<=xmax0))[0]
        xnew = numpy.take(xRange, mask)
        ynew = numpy.zeros(len(xnew))

        # Perform average
        for (x, y) in zip(xsort, ysort):
            if same:
                ynew += y  
            else:
                yinter = numpy.interp(xnew, x, y)
                ynew   += numpy.asarray(yinter)
        num = len(yarr)
        ynew /= num
        return xnew, ynew

    def extractLabels(self, info):
        xlabel = 'X'
        ylabel = 'Y'
        sel = info.get('selection', None)
        labelNames = info.get('LabelNames',[])
        if not len(labelNames):
            pass
        elif len(labelNames) == 2:
                [xlabel, ylabel] = labelNames
                print 'extractLabels -- using labelNames'
                print '\t',labelNames
        elif sel:
            xsel = sel.get('x',[])
            ysel = sel.get('y',[])
            if len(xsel) > 0:
                x = xsel[0]
                xlabel = labelNames[x]
            if len(ysel) > 0:
                y = ysel[0]
                ylabel = labelNames[y]
            print 'extractLabels -- xsel, ysel:'
            print '\txsel: ', xsel
            print '\tysel: ', ysel
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
        info = {'xlabel': xlabel, 'ylabel': ylabel}        
        self.addCurve(m.x[0], avg, xmcdLegend, info)
        self.xas = self.dataObjectsList[-1]

    def performXMCD(self):
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
        info = {'xlabel': xlabel, 'ylabel': ylabel}
        self.addCurve(p.x[0], diff, xmcdLegend, info)
        self.graph.mapToY2(' '.join([xmcdLegend, ylabel]))
        self._zoomReset()
        self.xmcd = self.dataObjectsList[-1]

    def selectionInfo(self, id, key):
        '''
        Convenience function to retrieve values
        from the info dictionaries of the curves
        stored selectionDict.
        '''
        sel = self.selectionDict[id]
        ret = '%s: '%id
        for legend in sel:
            curr = self.curvesDict[legend]
            value = curr.info.get(key, None)
            if value:
                ret = ' '.join([ret, value])
        return ret

    def _saveIconSignal(self):
        saveDir = PyMcaDirs.outputDir
        filter = ['spec File (*.spec)','Any File (*.*)']
        try:
            filename = PyMcaFileDialogs.\
                            getFileList(parent=self,
                                filetypelist=filter,
                                message='Save XMCD Analysis',
                                mode='SAVE',
                                single=True)[0]
        except IndexError:
            # Returned list is empty
            return

        ext = splitext(filename)[1]
        if not len(ext):
            ext = '.spec'
            filename += ext
        try:
            filehandle = open(filename, 'w')
        except IOError:
            msg = qt.QMessageBox(text="Unable to write to '%s'"%filename)
            msg.exec_()
            return
        # Keep plots in the order they were added!
        legends = self.dataObjectsList
        curves = [self.dataObjectsDict[leg] for leg in legends]
        yVals = [curve.y[0] for curve in curves]
        # xrange is the same for every curve
        xVals = [curves[0].x[0]]
        outArray = numpy.vstack([xVals, yVals]).T
        if len(outArray.shape) > 1:
            ncols = outArray.shape[1]
        else:
            ncols = 1   
        # Add Title to SaveFile
        # <selectionlegend>  <ScanNumber>
        title = ''
        tmpLegs = sorted(self.curvesDict.keys())
        if len(tmpLegs) > 0:
            title += self.curvesDict[tmpLegs[0]].info.get('selectionlegend','')
        for leg in tmpLegs[1:]:
            curr = self.curvesDict[leg]
            title += (' ' + curr.info.get('Key',''))
        title = 'XMCD Analysis ' + title
        delim = ' '
        header  = '#S 1 %s'%title + NEWLINE
        header += ('#U00 ' + self.selectionInfo('p', 'Key') + NEWLINE)
        header += ('#U01 ' + self.selectionInfo('m', 'Key') + NEWLINE)
        header += '#N %d'%ncols + NEWLINE
        if ext == '.spec':
            header += ('#L ' + self.getGraphXTitle() + '  ' + '  '.join(legends) + NEWLINE)
        else:
            header += ('#L ' + self.getGraphXTitle() + '  ' + delim.join(legends) + NEWLINE)

        filehandle.write(NEWLINE)
        filehandle.write(header)
        for line in outArray:
            tmp = delim.join(['%f'%num for num in line])
            filehandle.write(tmp + NEWLINE)
        filehandle.write(NEWLINE)
        filehandle.close()
        
        # Open filehandler for config file
        self.saveOptionsSignal.emit(splitext(filename)[0])
    
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
        

    
class XMCDMenu(qt.QMenu):
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

class XMCDTreeWidget(qt.QTreeWidget):

    selectionModifiedSignal = qt.pyqtSignal()

    def __init__(self,  parent, identifiers = ['p','m','d'], color=True):
        qt.QTreeWidget.__init__(self,  parent)
        self.identifiers = identifiers
        self.actionList  = []
        self.contextMenu = qt.QMenu('Perform',  self)
        self.color = color
        self.colorDict = {
            identifiers[0] : qt.QBrush(qt.QColor(220, 220, 255)),
            identifiers[1] : qt.QBrush(qt.QColor(255, 210, 210)),
            '': qt.QBrush(qt.QColor(255, 255, 255))
        }

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

    def getColumn(self, ncol, selectedOnly=False, convertType=str):
        '''
        Returns items in tree column ncol and converts them
        to convertType. If the conversion fails, the default
        type is a python string.
        
        If selectedOnly is set to True, only the selected
        the items of selected rows are returned.
        '''
        out = []
        convert = (convertType != str)
        if ncol > (self.columnCount()-1):
            if DEBUG:
                print 'getColum -- Selected column out of bounds'
            raise IndexError("Selected column '%d' out of bounds" % ncol)
            return out
        if selectedOnly:
            sel = self.selectedItems()
        else:
            root = self.invisibleRootItem()
            sel = [root.child(i) for i in range(root.childCount())]
        for item in sel:
            tmp = str(item.text(ncol))
            if convert:
                try:
                    tmp = convertType(tmp)
                except (TypeError, ValueError):
                    if convertType == float:
                        tmp = float('NaN')
                    else:
                        if DEBUG:
                            print 'getColum -- Conversion failed!'
                        raise TypeError
            out += [tmp]
        return out

    def build(self,  items,  headerLabels):
        '''
        (Re-) Builds the tree display

        headerLabels must be of type QStringList
        items must be of type [QStringList] (List of Lists)
        '''
        # Remember selection, then clear list
        sel = self.getColumn(1, True)
        self.clear()
        self.setHeaderLabels(headerLabels)
        for item in items:
#            treeItem = qt.QTreeWidgetItem(self,  item)
            treeItem = TreeWidgetItem(self,  item)
            if self.color:
                id = str(treeItem.text(0))
                for i in range(self.columnCount()):
                    treeItem.setBackground(i, self.colorDict[id])
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
            raise ValueError('XMCDTreeWidget: invalid identifer \'%s\''%id)
        sel = self.selectedItems()
        if id == self.identifiers[-1]:
            id = ''
        for item in sel:
            item.setText(0,id)
            if self.color:
                for i in range(self.columnCount()):
                    item.setBackground(i, self.colorDict[id])
        self.selectionModifiedSignal.emit()

    def setSelectionToSequence(self, seq=None, selectedOnly=False):
        '''
        Sets the id column (col 0) to seq. If
        sequence is None, a dialog window is 
        shown.
        '''
        chk = True
        if selectedOnly:
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
            if self.color:
                for i in range(self.columnCount()):
                    item.setBackground(i, self.colorDict[id])
        self.selectionModifiedSignal.emit()

    def clearSelection(self, selectedOnly=True):
        '''
        Empties the id column for the selected rows.
        '''
        if selectedOnly:
            sel = self.selectedItems()
        else:
            root = self.invisibleRootItem()
            sel = [root.child(i) for i in range(root.childCount())]
        for item in sel:
            item.setText(0,'')
            if self.color:
                for i in range(self.columnCount()):
                    item.setBackground(i, self.colorDict[''])
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


class XMCDWidget(qt.QWidget):

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
            XLD/XMCD Analysis
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
        qt.QWidget.__init__(self, parent)
        self.plotWindow = plotWindow
        self.legendList = []
        self.motorsList = []
        self.infoList   = []
        # Set self.plotWindow before calling self._setLists!
        self._setLists()
        self.motorNamesList = [''] + self._getAllMotorNames()
        self.motorNamesList.sort()
        self.numCurves = len(self.legendList)
        #self.cBoxList = []
        self.ScanWindow = XMCDScanWindow(origin=plotWindow, 
                                              parent=None)
        self.optsWindow = XMCDOptions(self, self.motorNamesList)
                                              
        self.selectionDict = {'d': [],
                              'p': [],
                              'm': []}
        self.beamline = beamline
        self.setSizePolicy(qt.QSizePolicy.MinimumExpanding, 
                           qt.QSizePolicy.Expanding)

        self.setWindowTitle("XLD/XMCD Analysis")
        updatePixmap  = qt.QPixmap(IconDict["reload"])
        buttonUpdate  = qt.QPushButton(qt.QIcon(updatePixmap), '', self)
        buttonOptions = qt.QPushButton('Options', self)
        #for i in range(nSelectors):
        #    cBox = qt.QComboBox(self)
        #    cBox.addItems(self.motorNamesList)
        #    cBox.currentIndexChanged['QString'].connect(self.updateTree)
        #    self.cBoxList += [cBox]

        self.ident = 'Key'
        for ddict in self.infoList:
            if self.ident not in ddict.keys():
                self.ident = 'selectionlegend'
                break
            elif not len(ddict[self.ident]):
                self.ident = 'selectionlegend'
                break

        self.list = XMCDTreeWidget(self)
        labels = ['Legend'] + nSelectors*['']
        ncols  = len(labels)
        self.list.setColumnCount(ncols)
        self.list.setHeaderLabels(labels)
        self.list.setSortingEnabled(True)
        self.list.setSelectionMode(
            qt.QAbstractItemView.ExtendedSelection)
        listContextMenu = XMCDMenu(None)
        listContextMenu.setActionList(
              [('Perform analysis', self.triggerXMCD),
               ('$SEPERATOR', None),
               ('Set as p', self.setAsP),
               ('Set as m', self.setAsM),
               ('Enter sequence', self.list.setSelectionToSequence),
               ('Remove selection', self.list.clearSelection),
               ('$SEPERATOR', None),
               ('Invert selection', self.list.invertSelection), 
               ('Remove curve(s)', self.removeCurve_)])
        self.list.setContextMenu(listContextMenu)
        self.expCBox = qt.QComboBox(self)
        self.expCBox.addItems(
                        ['Select Experiment',
                         'ID08: Linear Dichorism',
                         'ID08: Linear Dichorism (old)',
                         'ID08: XMCD',
                         'ID08: XMCD (old)',
                         'Add new experiment'])
        self.expCBox.insertSeparator(5)
        
        self.experimentsDict = {
            'Select Experiment': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': '',
                  'motor1': ''
            },
            'ID08: XMCD': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': 'phaseD',
                  'motor1': 'magnet'
            },
            'ID08: XMCD (old)': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': 'PhaseD',
                  'motor1': 'oxPS'
            },
            'ID08: Linear Dichorism (old)': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': 'PhaseD',
                  'motor1': ''                     
            },
            'ID08: Linear Dichorism': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': 'phaseD',
                  'motor1': ''                     
            }
        }
        
        #cBoxLayout = qt.QHBoxLayout(None)
        #self.cBoxWidget = qt.QWidget()
        #cBoxLayout.addWidget(qt.HorizontalSpacer(self))
        #cBoxLayout.addWidget(
        #        qt.QLabel('Selected motor(s):',  self))
        #for cBox in self.cBoxList:
        #    cBoxLayout.addWidget(cBox) 
        #self.cBoxWidget.setLayout(cBoxLayout)
        #cBoxLayout.setContentsMargins(0,0,0,0)
        
        topLayout  = qt.QHBoxLayout()
        topLayout.addWidget(buttonUpdate)
        topLayout.addWidget(buttonOptions)
        topLayout.addWidget(qt.HorizontalSpacer(self))
        topLayout.addWidget(self.expCBox)

        leftLayout = qt.QGridLayout()
        leftLayout.setContentsMargins(1, 1, 1, 1)
        leftLayout.setSpacing(2)
        #leftLayout.addWidget(self.cBoxWidget, 2, 0)
        leftLayout.addWidget(self.list, 1, 0)
        leftLayout.addLayout(topLayout, 0, 0)
        leftWidget = qt.QWidget(self)
        leftWidget.setLayout(leftLayout)
        
        self.splitter = qt.QSplitter(qt.Qt.Horizontal, self)
        self.splitter.addWidget(leftWidget)
        self.splitter.addWidget(self.ScanWindow)
        
        mainLayout = qt.QVBoxLayout()
        mainLayout.setContentsMargins(0,0,0,0)
        mainLayout.addWidget(self.splitter)
        self.setLayout(mainLayout)

        # Connects
        self.expCBox.currentIndexChanged['QString'].connect(self.updateTree)
        self.expCBox.currentIndexChanged['QString'].connect(self.selectExperiment)
        self.list.selectionModifiedSignal.connect(self.updateSelectionDict)
        self.setSelectionSignal.connect(self.ScanWindow.processSelection)
        self.ScanWindow.saveOptionsSignal.connect(self.optsWindow.saveOptions)
        self.optsWindow.accepted[()].connect(self.updateTree)
        buttonUpdate.clicked.connect(self.updatePlots)
        buttonOptions.clicked.connect(self.showOptionsWindow)

        self.updateTree()
        self.list.sortByColumn(1, qt.Qt.AscendingOrder)
        self._setBeamlineSpecific(self.beamline)

    def addExperiment(self):
        exp, chk = qt.QInputDialog.\
                        getText(self,
                                'Configure new experiment',
                                'Enter experiment title',
                                qt.QLineEdit.Normal, 
                                'ID00: <Title>')
        if chk and (not exp.isEmpty()):
            print 'Here'
            exp = str(exp)
            opts = XMCDOptions(self, self.motorNamesList, False)
            if opts.exec_():
                self.experimentsDict[exp] = opts.getOptions()
                cBox = self.expCBox
                new = [cBox.itemText(i) for i in range(cBox.count())][0:-3]
                new += [exp]
                new.append('Add new experiment')
                cBox.clear()
                cBox.addItems(new)
                cBox.insertSeparator(len(new)-1)
                idx = cBox.findText([exp][0])
                if idx < 0:
                    cBox.setCurrentIndex(0)
                else:
                    cBox.setCurrentIndex(idx)
            # What if save is canceled?
            opts.destroy()
            idx = self.expCBox.findText(exp)
            if idx < 0:
                idx = 0
            self.expCBox.setCurrentIndex(idx)
        
    def showOptionsWindow(self):
        if self.optsWindow.exec_():
            options = self.optsWindow.getOptions()
            self.ScanWindow.processOptions(options)

# Implement new assignment routines here BEGIN
    def selectExperiment(self, exp):
        exp = str(exp)
        print 'selectExperiment -- "%s"'%exp
        try:
            print '\t',self.experimentsDict[exp]
        except:
            print '\t','no experimentsDict'
        if exp == 'Add new experiment':
            self.addExperiment()
            self.updateTree()
            print '\tHere!'
        elif exp in self.experimentsDict:
            try:
                self.optsWindow.setOptions(self.experimentsDict[exp])
            except ValueError:
                self.optsWindow.setOptions(
                        self.experimentsDict['Select Experiment'])
                return
            self.updateTree()
            # Get motor values from tree
            values0 = numpy.array(
                        self.list.getColumn(2, convertType=float))
            values1 = numpy.array(
                        self.list.getColumn(3, convertType=float))
            # Calculate p/m selection
            if exp.startswith('ID08: Linear Dichorism'):
                values = values0
                mask = numpy.where(numpy.isfinite(values))[0]
                minmax = values.take(mask)
                if len(minmax):
                    vmin = minmax.min()
                    vmax = minmax.max()
                    vpivot = .5 * (vmax + vmin)
                else:
                    values = numpy.array(
                                [float('NaN')]*len(self.legendList))
            elif exp.startswith('ID08: XMCD'):
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
# Implement new assignment routines here END

    def triggerXMCD(self):
        msel = self.selectionDict['m']
        psel = self.selectionDict['p']
        self.ScanWindow.processSelection(msel, psel)

    def removeCurve_(self):
        sel = self.list.getColumn(1, 
                                  selectedOnly=True,
                                  convertType=str)
        # Convert from scan number to legend if needed
        if self.ident == 'Key':
            legends = []
            for item in sel:
                for (idx, info) in enumerate(self.infoList):
                    if item == info['Key']:
                        legends += [self.legendList[idx]]
        else:
            legends = sel
            
        if DEBUG:
            print 'removeCurve_ -- sel(ection):'
            print '\t', sel
            print 'removeCurve_ -- legends:'
            print '\t', legends
        for legend in legends:
            self.plotWindow.removeCurve(legend)
            for selection in self.selectionDict.values():
                if legend in selection:
                    selection.remove(legend)
            # Remove from XMCDScanWindow.curvesDict
            if legend in self.ScanWindow.curvesDict.keys():
                del(self.ScanWindow.curvesDict[legend])
            # Remove from XMCDScanWindow.selectionDict
            for selection in self.ScanWindow.selectionDict.values():
                if legend in selection:
                    selection.remove(legend)
        self.updatePlots()

    def updateSelectionDict(self):
        selDict = self.list.getSelection()
        # self.selectionDict -> Uses ScanNumbers instead of legends...
        newDict = {}
        if self.ident == 'Key':
            scanNumberList = [info['Key'] for info in self.infoList]
        for (id, selList) in selDict.items():
            if id not in newDict.keys():
                newDict[id] = []
            if self.ident == 'selectionlegend':
                for legend in selList:
                    newDict[id] += [legend]
            else:
                for scanNumber in selList:
                    idx = scanNumberList.index(scanNumber)
                    legend = self.legendList[idx]
                    newDict[id] += [legend]
        self.selectionDict = newDict
        self.setSelectionSignal.emit(self.selectionDict['m'],
                                     self.selectionDict['p'])

    def updatePlots(self,
                    newLegends = None,  
                    newMotorValues = None):
        # Check if curves in plotWindow changed..
        curves = self.plotWindow.getAllCurves(just_legend=True)
        if curves == self.legendList:
            # ..if not, just replot to account for zoom
            self.triggerXMCD()
            return
        self._setLists()
        self.motorNamesList = [''] + self._getAllMotorNames()
        self.motorNamesList.sort()
        self.optsWindow.updateMotorList(self.motorNamesList)
        self.updateTree()
        experiment = str(self.expCBox.currentText())
        if experiment != 'Select Experiment':
            self.selectExperiment(experiment)
        return

    def updateTree(self):
        mList  = self.optsWindow.getMotors()
        if self.ident == 'Key':
            labels = ["Sel",'S#'] + mList
        else:
            labels = ["Sel",'Legend'] + mList
        items  = []
        for i in range(len(self.legendList)):
            legend = self.legendList[i]
            values = self.motorsList[i]
            info = self.infoList[i]
            selection = ''
            for (id,v) in self.selectionDict.items():
                if (legend in v) and (id != 'd'):
                    selection = id
                    break
            if self.ident == 'Key':
                tmp = qt.QStringList([selection, info['Key']])
            else:
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
        self.legendList = [leg for (xvals, yvals,  leg,  info) in curves] 
        self.infoList = [info for (xvals, yvals,  leg,  info) in curves] 
        self.motorsList = self._convertInfoDictionary(self.infoList)

    def _setBeamlineSpecific(self, beamline):
        '''
        beamline : python str
            Beamline identifier, all upper case
            (e.g. ID08, ID12, ..)
        '''
        options = [str(self.expCBox.itemText(i))\
                   for i in range(self.expCBox.count())]
        for (i, option) in enumerate(options):
            if option.startswith(beamline):
                self.expCBox.setCurrentIndex(i)
                self.expCBox.activated[qt.QString].emit(option)
                break

def main():
    import sys,  numpy
    app = qt.QApplication(sys.argv)
    swin = sw.ScanWindow()
    swin.setWindowTitle('TESTSCANWINDOW')
    swin.setGraphXTitle('Blubb')
    info0 = {'xlabel': 'foo', 'ylabel': 'arb', 'MotorNames': 'oxPS Motor11 Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3',  'MotorValues': '1 8.69271399699 21.9836418539 0.198068826612 0.484475455792 0.350252217264 0.663925270933 0.813033264421 0.221149410218 0.593188258866 0.678010392881 0.267389247833 0.677890617858'}
    info1 = {'MotorNames': 'PhaseD oxPS Motor16 Motor15 Motor14 Motor13 Motor12 Motor11 Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3', 'MotorValues': '0.470746882688 -0.695816070299 0.825780811755 0.25876374531 0.739264467436 0.090842892619 2 0.213445659833 0.823400550314 0.020278096857 0.568744021322 0.85378115537 0.696730386891 0.269196313956 0.793293334395 0.769216567757 0.959092709527 0.0109264683697 0.538264972553'}
    info2 = {'MotorNames': 'PhaseD oxPS Motor10 Motor8 Motor9 Motor4 Motor5 Motor6 Motor7 Motor0 Motor1 Motor2 Motor3',  'MotorValues': '2 0.44400576644 0.613870067852 0.901968648111 0.319768771085 0.571432278628 0.278675836163 0.154436774878 0.416231999332 0.294201017231 0.813913587748 0.577572903105 0.869045182568'}
    x = numpy.arange(100.,1100.)
    y0 =  10 * x + 10000. * numpy.exp(-0.5*(x-500)*(x-500)/400) + 1500 * numpy.random.random(1000.)
    y1 =  10 * x + 10000. * numpy.exp(-0.5*(x-600)*(x-600)/400) + 1500 * numpy.random.random(1000.)
    y2 =  10 * x + 10000. * numpy.exp(-0.5*(x-400)*(x-400)/400) + 1500 * numpy.random.random(1000.)
    y2[320:322] = 50000.
    
    swin.newCurve(x, y2, legend="Curve2", xlabel='ene_st2', ylabel='zratio2', info=info2, replot=False, replace=False)
    swin.newCurve(x, y0, legend="Curve0", xlabel='ene_st0', ylabel='zratio0', info=info0, replot=False, replace=False)
    swin.newCurve(x, y1, legend="Curve1", xlabel='ene_st1', ylabel='zratio1', info=info1, replot=False, replace=False)
    
    # info['Key'] is overwritten when using newCurve
    swin.dataObjectsDict['Curve2 zratio2'].info['Key'] = '1.1'
    swin.dataObjectsDict['Curve0 zratio0'].info['Key'] = '34.1'
    swin.dataObjectsDict['Curve1 zratio1'].info['Key'] = '123.1'

    w = XMCDWidget(None, swin, 'ID08', nSelectors = 2)
    w.show()
    app.exec_()

if __name__ == '__main__':
    main()
