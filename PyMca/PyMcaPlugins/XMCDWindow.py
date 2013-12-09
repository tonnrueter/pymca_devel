import numpy, copy
from os.path import splitext, basename, dirname, exists, join as pathjoin
from re import search as regexpSearch
from PyMca.PyMca_Icons import IconDict
from PyMca import PyMcaDirs, PyMcaFileDialogs
from PyMca import ConfigDict
from PyMca import PyMcaQt as qt
from PyMca import specfilewrapper as specfile

from PyMca import ScanWindow as sw

DEBUG = 0
if DEBUG:
    numpy.set_printoptions(threshold=50)

NEWLINE = '\n'
class TreeWidgetItem(qt.QTreeWidgetItem):
    
    __legendColumn = 1
    
    def __init__(self, parent, itemList):
        qt.QTreeWidgetItem.__init__(self, parent, itemList)
        
    def __lt__(self, other):
        col = self.treeWidget().sortColumn()
        val      = self.text(col)
        valOther = other.text(col)
        if val == '---':
                ret = True
        elif col > self.__legendColumn:
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
        buttonOK.setToolTip('Accept the configuration')
        buttonCancel = qt.QPushButton('Cancel')
        buttonCancel.setToolTip('Return to XMCD Analysis\nwithout changes')
        if full:
            buttonSave = qt.QPushButton('Save')
            buttonSave.setToolTip('Save configuration to *.cfg-File')
        buttonLoad = qt.QPushButton('Load')
        buttonLoad.setToolTip('Load existing configuration from *.cfg-File')

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
        normMeth = qt.QComboBox()
        normMeth.addItems(['(y-min(y))/trapz(max(y)-min(y),x)',
                           'y/max(y)',
                           '(y-min(y))/(max(y)-min(y))',
                           '(y-min(y))/sum(max(y)-min(y))'])
        normMeth.setEnabled(False)
        self.optsDict = {
            'normalization' : normBG,
            'normalizationMethod' : normMeth,
            'xrange' : xrangeBG
        }
        for idx in range(5):
            # key: motor0, motor1, ...
            key = 'motor%d'%idx
            tmp = qt.QComboBox()
            tmp.addItems(mList)
            self.optsDict[key] = tmp
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
        for idx in range(5):
            label = qt.QLabel('Motor %d:'%(idx+1))
            cbox  = self.optsDict['motor%d'%idx]
            motorLayout.addWidget(label,idx,0)
            motorLayout.addWidget(cbox,idx,1)
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
        # Plugin does not destroy Options Window when accepted
        # Reset self.saved manually
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
        for (idx, radioText) in enumerate(optionList):
            radio = qt.QRadioButton(radioText)
            gbLayout.addWidget(radio)
            if buttongroup:
                buttongroup.addButton(radio, idx)
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
            if DEBUG:
                print('saveOptions -- Filename: "%s"' % filename)
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
                print('loadOptions -- int conversion failed:',)
                print('Invalid value for option \'%s\'' % e)
            else:
                msg = qt.QMessageBox()
                msg.setWindowTitle('XMCD Options Error')
                msg.setText('Configuration file \'%s\' corruted' % filename)
                msg.exec_()
                return
        except KeyError as e:
            if DEBUG:
                print('loadOptions -- invalid identifier:',)
                print('option \'%s\' not found' % e)
            else:
                msg = qt.QMessageBox()
                msg.setWindowTitle('XMCD Options Error')
                msg.setText('Configuration file \'%s\' corruted' % filename)
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
        motors = sorted([key for key in self.optsDict.keys()\
                         if key.startswith('motor')])
        return [str(self.optsDict[motor].currentText()) \
                for motor in motors]

    def setOptions(self, ddict):
        for option in ddict.keys():
            obj = self.optsDict[option]
            if isinstance(obj, qt.QComboBox):
                name = ddict[option]
                if option == 'normalizationMethod':
                    name = self.normalizationMethod(name)
                if option.startswith('Motor') and name == 'None':
                    name = ''
                idx = obj.findText(qt.QString(name))
                obj.setCurrentIndex(idx)
            elif isinstance(obj, qt.QButtonGroup):
                try:
                    idx = int(ddict[option])
                except ValueError:
                    raise ValueError(option)
                button = self.optsDict[option].button(idx)
                if type(button) == type(qt.QRadioButton()):
                        button.setChecked(True)

class XMCDScanWindow(sw.ScanWindow):

    plotModifiedSignal = qt.pyqtSignal()
    saveOptionsSignal  = qt.pyqtSignal('QString')

    def __init__(self,
                 origin,
                 parent=None):
        '''
        Input
        -----
        
        origin : ScanWindow instance
            Plot window containing the data on which
            the analysis is performed
        parent : QWidgit instance
        '''
        sw.ScanWindow.__init__(self, 
                               parent, 
                               name='XLD/XMCD Analysis', 
                               specfit=None)
        self.plotWindow = origin
        self.scanWindowInfoWidget.hide()

        # Buttons to push spectra to main Window
        buttonWidget = qt.QWidget()
        buttonAdd = qt.QPushButton('Add', self)
        buttonAdd.setToolTip('Add active curve to main window') 
        buttonReplace = qt.QPushButton('Replace', self)
        buttonReplace.setToolTip(
            'Replace all curves in main window '
           +'with active curve in analysis window')
        buttonAddAll = qt.QPushButton('Add all', self)
        buttonAddAll.setToolTip(
            'Add all curves in analysis window '
           +'to main window')
        buttonReplaceAll = qt.QPushButton('Replace all', self)
        buttonReplaceAll.setToolTip(
            'Replace all curves in main window '
           +'with all curves from analysis window')
        self.graphBottomLayout.addWidget(qt.HorizontalSpacer())
        self.graphBottomLayout.addWidget(buttonAdd)
        self.graphBottomLayout.addWidget(buttonAddAll)
        self.graphBottomLayout.addWidget(buttonReplace)
        self.graphBottomLayout.addWidget(buttonReplaceAll)
        
        buttonAdd.clicked.connect(self.add)
        buttonReplace.clicked.connect(self.replace)
        buttonAddAll.clicked.connect(self.addAll)
        buttonReplaceAll.clicked.connect(self.replaceAll)

        # Copy spectra from origin
        self.selectionDict = {'A':[], 'B':[]}
        self.curvesDict = {}
        self.optsDict = {
            'normAfterAvg'  : False,
            'normBeforeAvg' : False,
            'useActive'     : False,
            'equidistant'   : False,
            'normalizationMethod' : self.NormOffsetAndArea
        }
        self.xRange = None
        
        # Keep track of Averages, XMCD and XAS curves by label
        self.avgA = None
        self.avgB = None
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
        groupA = self.selectionDict['A']
        groupB = self.selectionDict['B']
        self.processSelection(groupA, groupB)

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
        ynorm = y - numpy.min(y)
        ymax  = numpy.sum(ynorm)
        ynorm /= ymax
        return ynorm

    def NormOffsetAndArea(self,  x, y):
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
        xRange : ndarray
            x-range on which all curves are interpolated
            If set to none, the first curve in xRangeList
            is used
        equidistant : Bool
            Determines equidistant xrange on which
            all curves are interpolated
        xRangeList : List
            List of ndarray from which the overlap
            is determined. If set to none, self.curvesDict
            is used.
        
        Determines the interpolation x-range used for XMCD
        Analysis specified by the provided parameters. The
        interpolation x-range is limited by the maximal 
        overlap between the curves present in the plot window.
        
        Returns
        -------
        out : numpy array
            x-range between xmin and xmax containing n points
        '''
        if not xRangeList:
            # Default xRangeList: curvesDict sorted for legends
            keys = sorted(self.curvesDict.keys())
            xRangeList = [self.curvesDict[k].x[0] for k in keys]
        if not len(xRangeList):
            if DEBUG:
                print('interpXRange -- Nothing to do')
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
            # Ensure monotonically increasing x-range
            if not numpy.all(numpy.diff(x)>0.):
                mask = numpy.nonzero(numpy.diff(x)>0.)[0]
                x = numpy.take(x, mask)
            # Exclude the endpoints
            mask = numpy.nonzero((x > xmin) &
                                 (x < xmax))[0]
            out = numpy.sort(numpy.take(x, mask))
        if DEBUG:
            print('interpXRange -- Resulting xrange:')
            print('\tmin = %f' % out.min())
            print('\tmax = %f' % out.max())
            print('\tnum = %f' % len(out))
        return out

    def processSelection(self, groupA, groupB):
        '''
        Input
        -----
        groupA, groupB : Lists of strings
            Contain the legends of curves selected
            to Group A resp. B
        '''
        # Clear analysis window
        all = self.getAllCurves(just_legend=True)
        self.removeCurves(all)
        self.avgB, self.avgA = None, None
        self.xas, self.xmcd  = None, None
        
        self.selectionDict['A'] = groupA[:]
        self.selectionDict['B'] = groupB[:]
        self.curvesDict = self.copyCurves(groupA + groupB)
        
        if (len(self.curvesDict) == 0) or\
           ((len(self.selectionDict['A']) == 0) and\
           (len(self.selectionDict['B']) == 0)):
            # Nothing to do
            return

        # Make sure to use active curve when specified
        if self.optsDict['useActive']:
            # Get active curve
            active = self.plotWindow.getActiveCurve()
            if active:
                if DEBUG:
                    print('processSelection -- xrange: use active')
                x, y, leg, info = active[0:4]
                xRange = self.interpXRange(xRange=x)
            else:
                return
        elif self.optsDict['equidistant']:
            if DEBUG:
                print('processSelection -- xrange: use equidistant')
            xRange = self.interpXRange(equidistant=True)
        else:
            if DEBUG:
                print('processSelection -- xrange: use first')
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
        for idx in ['A','B']:
            sel = self.selectionDict[idx]
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
            avgName = 'avg_' + idx
            info = {'xlabel': xlabel, 'ylabel': ylabel}
            self.addCurve(avg_x, avg_y, avgName, info)
            if idx == 'A':
                self.avgA = self.dataObjectsList[-1]
            if idx == 'B':
                self.avgB = self.dataObjectsList[-1]
            
        if (self.avgA and self.avgB):
            self.performXMCD()
            self.performXAS()

    def copyCurves(self, selection):
        '''
        Input
        -----
        selection : List
            Contains legends of curves to be processed
        
        Creates a deep copy of the curves present in the
        plot window. In order to avoid interpolation
        errors later on, it is ensured that the xranges
        of the data is strictly monotonically increasing.
        
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
                tmp = copy.deepcopy(tmp)
                xarr, yarr = tmp.x, tmp.y
                #if len(tmp.x) == len(tmp.y):
                xprocArr, yprocArr = [], []
                for (x,y) in zip(xarr,yarr):
                    # Sort
                    idx = numpy.argsort(x, kind='mergesort')
                    xproc = numpy.take(x, idx)
                    yproc = numpy.take(y, idx)
                    # Ravel, Increase
                    xproc = xproc.ravel()
                    idx = numpy.nonzero((xproc[1:] > xproc[:-1]))[0]
                    xproc = numpy.take(xproc, idx)
                    yproc = numpy.take(yproc, idx)
                    xprocArr += [xproc]
                    yprocArr += [yproc]
                tmp.x = xprocArr
                tmp.y = yprocArr
                out[legend] = tmp
            else:
                # TODO: Errorhandling, curve not found
                if DEBUG:
                    print("copyCurves -- Retrieved none type curve")
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
                print('specAverage -- invalid input!')
                print('Array lengths do not match or are 0')
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
                    print('specAverage -- ')
                    print('No overlap between spectra!')
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
        elif sel:
            xsel = sel.get('x',[])
            ysel = sel.get('y',[])
            if len(xsel) > 0:
                x = xsel[0]
                xlabel = labelNames[x]
            if len(ysel) > 0:
                y = ysel[0]
                ylabel = labelNames[y]
        return xlabel, ylabel

    def performXAS(self):
        keys = self.dataObjectsDict.keys()
        if (self.avgA in keys) and (self.avgB in keys):
            a = self.dataObjectsDict[self.avgA]
            b = self.dataObjectsDict[self.avgB]
        else:
            if DEBUG:
                print('performXAS -- Data not found: ')
                print('\tavg_m = %f' % self.avgA)
                print('\tavg_p = %f' % self.avgB)
            return
        if numpy.all( a.x[0] == b.x[0] ):
            avg = .5*(b.y[0] + a.y[0])
        else:
            if DEBUG:
                print('performXAS -- x ranges are not the same! ')
                print('Force interpolation')
            avg = self.performAverage([a.x[0], b.x[0]],
                                      [a.y[0], b.y[0]],
                                       b.x[0])
        xmcdLegend = 'XAS'
        xlabel, ylabel = self.extractLabels(a.info)
        info = {'xlabel': xlabel, 'ylabel': ylabel}        
        self.addCurve(a.x[0], avg, xmcdLegend, info)
        self.xas = self.dataObjectsList[-1]

    def performXMCD(self):
        keys = self.dataObjectsDict.keys()
        if (self.avgA in keys) and (self.avgB in keys):
            a = self.dataObjectsDict[self.avgA]
            b = self.dataObjectsDict[self.avgB]
        else:
            if DEBUG:
                print('performXMCD -- Data not found:')
            return
        if numpy.all( a.x[0] == b.x[0] ):
            diff = b.y[0] - a.y[0]
        else:
            if DEBUG:
                print('performXMCD -- x ranges are not the same! ')
                print('Force interpolation using p Average xrange')
            # Use performAverage d = 2 * avg(y1, -y2)
            # and force interpolation on p-xrange
            diff = 2. * self.performAverage([a.x[0], b.x[0]],
                                            [-a.y[0], b.y[0]],
                                            b.x[0])
        xmcdLegend = 'XMCD'
        xlabel, ylabel = self.extractLabels(a.info)
        info = {'xlabel': xlabel, 'ylabel': ylabel}
        self.addCurve(b.x[0], diff, xmcdLegend, info)
        self.graph.mapToY2(' '.join([xmcdLegend, ylabel]))
        self._zoomReset()
        self.xmcd = self.dataObjectsList[-1]

    def selectionInfo(self, idx, key):
        '''
        Convenience function to retrieve values
        from the info dictionaries of the curves
        stored selectionDict.
        '''
        sel = self.selectionDict[idx]
        ret = '%s: '%idx
        for legend in sel:
            curr = self.curvesDict[legend]
            value = curr.info.get(key, None)
            if value:
                ret = ' '.join([ret, value])
        return ret

    def _saveIconSignal(self):
        saveDir = PyMcaDirs.outputDir
        filter = 'spec File (*.spec);;Any File (*.*)'
        try:
            (filelist, append, comment) = getSaveFileName(parent=self,
                                                          caption='Save XMCD Analysis',
                                                          filter=filter,
                                                          directory=saveDir)
            filename = filelist[0]
        except IndexError:
            # Returned list is empty
            return
        except ValueError:
            # Returned list is empty
            return

        if append:
            specf  = specfile.Specfile(filename)
            scanno = specf.scanno() + 1
        else:
            scanno = 1

        ext = splitext(filename)[1]
        if not len(ext):
            ext = '.spec'
            filename += ext
        try:
            if append:
                sepFile = splitext(basename(filename))
                sepFileName = sepFile[0] + '_%.2d'%scanno + sepFile[1]
                sepFileName = pathjoin(dirname(filename),sepFileName)
                if scanno == 2:
                    # Case: Scan appended to file containing
                    # a single scan. Make sure, that the first 
                    # scan is also written to seperate file and
                    # the corresponding cfg-file is copied
                    # 1. Create filename of first scan
                    sepFirstFileName = sepFile[0] + '_01' + sepFile[1]
                    sepFirstFileName = pathjoin(dirname(filename),sepFirstFileName)
                    # 2. Guess filename of first config
                    confname = sepFile[0] + '.cfg'
                    confname = pathjoin(dirname(filename),confname)
                    # 3. Create new filename of first config
                    sepFirstConfName = sepFile[0] + '_01' + '.cfg'
                    sepFirstConfName = pathjoin(dirname(filename),sepFirstConfName)
                    # Copy contents
                    firstSeperateFile = open(sepFirstFileName, 'wb')
                    firstSeperateConf = open(sepFirstConfName, 'wb')
                    filehandle = open(filename, 'rb')
                    confhandle = open(confname, 'rb')
                    firstFile = filehandle.read()
                    firstConf = confhandle.read()
                    firstSeperateFile.write(firstFile)
                    firstSeperateConf.write(firstConf)
                    firstSeperateFile.close()
                    firstSeperateConf.close()
                filehandle = open(filename, 'ab')
                seperateFile = open(sepFileName, 'wb')
            else:
                filehandle = open(filename, 'wb')
                seperateFile = None
        except IOError:
            msg = qt.QMessageBox(text="Unable to open '%s'"%filename)
            msg.exec_()
            return
            
        title = ''
        legends = self.dataObjectsList
        tmpLegs = sorted(self.curvesDict.keys())
        if len(tmpLegs) > 0:
            title += self.curvesDict[tmpLegs[0]].info.get('selectionlegend','')
            # Keep plots in the order they were added!
            curves = [self.dataObjectsDict[leg] for leg in legends]
            yVals = [curve.y[0] for curve in curves]
            # xrange is the same for every curve
            xVals = [curves[0].x[0]]
        else:
            yVals = []
            xVals = []
        outArray = numpy.vstack([xVals, yVals]).T
        if not len(outArray):
            ncols = 0
        elif len(outArray.shape) > 1:
            ncols = outArray.shape[1]
        else:
            ncols = 1
        delim = ' '
        title = 'XMCD Analysis ' + title
        header  = '#S %d %s'%(scanno,title) + NEWLINE
        header += ('#U00 Selected in Group ' +\
                   self.selectionInfo('A', 'Key') + NEWLINE)
        header += ('#U01 Selected in Group ' +\
                   self.selectionInfo('B', 'Key') + NEWLINE)
        # Write Comments
        if len(comment) > 0:
            header += ('#U02 User commentary:' + NEWLINE)
            lines = comment.splitlines()[:97]
            for (idx, line) in enumerate(lines):
                header += ('#U%.2d %s'%(idx+3, line) + NEWLINE)
        header += '#N %d'%ncols + NEWLINE
        if ext == '.spec':
            header += ('#L ' + self.getGraphXTitle() + '  ' + '  '.join(legends) + NEWLINE)
        else:
            header += ('#L ' + self.getGraphXTitle() + '  ' + delim.join(legends) + NEWLINE)

        for fh in [filehandle, seperateFile]:
            if fh is not None:
                fh.write(NEWLINE)
                fh.write(header)
                for line in outArray:
                    tmp = delim.join(['%f'%num for num in line])
                    fh.write(tmp + NEWLINE)
                fh.write(NEWLINE)
                fh.close()
        
        # Emit saveOptionsSignal to save config file
        self.saveOptionsSignal.emit(splitext(filename)[0])
        if seperateFile is not None:
            self.saveOptionsSignal.emit(splitext(sepFileName)[0])
    
    def add(self):
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
        allCurves = self.getAllCurves()
        for (idx, (xVal,  yVal,  legend,  info)) in enumerate(allCurves):
            if 'selectionlegend' in info:
                newLegend = info['selectionlegend']
            elif 'operation' in info:
                newLegend = (str(operation) + ' ' + self.title)
            else:
                newLegend = (legend + ' ' + self.title)
            if idx == 0:
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
            if name == '$SEPARATOR':
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

    __colGroup   = 0
    __colLegend  = 1
    __colScanNo  = 2
    __colCounter = 3
    selectionModifiedSignal = qt.pyqtSignal()

    def __init__(self,  parent, groups = ['B','A','D'], color=True):
        qt.QTreeWidget.__init__(self,  parent)
        # Last identifier in groups is the ignore instruction
        self.groupList = groups
        self.actionList  = []
        self.contextMenu = qt.QMenu('Perform',  self)
        self.color = color
        self.colorDict = {
            groups[0] : qt.QBrush(qt.QColor(220, 220, 255)),
            groups[1] : qt.QBrush(qt.QColor(255, 210, 210)),
            '': qt.QBrush(qt.QColor(255, 255, 255))
        }

    def sizeHint(self):
        vscrollbar = self.verticalScrollBar()
        width = vscrollbar.width()
        for i in range(self.columnCount()):
            width += (2 + self.columnWidth(i))
        return qt.QSize( width, 20*22)

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
                print('getColum -- Selected column out of bounds')
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
                            print('getColum -- Conversion failed!')
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
        sel = self.getColumn(self.__colLegend, True)
        self.clear()
        self.setHeaderLabels(headerLabels)
        for item in items:
            treeItem = TreeWidgetItem(self,  item)
            if self.color:
                idx = str(treeItem.text(self.__colGroup))
                for i in range(self.columnCount()):
                    treeItem.setBackground(i, self.colorDict[idx])
            if treeItem.text(self.__colLegend) in sel:
                treeItem.setSelected(True)

    def setSelectionAs(self, idx):
        '''
        Sets the items currently selected to 
        the identifier given in idx.
        '''
        if idx not in self.groupList:
            raise ValueError('XMCDTreeWidget: invalid identifer \'%s\'' % idx)
        sel = self.selectedItems()
        if idx == self.groupList[-1]:
            # Last identifier in self.groupList
            # is the dummy identifier
            idx = ''
        for item in sel:
            item.setText(self.__colGroup, idx)
            if self.color:
                for i in range(self.columnCount()):
                    item.setBackground(i, self.colorDict[idx])
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
        # Try to sort for scanNo
        #self.sortItems(self.__colLegend, qt.Qt.AscendingOrder)
        self.sortItems(self.__colScanNo, qt.Qt.AscendingOrder)
        if not seq:
            seq, chk = qt.QInputDialog.\
                getText(None, 
                        'Sequence Dialog', 
                        'Valid identifiers are: ' + ', '.join(self.groupList),
                        qt.QLineEdit.Normal,
                        'Enter sequence')
        seq = str(seq).upper()
        if not chk:
            return
        for idx in seq:
            if idx not in self.groupList:
                invalidMsg = qt.QMessageBox(None)
                invalidMsg.setText('Invalid identifier. Try again.')
                invalidMsg.setStandardButtons(qt.QMessageBox.Ok)
                invalidMsg.exec_()
                return
        if len(sel) != len(seq):
            # Assume pattern and repeat
            seq = seq * (len(sel)//len(seq) + 1)
            #invalidMsg = qt.QMessageBox(None)
            #invalidMsg.setText('Sequence length does not match item count.')
            #invalidMsg.setStandardButtons(qt.QMessageBox.Ok)
            #invalidMsg.exec_()
        for (idx, item) in zip(seq, sel):
            if idx == self.groupList[-1]:
                idx = ''
            item.setText(self.__colGroup, idx)
            if self.color:
                for i in range(self.columnCount()):
                    item.setBackground(i, self.colorDict[idx])
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
        are the identifiers ('D', 'A', 'B') and 
        the values are (sorted) lists containing
        legends to which the respective identifier is
        assigned to.
        '''
        out = dict((group, []) for group in self.groupList)
        root = self.invisibleRootItem()
        for i in range(root.childCount()):
            item   = root.child(i)
            group  = str(item.text(0))
            legend = str(item.text(1))
            #nCols  = item.columnCount()
            #legend = str(item.text(nCols-1))
            if len(group) == 0:
                group = self.groupList[-1]
            out[group] += [legend]
        for value in out.values():
            value.sort()
        return out

class XMCDWidget(qt.QWidget):

    setSelectionSignal = qt.pyqtSignal(object, object)

    def __init__(self,  parent,
                        plotWindow,
                        beamline,
                        nSelectors = 5):
        """
        Input
        -----
        plotWindow : ScanWindow instance
            ScanWindow from which curves are passed for
            XLD/XMCD Analysis
        nSelectors : Int
            Number of columns show in the widget. Per default
            these are
            <Group> <Legend> <ScanNo> <Counter> <Motor 1> ... <Motor5>
        """
        qt.QWidget.__init__(self, parent)
        self.setWindowIcon(qt.QIcon(IconDict['peak']))
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
        self.analysisWindow = XMCDScanWindow(origin=plotWindow, 
                                              parent=None)
        self.optsWindow = XMCDOptions(self, self.motorNamesList)
                                              
        self.selectionDict = {'D': [],
                              'B': [],
                              'A': []}
        self.setSizePolicy(qt.QSizePolicy.MinimumExpanding, 
                           qt.QSizePolicy.Expanding)

        self.setWindowTitle("XLD/XMCD Analysis")
        
        buttonOptions = qt.QPushButton('Options', self)
        buttonOptions.setToolTip(
            'Set normalization and interpolation\n'
           +'method and motors shown')
                                
        updatePixmap  = qt.QPixmap(IconDict["reload"])
        buttonUpdate  = qt.QPushButton(
                            qt.QIcon(updatePixmap), '', self)
        buttonUpdate.setIconSize(qt.QSize(21,21))
        buttonUpdate.setToolTip(
            'Update curves in XMCD Analysis\n'
           +'by checking the plot window')

        self.list = XMCDTreeWidget(self)
        labels = ['Group', 'Legend', 'S#','Counter']+\
            (['']*nSelectors)
        ncols  = len(labels)
        self.list.setColumnCount(ncols)
        self.list.setHeaderLabels(labels)
        self.list.setSortingEnabled(True)
        self.list.setSelectionMode(
            qt.QAbstractItemView.ExtendedSelection)
        listContextMenu = XMCDMenu(None)
        listContextMenu.setActionList(
              [('Perform analysis', self.triggerXMCD),
               ('$SEPARATOR', None),
               ('Set as A', self.setAsA),
               ('Set as B', self.setAsB),
               ('Enter sequence', self.list.setSelectionToSequence),
               ('Remove selection', self.list.clearSelection),
               ('$SEPARATOR', None),
               ('Invert selection', self.list.invertSelection), 
               ('Remove curve(s)', self.removeCurve_)])
        self.list.setContextMenu(listContextMenu)
        self.expCBox = qt.QComboBox(self)
        self.expCBox.setToolTip('Select configuration of predefined\n'
                               +'experiment or configure new experiment')
        
        self.experimentsDict = {
            'Generic Dichorism': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': '',
                  'motor1': '',
                  'motor2': '',
                  'motor3': '',
                  'motor4': ''
            },
            'ID08: XMCD 9 Tesla Magnet': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': 'phaseD',
                  'motor1': 'magnet',
                  'motor2': '',
                  'motor3': '',
                  'motor4': ''
            },
            'ID08: XMCD 5 Tesla Magnet': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': 'PhaseD',
                  'motor1': 'oxPS',
                  'motor2': '',
                  'motor3': '',
                  'motor4': ''
            },
            'ID08: XLD 5 Tesla Magnet': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': 'PhaseD',
                  'motor1': '',                    
                  'motor2': '',                     
                  'motor3': '',                    
                  'motor4': ''                    
            },
            'ID08: XLD 9 Tesla Magnet': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': 'phaseD',
                  'motor1': '',                  
                  'motor2': '',                     
                  'motor3': '',                     
                  'motor4': ''                     
            },
            'ID12: XMCD (Flipper)': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': 'BRUKER',
                  'motor1': 'OXFORD',
                  'motor2': 'CRYO',
                  'motor3': '',
                  'motor4': ''
            },
            'ID12: XMCD': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': 'Phase',
                  'motor1': 'PhaseA',
                  'motor2': 'BRUKER',
                  'motor3': 'OXFORD',
                  'motor4': 'CRYO'
            },
            'ID12: XLD (quater wave plate)': {
                  'xrange': 0,
                  'normalization': 0,
                  'normalizationMethod': 'offsetAndArea',
                  'motor0': '',
                  'motor1': '',
                  'motor2': '',
                  'motor3': '',
                  'motor4': ''
            }
        }
        self.expCBox.addItems(
                        ['Generic Dichorism',
                         'ID08: XLD 9 Tesla Magnet',
                         'ID08: XLD 5 Tesla Magnet',
                         'ID08: XMCD 9 Tesla Magnet',
                         'ID08: XMCD 5 Tesla Magnet',
                         'ID12: XLD (quater wave plate)',
                         'ID12: XMCD (Flipper)',
                         'ID12: XMCD',
                         'Add new configuration'])
        self.expCBox.insertSeparator(len(self.experimentsDict))
        
        topLayout  = qt.QHBoxLayout()
        topLayout.addWidget(buttonUpdate)
        topLayout.addWidget(buttonOptions)
        topLayout.addWidget(qt.HorizontalSpacer(self))
        topLayout.addWidget(self.expCBox)

        leftLayout = qt.QGridLayout()
        leftLayout.setContentsMargins(1, 1, 1, 1)
        leftLayout.setSpacing(2)
        leftLayout.addLayout(topLayout, 0, 0)
        leftLayout.addWidget(self.list, 1, 0)
        leftWidget = qt.QWidget(self)
        leftWidget.setLayout(leftLayout)
        
        self.analysisWindow.setSizePolicy(qt.QSizePolicy.Minimum, qt.QSizePolicy.Minimum)
        #self.splitter = qt.QSplitter(qt.Qt.Horizontal, self)
        self.splitter = qt.QSplitter(qt.Qt.Vertical, self)
        self.splitter.addWidget(leftWidget)
        self.splitter.addWidget(self.analysisWindow)
        stretch = int(leftWidget.width())
        # If window size changes, only the scan window size changes
        self.splitter.setStretchFactor(
                    self.splitter.indexOf(self.analysisWindow),1)
        
        mainLayout = qt.QVBoxLayout()
        mainLayout.setContentsMargins(0,0,0,0)
        mainLayout.addWidget(self.splitter)
        self.setLayout(mainLayout)

        # Connects
        self.expCBox.currentIndexChanged['QString'].connect(self.updateTree)
        self.expCBox.currentIndexChanged['QString'].connect(self.selectExperiment)
        self.list.selectionModifiedSignal.connect(self.updateSelectionDict)
        self.setSelectionSignal.connect(self.analysisWindow.processSelection)
        self.analysisWindow.saveOptionsSignal.connect(self.optsWindow.saveOptions)
        self.optsWindow.accepted[()].connect(self.updateTree)
        buttonUpdate.clicked.connect(self.updatePlots)
        buttonOptions.clicked.connect(self.showOptionsWindow)

        self.updateTree()
        self.list.sortByColumn(1, qt.Qt.AscendingOrder)

    def sizeHint(self):
        return self.list.sizeHint() + self.analysisWindow.sizeHint()

    def addExperiment(self):
        exp, chk = qt.QInputDialog.\
                        getText(self,
                                'Configure new experiment',
                                'Enter experiment title',
                                qt.QLineEdit.Normal, 
                                'ID00: <Title>')
        if chk and (not exp.isEmpty()):
            exp = str(exp)
            opts = XMCDOptions(self, self.motorNamesList, False)
            if opts.exec_():
                self.experimentsDict[exp] = opts.getOptions()
                cBox = self.expCBox
                new = [cBox.itemText(i) for i in range(cBox.count())][0:-2]
                new += [exp]
                new.append('Add new configuration')
                cBox.clear()
                cBox.addItems(new)
                cBox.insertSeparator(len(new)-1)
                idx = cBox.findText([exp][0])
                if idx < 0:
                    cBox.setCurrentIndex(0)
                else:
                    cBox.setCurrentIndex(idx)
            opts.destroy()
            idx = self.expCBox.findText(exp)
            if idx < 0:
                idx = 0
            self.expCBox.setCurrentIndex(idx)
        
    def showOptionsWindow(self):
        if self.optsWindow.exec_():
            options = self.optsWindow.getOptions()
            self.analysisWindow.processOptions(options)

# Implement new assignment routines here BEGIN
    def selectExperiment(self, exp):
        exp = str(exp)
        if exp == 'Add new configuration':
            self.addExperiment()
            self.updateTree()
        elif exp in self.experimentsDict:
            try:
                # Sets motors 0 to 4 in optsWindow
                self.optsWindow.setOptions(self.experimentsDict[exp])
            except ValueError:
                self.optsWindow.setOptions(
                        self.experimentsDict['Generic Dichorism'])
                return
            # Get motor values from tree
            self.updateTree()
            values0 = numpy.array(
                        self.list.getColumn(4, convertType=float))
            values1 = numpy.array(
                        self.list.getColumn(5, convertType=float))
            values2 = numpy.array(
                        self.list.getColumn(6, convertType=float))
            values3 = numpy.array(
                        self.list.getColumn(7, convertType=float))
            values4 = numpy.array(
                        self.list.getColumn(8, convertType=float))
            # Determine p/m selection
            if exp.startswith('ID08: XLD'):
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
            elif exp.startswith('ID12: XLD (quater wave plate)'):
                # Extract counters from third column
                counters = self.list.getColumn(3, convertType=str)
                polarization = []
                for counter in counters:
                    # Relevant counters Ihor, Iver resp. Ihor0, Iver0, etc.
                    if 'hor' in counter:
                        pol = -1.
                    elif 'ver' in counter:
                        pol =  1.
                    else:
                        pol = float('nan')
                    polarization += [pol]
                values = numpy.asarray(polarization, dtype=float)
                vpivot = 0.
            elif exp.startswith('ID12: XMCD (Flipper)'):
                # Extract counters from third column
                counters = self.list.getColumn(1, convertType=str)
                polarization = []
                for counter in counters:
                    # Relevant counters: Fminus/Fplus resp. Rminus/Rplus
                    if 'minus' in counter:
                        pol = 1.
                    elif 'plus' in counter:
                        pol = -1.
                    else:
                        pol = float('nan')
                    polarization += [pol]
                magnets = values0 + values1 + values2
                values = numpy.asarray(polarization, dtype=float)*\
                            magnets
                vpivot = 0.
            elif exp.startswith('ID12: XMCD'):
                # Sum over phases..
                polarization = values0 + values1 
                # ..and magnets
                magnets = values2 + values3 + values4
                signMagnets = numpy.sign(magnets)
                if numpy.all(signMagnets==0.):
                    values = polarization
                else:
                    values = numpy.sign(polarization)*\
                                numpy.sign(magnets)
                vpivot = 0.
            else:
                values = numpy.array([float('NaN')]*len(self.legendList))
                vpivot = 0.
            # Sequence is generate according to values and vpivot
            seq = ''
            for x in values:
                if str(x) == 'nan':
                    seq += 'D'
                elif x>vpivot:
                    seq += 'A'
                else:
                    seq += 'B'
            self.list.setSelectionToSequence(seq)
# Implement new assignment routines here END

    def triggerXMCD(self):
        groupA = self.selectionDict['A']
        groupB = self.selectionDict['B']
        self.analysisWindow.processSelection(groupA, groupB)

    def removeCurve_(self):
        sel = self.list.getColumn(1,
                                  selectedOnly=True,
                                  convertType=str)
        for legend in sel:
            self.plotWindow.removeCurve(legend)
            for selection in self.selectionDict.values():
                if legend in selection:
                    selection.remove(legend)
            # Remove from XMCDScanWindow.curvesDict
            if legend in self.analysisWindow.curvesDict.keys():
                del(self.analysisWindow.curvesDict[legend])
            # Remove from XMCDScanWindow.selectionDict
            for selection in self.analysisWindow.selectionDict.values():
                if legend in selection:
                    selection.remove(legend)
        self.updatePlots()

    def updateSelectionDict(self):
        # Get selDict from self.list. It consists of tree items:
        # {GROUP0: LIST_OF_LEGENDS_IN_GROUP0, 
        #  GROUP1: LIST_OF_LEGENDS_IN_GROUP1,
        #  GROUP2: LIST_OF_LEGENDS_IN_GROUP2}
        selDict = self.list.getSelection()
        # self.selectionDict -> Uses ScanNumbers instead of legends...
        newDict = {}
        for (idx, selList) in selDict.items():
            if idx not in newDict.keys():
                newDict[idx] = []
            for legend in selList:
                newDict[idx] += [legend]
        self.selectionDict = newDict
        self.setSelectionSignal.emit(self.selectionDict['A'],
                                     self.selectionDict['B'])

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
        if experiment != 'Generic Dichorism':
            self.selectExperiment(experiment)
        return

    def updateTree(self):
        mList  = self.optsWindow.getMotors()
        labels = ["Group",'Legend','S#','Counter'] + mList
        items  = []
        for i in range(len(self.legendList)):
            # Loop through rows
            # Each row is represented by QStringList
            legend = self.legendList[i]
            values = self.motorsList[i]
            info = self.infoList[i]
            selection = ''
            # Determine Group from selectionDict
            for (idx, v) in self.selectionDict.items():
                if (legend in v) and (idx != 'D'):
                    selection = idx
                    break
            # Add filename, scanNo, counter
            #sourceName = info.get('SourceName','')
            #if isinstance(sourceName,list):
            #    filename = basename(sourceName[0])
            #else:
            #    filename = basename(sourceName)
            filename = legend
            scanNo = info.get('Key','')
            counter = info.get('ylabel',None)
            if counter is None:
                selDict = info.get('selection',{})
                if len(selDict) == 0:
                    counter = ''
                else:
                    # When do multiple selections occur?
                    try:
                        yIdx = selDict['y'][0]
                        cntList = selDict['cnt_list']
                        counter = cntList[yIdx]
                    except Exception:
                        counter = ''
            tmp = qt.QStringList([selection, filename, scanNo, counter])
            # Determine value for each motor
            for m in mList:
                if len(m) == 0:
                    tmp.append('')
                else:
                    tmp.append(str(values.get(m, '---')))
            items.append(tmp)
        self.list.build(items,  labels)
        for idx in range(self.list.columnCount()):
            self.list.resizeColumnToContents(idx)

    def setAsA(self):
        self.list.setSelectionAs('A')

    def setAsB(self):
        self.list.setSelectionAs('B')

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
        '''
        Curves retrieved from the main plot window using the
        Plugin1DBase getActiveCurve() resp. getAllCurves()
        member functions are tuple resp. a list of tuples
        containing x-data, y-data, legend and the info dictionary.
        
        _setLists splits these tuples into lists, thus setting
        the attributes
        
            self.legendList
            self.infoList
            self.motorsList
        '''
        if self.plotWindow is not None:
            curves = self.plotWindow.getAllCurves()
        else:
            if DEBUG:
                print('_setLists -- Set self.plotWindow before calling self._setLists')
            return
        # nCurves = len(curves)
        self.legendList = [leg for (xvals, yvals,  leg,  info) in curves] 
        self.infoList   = [info for (xvals, yvals,  leg,  info) in curves]
        # Try to recover the scan number from the legend, if not set
        #for ddict in self.infoList:
        #    key = ddict.get('Key','')
        #    if len(key)== 0:
        #        selectionlegend = ddict['selectionlegend']
        #        match = regexpSearch(r'(?<= )\d{1,5}\.\d{1}',selectionlegend)
        #        if match:
        #            scanNo = match.group(0)
        #            ddict['Key'] = scanNo
        self.motorsList = self._convertInfoDictionary(self.infoList)

class XMCDFileDialog(qt.QFileDialog):
    def __init__(self, parent, caption, directory, filter):
        qt.QFileDialog.__init__(self, parent, caption, directory, filter)
        
        saveOptsGB = qt.QGroupBox('Save options', self)
        self.appendBox = qt.QCheckBox('Append to existing file', self)
        self.commentBox = qt.QTextEdit('Enter comment', self)
        
        mainLayout = self.layout()
        optsLayout = qt.QGridLayout()
        optsLayout.addWidget(self.appendBox, 0, 0)
        optsLayout.addWidget(self.commentBox, 1, 0)
        saveOptsGB.setLayout(optsLayout)
        mainLayout.addWidget(saveOptsGB, 4, 0, 1, 3)
        
        self.appendBox.stateChanged.connect(self.appendChecked)

    def appendChecked(self, state):
        if state == qt.Qt.Unchecked:
            self.setConfirmOverwrite(True)
            self.setFileMode(qt.QFileDialog.AnyFile)
        else:
            self.setConfirmOverwrite(False)
            self.setFileMode(qt.QFileDialog.ExistingFile)

def getSaveFileName(parent, caption, directory, filter):
    dial = XMCDFileDialog(parent, caption, directory, filter)
    dial.setAcceptMode(qt.QFileDialog.AcceptSave)
    append = None
    comment = None
    files = []
    if dial.exec_():
        append  = dial.appendBox.isChecked()
        comment = str(dial.commentBox.toPlainText())
        if comment == 'Enter comment':
            comment = ''
        files = [qt.safe_str(fn) for fn in dial.selectedFiles()]
    return (files, append, comment)

def main():
    app = qt.QApplication([])
    
    # Create dummy ScanWindow
    swin = sw.ScanWindow()
    info0 = {'xlabel': 'foo',
             'ylabel': 'arb',
             'MotorNames': 'oxPS PhaseA Phase BRUKER CRYO OXFORD', 
             'MotorValues': '1 -6.27247094 -3.11222732 6.34150808 -34.75892563 21.99607165'}
    info1 = {'MotorNames': 'PhaseD oxPS PhaseA Phase BRUKER CRYO OXFORD',
             'MotorValues': '0.470746882688 0.25876374531 -0.18515967 -28.31216591 18.54513221 -28.09735532 -26.78833172'}
    info2 = {'MotorNames': 'PhaseD oxPS PhaseA Phase BRUKER CRYO OXFORD',
             'MotorValues': '-9.45353059 -25.37448851 24.37665651 18.88048044 -0.26018745 2 0.901968648111 '}
    x = numpy.arange(100.,1100.)
    y0 =  10*x + 10000.*numpy.exp(-0.5*(x-500)**2/400) + 1500*numpy.random.random(1000.)
    y1 =  10*x + 10000.*numpy.exp(-0.5*(x-600)**2/400) + 1500*numpy.random.random(1000.)
    y2 =  10*x + 10000.*numpy.exp(-0.5*(x-400)**2/400) + 1500*numpy.random.random(1000.)
    
    swin.newCurve(x, y2, legend="Curve2", xlabel='ene_st2', ylabel='Ihor', info=info2, replot=False, replace=False)
    swin.newCurve(x, y0, legend="Curve0", xlabel='ene_st0', ylabel='Iver', info=info0, replot=False, replace=False)
    swin.newCurve(x, y1, legend="Curve1", xlabel='ene_st1', ylabel='Ihor', info=info1, replot=False, replace=False)
    
    # info['Key'] is overwritten when using newCurve
    swin.dataObjectsDict['Curve2 Ihor'].info['Key'] = '1.1'
    swin.dataObjectsDict['Curve0 Iver'].info['Key'] = '34.1'
    swin.dataObjectsDict['Curve1 Ihor'].info['Key'] = '123.1'

    w = XMCDWidget(None, swin, 'ID08', nSelectors = 5)
    w.show()
    app.exec_()

if __name__ == '__main__':
    main()
