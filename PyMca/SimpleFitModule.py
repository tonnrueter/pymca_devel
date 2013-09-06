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
import sys
import os
import numpy
import copy
import types
from PyMca import Gefit
from PyMca import SpecfitFuns

DEBUG = 0

class SimpleFit(object):
    def __init__(self):
        #no data available by default
        self._x0 = None
        self._y0 = None
        
        #get default configuration
        self.getDefaultConfiguration()

        #the list and dictionary of defined functions
        self._functionList = []
        self._functionDict = {}

        #the current fit function
        self._stripFunction = None

    def getDefaultConfiguration(self):
        self._fitConfiguration = {}
        self._fitConfiguration['fit'] = {}
        self._fitConfiguration['fit']['fit_function'] = "None"
        self._fitConfiguration['fit']['function_flag'] = 1
        self._fitConfiguration['fit']['background_function'] = "None"
        self._fitConfiguration['fit']['background_flag'] = 1
        self._fitConfiguration['fit']['strip_function'] = "Strip"
        self._fitConfiguration['fit']['stripalgorithm'] = 0
        self._fitConfiguration['fit']['strip_flag'] = 1        
        self._fitConfiguration['fit']['fit_algorithm'] = "Levenberg-Marquardt"
        self._fitConfiguration['fit']['weight'] = "NO Weight"
        self._fitConfiguration['fit']['maximum_fit_iterations'] = 10
        self._fitConfiguration['fit']['background_estimation_policy'] = "Estimate always"
        self._fitConfiguration['fit']['function_estimation_policy'] = "Estimate always"
        self._fitConfiguration['fit']['minimum_delta_chi'] = 0.0010
        self._fitConfiguration['fit']['use_limits'] = 0
        self._fitConfiguration['fit']['xmin'] =    0.
        self._fitConfiguration['fit']['xmax'] = 1000.
        self._fitConfiguration['fit']['functions'] = []
        #strip/snip background configuration related
        self._fitConfiguration['fit']['stripanchorsflag'] = 0
        self._fitConfiguration['fit']['stripanchorslist'] = []
        self._fitConfiguration['fit']['stripfilterwidth'] = 1
        self._fitConfiguration['fit']['snipwidth'] = 10
        self._fitConfiguration['fit']['stripwidth'] = 4
        self._fitConfiguration['fit']['stripiterations'] = 5000
        self._fitConfiguration['fit']['stripconstant'] = 1.0
        self._fitConfiguration['functions'] = {}
        
    def configure(self, ddict=None):
        if ddict is None:
            return self.getConfiguration()
        else:
            return self.setConfiguration(ddict)

    def setConfiguration(self, ddict, try_import=False):
        oldConfig = self.getConfiguration()
        if ddict is None:
            return oldConfig
        if 'fit' in ddict:
            givenKeys = ddict['fit'].keys()
            for key in self._fitConfiguration['fit'].keys():
                if key in givenKeys:
                    self._fitConfiguration['fit'][key] = ddict['fit'][key]
            for key in ddict.keys():
                if key in ['fit', 'functions']:
                    continue
                self._fitConfiguration[key] = ddict[key]

        if 'functions' in ddict:
            functionNames = ddict['functions'].keys()
            for fName in functionNames:
                if fName not in self._fitConfiguration['functions'].keys():
                    if try_import:
                        ffile = ddict['functions'][fName].get('file', None)
                        if ffile is not None:
                            self.importFunctions(ffile)
                    else:
                        print("WARNING:Function %s not among defined functions" % fName)
                        continue
                self._fitConfiguration['functions'][fName]['configuration']=\
                        ddict['functions'][fName]['configuration']
                configureMethod = self._fitConfiguration['functions'][fName]\
                                  ['configure']
                if configureMethod is not None:
                    configureMethod(ddict['functions'][fName]['configuration'])

        #if data are present, update strip background                    
        if (self._x0 is None) or (self._y0 is None):
            return
        if (oldConfig['fit']['xmin'] != self._fitConfiguration['fit']['xmin']) or\
           (oldConfig['fit']['xmax'] != self._fitConfiguration['fit']['xmax']):
            if DEBUG:
                print("SETTING DATA AGAIN")
            self.setData(self._x0, self._y0,
                         xmin=self._fitConfiguration['fit']['xmin'],
                         xmax=self._fitConfiguration['fit']['xmax'])
            return

        for key in ['strip_flag', 'stripanchorsflag', 'stripalgorithm',
                    'stripwidth', 'stripiterations', 'stripconstant']:
            if oldConfig['fit'][key] != self._fitConfiguration['fit'][key]:
                if DEBUG:
                    print("RECALCULATING STRIP")
                self._getStripBackground()
                break
            if key == 'stripanchorsflag':
                if len(oldConfig['fit']['stripanchorslist']) !=\
                   len(self._fitConfiguration['fit']['stripanchorslist']):
                    if DEBUG:
                        print("ANCHORS CHANGE, RECALCULATING STRIP")
                    self._getStripBackground()
                    break

    def getConfiguration(self):
        ddict = {}
        for key in self._fitConfiguration.keys():
            if key == 'functions':
                continue
            ddict[key] = copy.deepcopy(self._fitConfiguration[key])

        ddict['functions'] = {}
        for key in self._fitConfiguration['functions'].keys():
            ddict['functions'][key] = {}
            ddict['functions'][key]['configuration'] = copy.deepcopy(\
                self._fitConfiguration['functions'][key]['configuration'])
            configureMethod = self._fitConfiguration['functions']\
                                                      [key]['configure']
            if configureMethod is not None:
                currentFunctionConfig = configureMethod()
                for newKey in currentFunctionConfig.keys():
                    if newKey not in ['estimation']:
                        ddict['functions'][key]['configuration'][newKey] = currentFunctionConfig[newKey]
            parameters = self._fitConfiguration['functions'][key]['parameters']
            ddict['functions'][key]['parameters'] =  parameters
            widget = self._fitConfiguration['functions'][key]['widget']
            ddict['functions'][key]['widget'] = widget
            fname = self._fitConfiguration['functions'][key]['file']
            ddict['functions'][key]['file'] = fname
        return ddict

    def setData(self, x, y, sigma=None, xmin=None, xmax=None, **kw):
        # make sure last fit result is not used
        self._fitResult = None
        idx = numpy.argsort(x)
        if sigma is not None:
            self._sigma = sigma[idx]
        else:
            self._sigma = None
        self._y0 = y[idx]
        self._x0 = x[idx]
        if sigma is not None:
            self._sigma0 = sigma[idx]
        xmin, xmax = self._getLimits(self._x0, xmin, xmax)
        idx = (self._x0 >= xmin) & (self._x0 <= xmax)
        self._x = self._x0[idx]
        self._y = self._y0[idx]
        self._fitConfiguration['fit']['xmin'] = xmin * 1.0
        self._fitConfiguration['fit']['xmax'] = xmax * 1.0
        if sigma is not None:
            self._sigma = self._sigma0[idx]
        if DEBUG:
            print("TODO: Make sure we have something to fit")
        #get strip/SNIP background
        self._z = self._getStripBackground()

    def importFunctions(self, modname):
        #modname can be a module or a file
        if type(modname) == types.ModuleType:
            newfun = modname
        elif os.path.exists(modname):
            sys.path.append(os.path.dirname(modname))
            f=os.path.basename(os.path.splitext(modname)[0])
            newfun=__import__(f)
        else:
            raise ValueError("Cannot interprete/find %s" % modname)
        
        theory = newfun.THEORY
        function=newfun.FUNCTION
        parameters = newfun.PARAMETERS
        try:
            estimate=newfun.ESTIMATE
        except:
            estimate=None        
        try:
            derivative=newfun.DERIVATIVE
        except:
            derivative=None
        try:
            configure=newfun.CONFIGURE
        except:
            configure=None
        try:
            widget=newfun.WIDGET
        except:
            widget=None

        basename = os.path.basename(newfun.__file__)

        for i in range(len(theory)):
            ddict = {}
            functionName = theory[i]
            ddict['function'] = function[i]
            ddict['parameters'] = parameters[i]
            ddict['default_parameters'] = None
            ddict['estimate']   = None
            ddict['derivative'] = None
            ddict['configure']  = None
            ddict['widget']     = None
            ddict['file']       = newfun.__file__
            ddict['configuration'] = {} 
            if estimate is not None:
                ddict['estimate'] = estimate[i]
            if derivative is not None:
                ddict['derivative'] = derivative[i]
            if configure is not None:
                ddict['configure'] = configure[i]
                if ddict['configure'] is not None:
                    ddict['configuration'] = configure[i]()
                if ddict['estimate'] is None:
                    ddict['configuration']['estimation'] = None
            if widget is not None:
                ddict['widget'] = widget[i]
            self._fitConfiguration['functions'][functionName] = ddict
            self._fitConfiguration['fit']['functions'].append(functionName)

    def setFitFunction(self, name):
        if name in [None, "None", "NONE"]:
            self._fitConfiguration['fit']['fit_function'] = "None"
            return
        self._fitFunctionConfigured = False
        if name not in self._fitConfiguration['fit']['functions']:
            txt = "Function %s not among defined functions"  % name
            raise KeyError(txt)
        self._fitConfiguration['fit']['fit_function'] = name

    def getFitFunction(self):
        return "%s" % self._fitConfiguration['fit']['fit_function']

    def setBackgroundFunction(self, name):
        if name in [None, "None", "NONE"]:
            self._fitConfiguration['fit']['background_function'] = "None"
            return
        self._backgroundFunctionConfigured = False
        if name not in self._fitConfiguration['fit']['functions']:
            txt = "Function %s not among defined functions"  % name
            raise KeyError(txt)
        self._fitConfiguration['fit']['background_function'] = name

    def getBackgroundFunction(self):
        return "%s" % self._fitConfiguration['fit']['background_function']
    
    def _getLimits(self, x, xmin, xmax):
        if self._fitConfiguration['fit']['use_limits']:
            xmin = self._fitConfiguration['fit']['xmin']
            xmax = self._fitConfiguration['fit']['xmax']
            return xmin, xmax
        if xmin is None:
            xmin = x[0]
        if xmax is None:
            xmax = x[-1]
        return xmin, xmax
            
    def _getStripBackground(self, x=None, y=None):
        #this makes the assumption x are equally spaced
        #and I should build a spline if that is not the case
        #but I do not want to put a dependency on SciPy
        if y is not None:
            ywork = y
        else:
            ywork = self._y

        if x is not None:
            xwork = x
        else:
            xwork = self._x

        n=len(xwork)
        #loop for anchors
        anchorslist = []
        if self._fitConfiguration['fit']['stripanchorsflag']:
            if self._fitConfiguration['fit']['stripanchorslist'] is not None:
                oldShape = xwork.shape
                ravelled = xwork
                ravelled.shape = -1
                for channel in self._fitConfiguration['fit']['stripanchorslist']:
                    if channel <= ravelled[0]:continue
                    index = numpy.nonzero(ravelled >= channel)
                    if len(index):
                        index = min(index)
                        if index > 0:
                            anchorslist.append(index)
                ravelled.shape = oldShape

        #work with smoothed data
        ysmooth = self._getSmooth(xwork, ywork)
        
        #SNIP algorithm
        if self._fitConfiguration['fit']['stripalgorithm'] in ["SNIP", 1]:
            if DEBUG:
                print("CALCULATING SNIP")
            if len(anchorslist) == 0:
                anchorslist = [0, len(ysmooth)-1]
            anchorslist.sort()
            result = 0.0 * ysmooth
            lastAnchor = 0
            width = self._fitConfiguration['fit']['snipwidth']
            for anchor in anchorslist:
                if (anchor > lastAnchor) and (anchor < len(ysmooth)):
                    result[lastAnchor:anchor] =\
                            SpecfitFuns.snip1d(ysmooth[lastAnchor:anchor], width, 0)
                    lastAnchor = anchor
            if lastAnchor < len(ysmooth):                
                result[lastAnchor:] =\
                        SpecfitFuns.snip1d(ysmooth[lastAnchor:], width, 0)
            return result
        
        #strip background
        niter = self._fitConfiguration['fit']['stripiterations']
        if niter > 0:
            if DEBUG:
                print("CALCULATING STRIP")
                print("iterations = ", niter)
                print("constant   = ",
                      self._fitConfiguration['fit']['stripconstant'])
                print("width      = ",
                      self._fitConfiguration['fit']['stripwidth'])
                print("anchors    = ", anchorslist)
            result=SpecfitFuns.subac(ysmooth,
                                  self._fitConfiguration['fit']['stripconstant'],
                                  niter,
                                  self._fitConfiguration['fit']['stripwidth'],
                                  anchorslist)
            if niter > 1000:
                #make sure to get something smooth
                result = SpecfitFuns.subac(result,
                                  self._fitConfiguration['fit']['stripconstant'],
                                  500,1,
                                  anchorslist)
            else:
                #make sure to get something smooth but with less than
                #500 iterations
                result = SpecfitFuns.subac(result,
                                  self._fitConfiguration['fit']['stripconstant'],
                                  int(self._fitConfiguration['fit']['stripwidth']*2),
                                  1,
                                  anchorslist)
        else:
            if DEBUG:
                print("NO STRIP, NO SNIP")
            result     = numpy.zeros(ysmooth.shape, numpy.float) + min(ysmooth)

        return result

    def _getSmooth(self, x, y): 
        f=[0.25,0.5,0.25]
        try:
            if hasattr(y, "shape"):
                if len(y.shape) > 1:
                    result=SpecfitFuns.SavitskyGolay(numpy.ravel(y).astype(numpy.float), 
                                    self._fitConfiguration['fit']['stripfilterwidth'])
                else:                                
                    result=SpecfitFuns.SavitskyGolay(numpy.array(y).astype(numpy.float), 
                                    self._fitConfiguration['fit']['stripfilterwidth'])
            else:
                result=SpecfitFuns.SavitskyGolay(numpy.array(y).astype(numpy.float),
                                    self._fitConfiguration['fit']['stripfilterwidth'])
        except:
            err = sys.exc_info()[1]
            raise ValueError("Unsuccessful Savitsky-Golay smoothing: %s" % err)
            result=numpy.array(y).astype(numpy.float)
        if len(result) > 1:
            result[1:-1]=numpy.convolve(result,f,mode=0)
            result[0]=0.5*(result[0]+result[1])
            result[-1]=0.5*(result[-1]+result[-2])
        return result

    def fit(self):
        if self._y0 is None:
            self._setStatus("No data to be fitted")
            return
        self.estimate()
        self.startFit()
        return self.getResult()

    def estimate(self):
        self._fitResult = None
        if self._y0 is None:
            self._setStatus("No data to be fitted")
            return
        self._setStatus("Estimate started")
        backgroundDict  = {'parameters':[]}
        fitFunctionDict = {'parameters':[]}
        backgroundParameters, backgroundConstraints = [], [[],[],[]]
        backgroundFunction = self.getBackgroundFunction() 
        if self._fitConfiguration['fit']['background_flag']:
            if backgroundFunction not in [None, "None", "NONE"]:
                backgroundParameters, backgroundConstraints =\
                                      self.estimateBackground()
                backgroundDict = self._fitConfiguration['functions']\
                              [backgroundFunction]
        self._setStatus("Background estimation finished")
        functionParameters, functionConstraints = [], [[],[],[]]
        fitFunction = self._fitConfiguration['fit']['fit_function']
        if self._fitConfiguration['fit']['function_flag']:
            if fitFunction not in [None, "None", "NONE"]:
                functionParameters, functionConstraints=\
                                    self.estimateFunction()
                fitFunctionDict = self._fitConfiguration['functions']\
                                      [fitFunction]
        if DEBUG:
            print("ESTIMATION parameters  = ",functionParameters)
            print("ESTIMATION constraints = ",functionConstraints)
        self._setStatus("Fit function estimation finished")

        #estimations are made
        #Check if there can be conflicts between parameter names
        #because they can have same names in the background and
        #in the fit function
        conflict = False
        for parname in backgroundDict['parameters']:
            if parname in fitFunctionDict['parameters']:
                conflict = True
                break

        #build the parameter names
        self.final_theory=[]
        nBasePar   = len(backgroundDict['parameters'])
        nActualPar = len(backgroundParameters)
        self.__nBackgroundParameters = nActualPar
        if nActualPar:
            for i in range(nActualPar):
                parname = backgroundDict['parameters'][i%nBasePar]
                if conflict:
                    parname = "Bkg_"+parname
                if nBasePar < nActualPar:
                    parname = parname + ("_%d" % (1+int(i/nBasePar)))
                self.final_theory.append(parname)

        nBasePar   = len(fitFunctionDict['parameters'])
        nActualPar = len(functionParameters)
        if nActualPar:
            for i in range(nActualPar):
                parname = fitFunctionDict['parameters'][i%nBasePar]
                if nBasePar < nActualPar:
                    parname = parname + ("_%d" % (1+int(i/nBasePar)))
                self.final_theory.append(parname)

        CONS=['FREE',
            'POSITIVE',
            'QUOTED',
            'FIXED',
            'FACTOR',
            'DELTA',
            'SUM',
            'IGNORE']
        self.paramlist=[]
        param          = self.final_theory
        j=0
        i=0
        k=0
        xmin=self._x.min()
        xmax=self._x.max()
        #print "xmin = ",xmin,"xmax = ",xmax
        for pname in self.final_theory:
            if i < len(backgroundParameters):
                self.paramlist.append({'name':pname,
                        'estimation':backgroundParameters[i],
                        'group':0,
                        'code':CONS[int(backgroundConstraints[0][i])],
                        'cons1':backgroundConstraints[1][i],
                        'cons2':backgroundConstraints[2][i],
                        'fitresult':0.0,
                        'sigma':0.0,
                        'xmin':xmin,
                        'xmax':xmax})
                i=i+1
            else:
                if (j % len(fitFunctionDict['parameters'])) == 0:
                    k=k+1
                if (CONS[int(functionConstraints[0][j])] == "FACTOR") or \
                   (CONS[int(functionConstraints[0][j])] == "DELTA"):
                        functionConstraints[1][j] = functionConstraints[1][j] +\
                                                    len(backgroundParameters)
                self.paramlist.append({'name':pname,
                       'estimation':functionParameters[j],
                       'group':k,
                       'code':CONS[int(functionConstraints[0][j])],
                       'cons1':functionConstraints[1][j],
                       'cons2':functionConstraints[2][j],
                       'fitresult':0.0,
                       'sigma':0.0,
                       'xmin':xmin,
                       'xmax':xmax})
                j=j+1
        self._setStatus("Estimate finished")
        return self.paramlist

    def _setStatus(self, status):
        self.__status = status

    def getStatus(self):
        return self.__status

    def estimateBackground(self):
        fname = self.getBackgroundFunction()
        if fname is None:
            return [],[[],[],[]]
        ddict = self._fitConfiguration['functions'][fname]
        estimateFunction = ddict['estimate']
        
        if estimateFunction is None:
            parameters = []
            constraints = [[],[],[]]
            if ddict['configuration']['estimation'] is not None:
                estimation = ddict['configuration']['estimation']
                defaultPar = estimation['parameters']
                for parameter in defaultPar:
                    parameters.append(estimation[parameter]['estimation'])
                    constraints[0].append(estimation[parameter]['code'])
                    constraints[1].append(estimation[parameter]['cons1'])
                    constraints[2].append(estimation[parameter]['cons2'])
            else:
                defaultPar = ddict['parameters']
                for parameter in defaultPar:
                    parameters.append(0.0)
                    constraints[0].append(0)
                    constraints[1].append(0)
                    constraints[2].append(0)
            return parameters, constraints
        parameters, constraints = estimateFunction(self._x, self._y, self._z)
        return parameters, constraints

    def estimateFunction(self):
        self._z = self._getStripBackground()
        fname = self.getFitFunction()
        if fname is None:
            return [],[[],[],[]]
        ddict = self._fitConfiguration['functions'][fname]
        estimateFunction = ddict['estimate']
        if estimateFunction is None:
            parameters = []
            constraints = [[],[],[]]
            if ddict['configuration']['estimation'] is not None:
                estimation = ddict['configuration']['estimation']
                defaultPar = estimation['parameters']
                for parameter in defaultPar:
                    parameters.append(estimation[parameter]['estimation'])
                    constraints[0].append(estimation[parameter]['code'])
                    constraints[1].append(estimation[parameter]['cons1'])
                    constraints[2].append(estimation[parameter]['cons2'])
            else:
                defaultPar = ddict['parameters']
                for parameter in defaultPar:
                    parameters.append(0.0)
                    constraints[0].append(0)
                    constraints[1].append(0)
                    constraints[2].append(0)
            return parameters, constraints
        parameters, constraints = estimateFunction(self._x * 1,
                                                   self._y * 1,
                                                   self._z * 1)
        return parameters, constraints

    def startFit(self):
        if self._y0 is None:
            self._setStatus("No data to be fitted")
            return
        self._setStatus("Fit started")
        param_list = self.final_theory
        length      = len(param_list)
        param_val   = []
        param_constrains   = [[],[],[]]
        flagconstrained=0
        for param in self.paramlist:
            #print param['name'],param['group'],param['estimation']
            param_val.append(param['estimation'])
            if (param['code'] != 'FREE') & (param['code'] != 0) & \
               (param['code'] != 0.0) :
                flagconstrained=1
            param_constrains [0].append(param['code'])
            param_constrains [1].append(param['cons1'])
            param_constrains [2].append(param['cons2'])

        #weight handling
        if self._fitConfiguration['fit']['weight'] in ["NO Weight", 0]:
            weightflag = 0
        else:
            weightflag = 1
        if DEBUG:
            print("STILL TO HANDLE DERIVATIVES")
        model_deriv = self.modelFunctionDerivative
        if self._fitConfiguration['fit']['strip_flag']:
            y = self._y - self._z
        else:
            y = self._y
        self._fitResult = None
        if not flagconstrained:
            param_constrains = []
        if DEBUG:
            result = Gefit.LeastSquaresFit(self.modelFunction,param_val,
                    xdata=self._x,
                    ydata=y,
                    sigmadata=self._sigma,
                    constrains=param_constrains,
                    weightflag=weightflag,
                    model_deriv=model_deriv,
                    fulloutput=True)
        else:
            try:
                result = Gefit.LeastSquaresFit(self.modelFunction,param_val,
                        xdata=self._x,
                        ydata=y,
                        sigmadata=self._sigma,
                        constrains=param_constrains,
                        weightflag=weightflag,
                        model_deriv=model_deriv,
                        fulloutput=True)
            except:
                text = sys.exc_info()[1]
                if type(text) is not type(" "): 
                    text = text.args
                    if len(text):
                        text = text[0]
                    else:
                        text = ''
                self._setStatus('Fit error : %s' %text)
                raise

        self._fitResult = {}
        self._fitResult['fit_function'] = self.getFitFunction()
        self._fitResult['background_function'] = self.getBackgroundFunction()
        self._fitResult['fittedvalues'] = result[0]
        self._fitResult['chisq']        = result[1]
        self._fitResult['sigma_values'] = result[2]
        self._fitResult['niter']        = result[3]
        self._fitResult['lastdeltachi'] = result[4]
        self._fitResult['n_background_parameters'] = self.__nBackgroundParameters
        if DEBUG:
            print("Found parameters = ", self._fitResult['fittedvalues'])
        i=0
        self._fitResult['parameters'] = []
        for param in self.paramlist:
           if param['code'] != 'IGNORE':
              self._fitResult['parameters'].append(param['name'])  
              param['fitresult'] = result[0][i]
              param['sigma']= result[2][i]
           i = i + 1
        self._setStatus("Fit finished")           
        return result

    def modelFunction(self, pars, t):
        result = 0.0 * t

        nb = self.__nBackgroundParameters
        if nb:
            result += self._fitConfiguration['functions'][self.getBackgroundFunction()]\
                      ['function'](pars[0:nb], t)
        if len(self.paramlist) > nb:
            result += self._fitConfiguration['functions'][self.getFitFunction()]\
                      ['function'](pars[nb:], t)
        return result

    def numericDerivative(self, f, parameters, index, x):
        """
        numericDerivative(self, f, parameters, index, x)
        calculates the numeric derivate of f(parameters, x) respect
        to the parameter indexed by the given index
        """
        #numeric derivative
        x=numpy.array(x)
        delta = (parameters[index] + numpy.equal(parameters[index],0.0)) * 0.00001

        #make a copy of the parameters
        newpar = parameters * 1
        newpar[index] = parameters[index] + delta
        f1 = f(newpar, x)
        newpar[index] = parameters[index] - delta
        f2 = f(newpar, x)
        return (f1-f2) / (2.0 * delta)

    def modelFunctionDerivative(self, pars, index, x):
        return self.numericDerivative(self.modelFunction, pars, index, x)

    def getResult(self, configuration=False):
        #print " get results to be implemented"
        ddict = {}
        ddict['result'] = self._fitResult
        if configuration:
            ddict['configuration'] = self.getConfiguration()
        return ddict        

    def _evaluateBackground(self, x=None):
        if x is None:
            x = self._x
        pars = self._fitResult['fittedvalues']
        nb = self.__nBackgroundParameters
        if nb:
            y = self._fitConfiguration['functions'][self.getBackgroundFunction()]\
                      ['function'](pars[:nb], x)
            
        else:
            y = numpy.zeros(x.shape, numpy.float)
        if self._fitConfiguration['fit']['strip_flag']:
            #If the x is not self._x, how to add the strip?
            try:
                y += self._z
            except:
                print("Cannot add strip background")
        return y

    def _evaluateFunction(self, x=None):
        if x is None:
            x = self._x
        pars = self._fitResult['fittedvalues']
        nb = self.__nBackgroundParameters
        if len(self.paramlist) > nb:
            return self._fitConfiguration['functions'][self.getFitFunction()]\
                      ['function'](pars[nb:], x)
        else:
            return numpy.zeros(x.shape, numpy.float)

    def evaluateDefinedFunction(self, x=None):
        if x is None:
            x = self._x
        y = self._evaluateBackground(x)
        y += self._evaluateFunction(x)
        return y 

def test():
    from PyMca import SpecfitFunctions
    a=SpecfitFunctions.SpecfitFunctions()
    x = numpy.arange(1000).astype(numpy.float)
    p1 = numpy.array([1500,100.,50.0])
    p2 = numpy.array([1500,700.,50.0])
    y = a.gauss(p1, x)+1
    y = y + a.gauss(p2,x)

    fit = SimpleFit()
    fit.importFunctions(SpecfitFunctions)
    fit.setFitFunction('Gaussians')
    #fit.setBackgroundFunction('Gaussians')
    #fit.setBackgroundFunction('Constant')
    fit.setData(x, y)
    fit.fit()
    print("Expected parameters 1500,100.,50.0, 1500,700.,50.0")
    from PyMca import PyMcaQt as qt
    from PyMca import Parameters
    a = qt.QApplication(sys.argv)
    qt.QObject.connect(a,qt.SIGNAL("lastWindowClosed()"),a,qt.SLOT("quit()"))
    w =Parameters.Parameters()
    w.fillfromfit(fit.paramlist)
    w.show()
    a.exec_()

if __name__=="__main__":
    DEBUG = 1
    test()
