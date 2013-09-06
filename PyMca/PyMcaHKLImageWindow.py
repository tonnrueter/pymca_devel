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
__author__ = "V.A. Sole - ESRF Data Analysis"
import sys
import numpy
from PyMca import PyMcaImageWindow
from PyMca import SixCircle

arctan = numpy.arctan

DEBUG = 0

class PyMcaHKLImageWindow(PyMcaImageWindow.PyMcaImageWindow):
    def __init__(self, *var, **kw):
        PyMcaImageWindow.PyMcaImageWindow.__init__(self, *var, **kw)
        self._HKLOn = False

    def _graphSignal(self, ddict):
        if (ddict['event'] != "MouseAt") or (not self._HKLOn):
            return PyMcaImageWindow.PyMcaImageWindow._graphSignal(self, ddict)

        if self._imageData is None:
            self.graphWidget.setInfoText("    H = ???? K = ???? L = ???? I = ????")
            return

        #pixel coordinates
        x = round(ddict['y'])
        if x < 0: x = 0
        y = round(ddict['x'])
        if y < 0: y = 0
        limits = self._imageData.shape
        x = min(int(x), limits[0]-1)
        y = min(int(y), limits[1]-1)
        z = self._imageData[x, y]

        text = "   X = %d Y = %d Z = %.7g " % (y, x, z)

        info = self._getHKLInfoFromWidget()

        toDeg = 180.0/numpy.pi

        phi = info['phi']
        chi = info['chi']
        theta = info['theta']

        if 0:
            # delta in vertical (following BM28)
            # gamma in horizontal (following BM28)
            deltaH = toDeg * numpy.arctan((x - info['pixel_zero_h']) *\
                        (info['pixel_size_h']/info['distance']))
        
            deltaV = toDeg *arctan((y - info['pixel_zero_v'])*\
                        (info['pixel_size_v']/info['distance']))
            if 0:
               #original
               gamma = info['gamma'] + deltaH
               delta = info['delta'] - deltaV
            else:
               #MarCCD settings
               gamma = info['gamma'] - deltaV
               delta = info['delta'] - deltaH            
           #end of BM28 customization
        else:
	    #ID03
            deltaH = toDeg * numpy.arctan((x - info['pixel_zero_v']) *\
                        (info['pixel_size_v']/info['distance']))
        
            deltaV = toDeg *arctan((y - info['pixel_zero_h'])*\
                        (info['pixel_size_h']/info['distance']))
	    #delta in horizontal
	    #gamma in vertical
            gamma = info['gamma'] - deltaH
            delta = info['delta'] - deltaV
            if 0:
                #ID03 test for EH1
                wavelength = 1.03321027
                ub = [1.0, 0.0, 0.0, 
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0]
                ub[0] = 0.060082400000000001
                ub[1] = 0.054556500000000001
                ub[2] = -0.92985700000000004
                ub[3] =  -1.5089399999999999
                ub[4] =  -2.61991
                ub[5] =  -0.0203886
                ub[6] =  -2.1539600000000001
                ub[7] =  0.230518
                ub[8] = -0.011654299999999999
                delta, theta, chi, phi, mu, gamma = 44.0035, -92.968, 90.715,\
                                                    1.26, 0.3, 0.578
                print(" Expected value = ", 1, 1, 0.1)
        mu    = info['mu']
        wavelength = info['lambda']
        ub = info['ub']

        if 0:
            #This should always give 1 1 1 
            wavelength = 0.363504
            ub = [1.0, 0.0, 0.0, 
                  0.0, 1.0, 0.0,
                  0.0, 0.0, 1.0]
            ub[0] = -4.080
            ub[1] =  0.000
            ub[2] =  0.000
            ub[3] =  0.000
            ub[4] =  4.080
            ub[5] =  0.000
            ub[6] =  0.000
            ub[7] =  0.000
            ub[8] = -4.080
            delta, theta, chi, phi, mu, gamma = 23.5910, 47.0595, -135.,\
                                                0.0, 0.0, 0.0

        HKL = SixCircle.getHKL(wavelength,
                               ub,
                               phi=phi,
                               chi=chi,
                               theta=theta,
                               gamma=gamma,
                               delta=delta,
                               mu=mu)
        HKL.shape = -1
        text += "H = %.3f " % HKL[0]
        text += "K = %.3f " % HKL[1]
        text += "L = %.3f " % HKL[2]
        self.graphWidget.setInfoText(text)


    def _getHKLInfoFromWidget(self):
        ddict = {}
        ddict['lambda'] = 1.0           # In Angstroms   
        ddict['distance'] = 1000.       # Same units as pixel size
        ddict['pixel_size_h'] = 0.080   # Same units as distance
        ddict['pixel_size_v'] = 0.080   # Same units as distance
        ddict['pixel_zero_h'] = 1024.   # In pixel units (float)
        ddict['pixel_zero_v'] = 1024.   # In pixel units (float)
        ddict['orientation'] = 0
        ddict['ub'] = [1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0]
        ddict['phi'] = 0.0
        ddict['chi'] = 0.0 
        ddict['theta'] = 0.0
        ddict['gamma'] = 0.0
        ddict['delta'] = 0.0
        ddict['mu']    = 0.0

        legend = self.dataObjectsList[0]
        dataObject = self.dataObjectsDict[legend]
        info = dataObject.info

        #try to get the information from the motors
        motPos = info.get('motor_pos', "")
        motMne = info.get('motor_mne', "")
        motPos = motPos.split()
        motMne = motMne.split()
        if len(motPos) == len(motMne):
            idx = -1
            for mne in motMne:
                idx += 1
                if mne.upper() in ['ENERGY', 'NRJ']:
                    energy = float(motPos[idx])
                    ddict['lambda'] = 12.39842 / energy
                    continue
                if mne in ['phi', 'chi', 'mu']:
                    ddict[mne] = float(motPos[idx])
                    continue
                if mne in ['th', 'theta']:
                    ddict['theta'] = float(motPos[idx])
                    continue
                if mne in ['del', 'delta', 'tth', 'twotheta']:
                    ddict['delta'] = float(motPos[idx])
                    continue
                if mne in ['gam', 'gamma']:
                    ddict['gamma'] = float(motPos[idx])
                    continue

        #and update it from the counters
        cntPos = info.get('counter_pos', "").split()
        cntMne = info.get('counter_mne', "").split()
        cntInfo = {}
        if len(cntPos) == len(cntMne):
            for i in range(len(cntMne)):
                cntInfo[cntMne[i]] = cntPos[i]
            
        for key in cntInfo.keys():
            # diffractometer
            if key in ['phicnt']:
                ddict['phi'] = float(cntInfo[key])
                continue
            if key in ['chicnt']:
                ddict['chi'] = float(cntInfo[key])
                continue
            if key in ['thcnt', 'thetacnt']:
                ddict['theta'] = float(cntInfo[key])
                continue
            if key in ['tthcnt'] and ('delcnt' not in cntInfo.keys()):
                #Avoid ID03 trap because they have delcnt and tthcnt ...
                ddict['delta'] = float(cntInfo[key])
                continue
            if key in ['delcnt', 'deltacnt']:
                ddict['delta'] = float(cntInfo[key])
                continue
            if key in ['gamcnt', 'gammacnt']:
                ddict['gamma'] = float(cntInfo[key])
                continue
            if key in ['mucnt']:
                ddict['mu'] = float(cntInfo[key])
                continue

        for key in info.keys():
            # UB matrix
            if key.upper() in ['UB_POS']:
                ddict['ub'] = [float(x) for x in info[key].split()]
                continue

            # direct beam
            if key in ['beam_x', 'pixel_zero_x']:
                ddict['pixel_zero_h'] = float(info[key])
                continue
            if key in ['beam_y', 'pixel_zero_y']:
                ddict['pixel_zero_v'] = float(info[key])
                continue

            #sample to direct beam distance
            if key in ['detector_distance', 'd_sample_det']:
                ddict['distance'] = float(info[key])
                continue

            #pixel sizes
            if key in ['pixel_size_x']:
                ddict['pixel_size_h'] = float(info[key])
                continue

            if key in ['pixel_size_y']:
                ddict['pixel_size_v'] = float(info[key])
                continue

            #wave length
            if key in ['source_wavelength']:
                ddict['lambda'] = float(info[key])
                continue

        if DEBUG:
            for key in ddict.keys():
                print(key, ddict[key])

        return ddict
