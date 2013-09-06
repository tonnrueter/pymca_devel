#/*##########################################################################
# Copyright (C) 2004-2012 European Synchrotron Radiation Facility
#
# This file is part of the PyMca X-ray Fluorescence Toolkit developed at
# the ESRF by the Software group.
#
# This file is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# This file is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
#############################################################################*/
__author__ = "V.A. Sole - ESRF Data Analysis"
import os
import numpy
from PyMca import ConfigDict
from PyMca import PyMcaDataDir

ElementList= ['H','He','Li','Be','B','C','N','O','F','Ne',
              'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     'Ga','Ge','As','Se','Br','Kr',
     'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
     'In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd',
     'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
     'Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At',
     'Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf',
     'Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt']

dirmod = PyMcaDataDir.PYMCA_DATA_DIR
ffile   = os.path.join(dirmod,"attdata")
ffile   = os.path.join(ffile,"incoh.dict")
if not os.path.exists(ffile):
    #freeze does bad things with the path ...
    dirmod = os.path.dirname(dirmod)
    ffile = os.path.join(dirmod, "attdata")
    ffile = os.path.join(ffile, "incoh.dict")
    if not os.path.exists(ffile):
        if dirmod.lower().endswith(".zip"):
            dirmod = os.path.dirname(dirmod)
            ffile = os.path.join(dirmod,"attdata")
            ffile = os.path.join(ffile, "incoh.dict")
    if not os.path.exists(ffile):
        print("Cannot find file ", ffile)
        raise IOError("Cannot find file %s" % ffile)

COEFFICIENTS = ConfigDict.ConfigDict()
COEFFICIENTS.read(ffile)
xvalues = COEFFICIENTS['ISCADT']['XSVAL']
svalues = numpy.reshape(COEFFICIENTS['ISCADT']['SCATF'], (100, len(xvalues)))
#svalues = COEFFICIENTS['ISCADT']['SCATF']
#print svalues[100:110]
KEVTOANG = 12.39852000
R0 = 2.82E-13 #electron radius in cm

def getZ(ele):
    if ele in ElementList:
        return float(ElementList.index(ele)+1)
    else:
        return None

def getElementComptonFormFactor(ele, theta, energy):
    return getElementIncoherentScatteringFunction(ele, theta, energy)


def getComptonScatteringEnergy(energy, theta):
    return energy/(1.0 + \
            (energy/511.) * (1 - numpy.cos(theta*(numpy.pi / 180.0))))

def getElementIncoherentScatteringFunction(ele, theta, energy):
    """
    Usage: 
        getIncoherentScatteringFunction(ele,theta, energy):
    
    ele   - Element
    theta - Scattering angle in degrees
    energy- Photon Energy in keV
    
    This routine calculates the incoherent scattering function 
    in electron units an interpolation to EGS4 tabulation of S(x,Z)/Z
    """
    if ele in ElementList:
        z = getZ(ele)
    else:
        z = float(ele)
    wavelength = KEVTOANG / energy
    sinhalftheta = numpy.sin(theta * (numpy.pi / 360.0))
    #Hubbel just give this term
    x =  sinhalftheta / wavelength
    
    #print "x old = ",x
    e = energy/511.0
    #Fajardo uses:
    x = x * numpy.sqrt(1.0 + e* (e+2.0)* pow(sinhalftheta, 2))/ \
            (1.0 + 2.0 * e * pow(sinhalftheta, 2))
    #print "x new = ",x
    
    ilow  = 0
    ihigh = 44
    i     = 22
    while (ihigh - ilow) > 1:
        if x < xvalues[i]:ihigh = i
        else:ilow =i
        i = int((ihigh+ilow)/2)

    if z > 100:
        if ihigh == ilow:
            value = svalues[int(99),ilow]  
        else:
            A = (x - xvalues[ilow])/(xvalues[ihigh]-xvalues[ilow])
            value = ((1.0 - A ) * svalues[int(99),ilow] + \
                    A * svalues[int(99),ihigh])
        value = value * (z/100.)
    else:
        if ihigh == ilow:
            value = svalues[int(z-1),ilow]  
        else:
            A = (x - xvalues[ilow])/(xvalues[ihigh]-xvalues[ilow])
            value = ((1.0 - A ) * svalues[int(z-1),ilow] + \
                    A * svalues[int(z-1),ihigh])
    
    return value



def getElementComptonDifferentialCrossSection(ele, theta, energy, p1=None):
    if p1 is None:p1=0.0
    if (p1 > 1.0) or (p1 < -1):
        raise ValueError(\
        "Invalid degree of linear polarization respect to the scattering plane")
    thetasin2 = pow(numpy.sin(theta * numpy.pi / 180.0), 2)
    thetacos  =  numpy.cos(theta * numpy.pi/180.0)
    e = energy/(1.0 + (energy/511.) * (1.0 - thetacos))
    return 0.5 * ((e/energy) + (energy/e) + (p1-1.0) * thetasin2) * \
           pow(R0*(e/energy)*getElementIncoherentScatteringFunction(ele, theta, energy),2)

getElementIncoherentDifferentialCrossSection=\
            getElementComptonDifferentialCrossSection

if __name__ == "__main__":
    import sys
    if len(sys.argv) >  3:
        ele   = sys.argv[1]
        theta = float(sys.argv[2])
        energy= float(sys.argv[3])
        print(getElementComptonFormFactor(ele, theta, energy))    
    else:
        print("Usage:")
        print("python IncoherentScatteringFunction.py Element Theta(deg) Energy(kev)")
