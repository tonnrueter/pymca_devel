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

DEBUG = 0
inputDir  = None
outputDir = None
nativeFileDialogs = False

class __ModuleWrapper:
  def __init__(self, wrapped):
    self.__dict__["_ModuleWrapper__wrapped"] = wrapped

  def __getattr__(self, name):
    if DEBUG:
        print("getting ", name)
    if name == "inputDir":
        if self.__wrapped.__dict__[name] is None:
            if self.__wrapped.__dict__['outputDir'] is not None:
                value = self.__wrapped.__dict__['outputDir']
            else:
                value = os.getcwd()
            if not os.path.isdir(value):
                value = os.getcwd()
            self.__setattr__('inputDir', value)
    elif name == "outputDir":
        if self.__wrapped.__dict__[name] is None:
            if self.__wrapped.__dict__['inputDir'] is not None:
                value = self.__wrapped.__dict__['inputDir']
            else:
                value = os.getcwd()
            if not os.path.isdir(value):
                value = os.getcwd()
            self.__setattr__('outputDir', value)
    if DEBUG:
        print("got ", name, getattr(self.__wrapped, name))
    return getattr(self.__wrapped, name)

  def __setattr__(self, name, value):
    if DEBUG:
        print("setting ", name, value)
    if name == "inputDir":
        if os.path.isdir(value):
            self.__wrapped.__dict__[name]=value
        else:
            if not len("%s" % value):
                self.__wrapped.__dict__[name] = os.getcwd()
            else:  
                raise ValueError("Non existing directory %s" % value)
    elif name == "outputDir":
        if os.path.isdir(value):
            self.__wrapped.__dict__[name]=value
        else:
            if not len("%s" % value):
                self.__wrapped.__dict__[name] = os.getcwd()
            else:  
                raise ValueError("Non existing directory %s" % value)
    elif name == "nativeFileDialogs":
        self.__wrapped.__dict__[name]=value
    elif name.startswith("__"):
        self.__dict__[name]=value
    else:
        raise AttributeError("Invalid attribute %s" % name)
        #self.__wrapped.__dict__[name]=value

sys.modules[__name__]=__ModuleWrapper(sys.modules[__name__])


