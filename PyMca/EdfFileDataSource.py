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
from PyMca import DataObject
from PyMca import EdfFile
import types
import sys
import os
import numpy

SOURCE_TYPE = "EdfFile"
DEBUG = 0

class EdfFileDataSource(object):
    def __init__(self,nameInput, fastedf=False):
        if type(nameInput) == list:
            nameList = nameInput
        else:
            nameList = [nameInput]
        if sys.version < '3.0':
            stringTypes = [types.StringType, types.UnicodeType]
        else:
            stringTypes = [type("a"), type(eval('b"a"'))]
        for name in nameList:
            if type(name) not in stringTypes:
                raise TypeError("Constructor needs string as first argument")           
        self.sourceName   = nameInput
        self.sourceType = SOURCE_TYPE
        self.__sourceNameList = nameList
        #this is to be added
        #self._fastedf = True
        self._fastedf = fastedf
        if fastedf:
            print("fastedf is unsafe!")
        self.refresh()

    def refresh(self):
        self._sourceObjectList=[]
        for name in self.__sourceNameList:
            self._sourceObjectList.append(EdfFile.EdfFile(name, access='rb', fastedf=self._fastedf))
        self.__lastKeyInfo = {}

    def getSourceInfo(self):
        """
        Returns information about the EdfFile object created by
        SetSource, to give application possibility to know about
        it before loading.
        Returns a dictionary with the key "KeyList" (list of all available keys
        in this source). Each element in "KeyList" has the form 'n1.n2' where
        n1 is the source number and n2 image number in file starting at 1.
        """        
        return self.__getSourceInfo()
        
        
    def __getSourceInfo(self):
        SourceInfo={}
        SourceInfo["SourceType"]=SOURCE_TYPE
        SourceInfo["KeyList"]=[]
        i = 0
        for sourceObject in self._sourceObjectList:
            i+=1
            nimages = sourceObject.GetNumImages()
            for n in range(nimages):
                SourceInfo["KeyList"].append("%d.%d" % (i,n+1))   
        SourceInfo["Size"]=len(SourceInfo["KeyList"])
        return SourceInfo
        
    def getKeyInfo(self,key):
        if key in self.getSourceInfo()['KeyList']:
            return self.__getKeyInfo(key)
        else:
            #should we raise a KeyError?
            if DEBUG:
                print("Error key not in list ")
            return {}
    
    def __getKeyInfo(self,key):
        try:
            index,image = key.split(".")
            index = int(index)-1
            image = int(image)-1            
        except:
            #should we rise an error?
            if DEBUG:
                print("Error trying to interpret key =",key)
            return {}

        sourceObject = self._sourceObjectList[index]
        info= sourceObject.GetStaticHeader(image)
        info.update(sourceObject.GetHeader(image))
        info["SourceType"]  = SOURCE_TYPE
        #doubts about if refer to the list or to the individual file
        info["SourceName"]  = self.sourceName
        info["Key"]         = key
        #specific info of interest
        info['FileName'] = sourceObject.FileName
        info["rows"]     = info['Dim_2']
        info["cols"]     = info['Dim_1']
        info["type"]     = info['DataType']
        if 'MCA start ch' in info:
            info['Channel0'] = int(info['MCA start ch'])
        else:
            info['Channel0'] = 0

        if not ( 'McaCalib' in info):
            if ('MCA a' in info) and\
               ('MCA b' in info) and\
               ('MCA c' in info):
                info['McaCalib'] = [float(info['MCA a']),
                                    float(info['MCA b']),
                                    float(info['MCA c'])]
            else:
                info['McaCalib'] = [ 0.0, 1.0, 0.0]
        else:
            if type(info['McaCalib']) in [type(" ")]:
                info['McaCalib'] = info['McaCalib'].replace("[","")
                info['McaCalib'] = info['McaCalib'].replace("]","")
                cala, calb, calc = info['McaCalib'].split(",")
                info['McaCalib'] = [float(cala),
                                    float(calb),
                                    float(calc)]

        self.__lastKeyInfo[key] = os.path.getmtime(sourceObject.FileName)
        return info

    def getDataObject(self,key,selection=None):
        """
        selection: a dictionnary with the keys pos and size: (x), (x,y) or (x,y,z)
                   tuples defining a roi
                   If not defined, takes full array
        """

        sourcekeys = self.getSourceInfo()['KeyList']
        #a key corresponds to an image        
        key_split= key.split(".")
        image_key= key_split[0]+"."+key_split[1]
        if image_key not in sourcekeys:
            #if image_key == "0.0":
            #    #this is in fact a special selection: The SUM
            #    pass
            #else:
                raise KeyError("Key %s not in source keys" % image_key)
        #create data object
        data = DataObject.DataObject()
        data.info = self.__getKeyInfo(image_key)
        data.info ['selection'] = selection
        data.info['selectiontype'] = "2D"
        index = key_split[0]
        image = key_split[1]
        index = int(index)-1
        image = int(image)-1
        MCAIMP = 0
        if len(key_split) == 4:
            if DEBUG:
                print("mca like selection")
            #print data.info
            if 1:
                MCAIMP = 1
                if key_split[2].upper() == 'R':
                    pos  = (0, int(key_split[3]))
                    size = (int(data.info['Dim_1']), 1) 
                elif key_split[2].upper() == 'C':
                    pos  = (int(key_split[3]), 0)
                    size = (1,int(data.info['Dim_2']))
                data.info['selectiontype'] = "1D"
            else:
                if DEBUG:
                    print("mca like selection not yet implemented")
                pos = None
                size = None
                data.info['selectiontype'] = "1D"
           
        elif selection is None:
            pos = None
            size = None
        else:
            if "pos" in selection:
                data.info["pos"]=selection['pos']
            else:
                data.info['pos']=None
            if "size" in selection:
                data.info["size"]=selection['size']
            else:
                data.info['size']=None
            pos  = data.info['pos']
            size = data.info['size']
        sourceObject = self._sourceObjectList[index]
        data.data=sourceObject.GetData(image,Pos=pos,Size=size)
        data.info['rows'], data.info['cols'] = data.data.shape[0:2]
        if data.info['selectiontype'] == "1D":
            if MCAIMP:
                data.y = [numpy.ravel(data.data[:]).astype(numpy.float)]
            else:
                if key_split[2].upper() == 'C':
                    data.y=[data.data[:,int(key_split[3])-1].astype(numpy.float)]
                elif key_split[2].upper() == 'R':
                    data.y=[data.data[int(key_split[3])-1, :].astype(numpy.float)]
                else:
                    raise ValueError("Unknown key %s" % key)
            ch0 = int(data.info['Channel0'])
            data.x = [ch0+numpy.arange(len(data.y[0])).astype(numpy.float)]
            data.m = None
            data.data = None
            #print "data.x.shape ", data.x[0].shape
            #print "data.y.shape ", data.y[0].shape
        return data

    def isUpdated(self, sourceName, key):
        #sourceName is redundant?
        index,image = key.split(".")
        index = int(index)-1
        lastmodified = os.path.getmtime(self.__sourceNameList[index])
        if lastmodified != self.__lastKeyInfo[key]:
            self.__lastKeyInfo[key] = lastmodified
            return True
        else:
            return False

source_types = { SOURCE_TYPE: EdfFileDataSource}

def DataSource(name="", source_type=SOURCE_TYPE):
  try:
     sourceClass = source_types[source_type]
  except KeyError:
     #ERROR invalid source type
     raise TypeError("Invalid Source Type, source type should be one of %s" %\
                     source_types.keys())
  
  return sourceClass(name)

        
if __name__ == "__main__":
    import time
    try:
        sourcename=sys.argv[1]
        key       =sys.argv[2]        
    except:
        print("Usage: EdfFileDataSource <file> <key>")
        sys.exit()
    #one can use this:
    obj = EdfFileDataSource(sourcename)
    #or this:
    obj = DataSource(sourcename)
    #data = obj.getData(key,selection={'pos':(10,10),'size':(40,40)})
    #data = obj.getDataObject(key,selection={'pos':None,'size':None})
    t0 = time.time()
    data = obj.getDataObject(key,selection=None)
    print("elapsed = ",time.time() - t0)
    print("info = ",data.info)
    if data.data is not None:
        print("data shape = ",data.data.shape)
        print(numpy.ravel(data.data)[0:10])
    else:
        print(data.y[0].shape)
        print(numpy.ravel(data.y[0])[0:10])
    data = obj.getDataObject('1.1',selection=None)
    r = int(key.split('.')[-1])
    print(" data[%d,0:10] = " % (r-1),data.data[r-1   ,0:10])
    print(" data[0:10,%d] = " % (r-1),data.data[0:10, r-1])
