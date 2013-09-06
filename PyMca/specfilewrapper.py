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
import os
import sys
import numpy
import re
try:
    from PyMca import specfile
except ImportError:
    print("specfilewrapper importting specfile directly")
    import specfile
from PyMca import Fit2DChiFileParser
from PyMca import APSMEDFileParser
from PyMca import SRSFileParser
try:
    from PyMca import SPXFileParser
    SPX = True
except:
    SPX = False

DEBUG = 0

if sys.version >= '2.6':
    def safe_str(bytesObject):
        try:
            return str(bytesObject, 'utf-8')
        except UnicodeDecodeError:
            try:
                return str(bytesObject, 'latin-1')
            except:
                try:
                    return str(bytesObject, 'utf-16')
                except:
                    return str(bytesObject)
else:
    def safe_str(*var, **kw):
        return str(var[0])

    #python 2.5 does not have bytes function
    def bytes(*var, **kw):
        return var[0]

def Specfile(filename):
    f = open(filename)
    line0  = f.readline()
    if filename.upper().endswith('DTA'):
        #TwinMic single column file
        line = line0 * 1
        line = line.replace('\r','')
        line = line.replace('\n','')
        line = line.replace('\t',' ')
        s = line.split(' ')
        if len(s) == 2:
            if len(s[-1]) == 0:
                try:
                    float(s[0])
                    f.close()
                    output = specfilewrapper(filename, dta=True)
                    f.close()
                    return output
                except:
                    #try to read in other way
                    pass
    line = line0
    while(len(line)):
        if len(line) > 1:
            if line[0:2] == '#S':
                break
        line = f.readline()
    f.close()
    amptek = False
    qxas   = False
    if len(line):
        #it is a Specfile
        output=specfile.Specfile(filename)
    elif SPX and filename.upper().endswith("SPX"):
        #spx file
        output = SPXFileParser.SPXFileParser(filename)
    else:
        if DEBUG:
            print("this does not look as a specfile")
        if len(line0) > 7:
            if line0.startswith('$SPEC_ID') or\
               line0.startswith('$DATE_MEA') or\
               line0.startswith('$MEAS_TIM'):
                qxas = True
        if (not qxas) and line0.startswith('<<'):
                amptek = True
        if (not qxas) and (not amptek) and Fit2DChiFileParser.isFit2DChiFile(filename):
            return Fit2DChiFileParser.Fit2DChiFileParser(filename)
        if (not qxas) and (not amptek) and APSMEDFileParser.isAPSMEDFile(filename):
            return APSMEDFileParser.APSMEDFileParser(filename)
        if (not qxas) and (not amptek) and SRSFileParser.isSRSFile(filename):
            return SRSFileParser.SRSFileParser(filename)
        output = specfilewrapper(filename, amptek=amptek, qxas=qxas)
    return output

class specfilewrapper(object):
    def __init__(self, filename, amptek=None, qxas=None, dta=None):
        if amptek is None:
            amptek = False
        if qxas is None:
            qxas = False
        if dta is None:
            dta = False
        self.amptek = amptek
        self.qxas   = qxas
        self.dta = dta
        self.header = []
        if self.dta:
            #TwinMic .dta files with only one spectrum
            if 0:
                f = open(filename, 'rb')
                raw_content = f.read()
                f.close()
                expr = '([-+]?\d+)\t\r\n'
                self.data = [float(i) for i in re.split(expr,raw_content) if i != '']
                self.data = numpy.array(self.data, numpy.float32)
            else:
                self.data = numpy.fromfile(filename,
                                           dtype=numpy.float32,
                                           sep='\t\r\n')
            self.header = ['#S1 %s' % os.path.basename(filename)]
            self.data.shape = -1, 1
            self.scandata=[myscandata(self.data,'MCA','1.1',
                                      scanheader=self.header)]
            return
        if self.qxas:
            f = open(filename)
        else:
            f = BufferedFile(filename)
        line = f.readline()
        outdata = []
        ncol0 = -1
        nlines= 0
        if amptek:
            if sys.version < '3.0':
                while "<<DATA>>" not in line:
                    self.header.append(line.replace("\n",""))
                    line = f.readline()
            else:
                while bytes("<<DATA>>", 'utf-8') not in line:
                    self.header.append(safe_str(line.replace(bytes("\n", 'utf-8'),
                                                    bytes("", 'utf-8'))))
                    line = f.readline()
        elif qxas:
            line.replace("\n","")
            line.replace("\x00","")
            self._qxasHeader = {}
            self._qxasHeader['S'] = '#S1 '+ " Unlabelled Spectrum"
            while 1:
                self.header.append(line)
                if line.startswith('$SPEC_ID:'):
                    line = f.readline().replace("\n","")
                    line.replace("\x00","")
                    self.header.append(line)
                    self._qxasHeader['S'] = '#S1 '+ line
                if line.startswith('$DATE_MEA'):
                    line = f.readline().replace("\n","")
                    self.header.append(line)
                    self._qxasHeader['D'] = line
                if line.startswith('$MEAS_TIM'):
                    line = f.readline().replace("\n","")
                    self.header.append(line)
                    tmpList = [float(i) for i in line.split()]
                    if len(tmpList) == 1:
                        preset = tmpList[0]
                        elapsed = preset
                    else:
                        preset, elapsed = tmpList[0:2]
                    self._qxasHeader['@CTIME'] = ['#@CTIME %f %f %f' % (preset, preset, elapsed)]
                if line.startswith('$MCA_CAL'):
                    try:
                       line = f.readline().replace("\n","")
                       self.header.append(line)
                       if line.startswith('$'):
                           continue
                       line = f.readline().replace("\n","")
                       self.header.append(line)
                       if line.startswith('$'):
                           continue
                       coefficients = [float(i) for i in line.split()]
                       if len(coefficients) == 2:
                           coefficients.append(0.0)
                       self._qxasHeader['@CALIB']=  ['#@CALIB %f  %f  %f' %\
                                            (coefficients[0], coefficients[1], coefficients[2])]
                    except:
                        pass
                if line.startswith('$DATA:'):
                    line = f.readline().replace("\n","")
                    self.header.append(line)
                    start, stop = [int(i) for i in line.split()]
                    self._qxasHeader['@CHANN'] = ['#@CHANN  %d  %d  %d  1' % (stop-start+1, start, stop)]
                    break
                line = f.readline().replace("\n","")
        if qxas:
            outdata = [float(x) for x in f.read().split()]
            nlines = len(outdata)
            f.close()
            self.data = numpy.resize(numpy.array(outdata).astype(numpy.float),(nlines,1))
        else:
            if sys.version < '3.0':
                line = line.replace(",","  ")
                line = line.replace(";","  ")
                line = line.replace("\t","  ")
                line = line.replace("\r","\n")
                line = line.replace('"',"")
                line = line.replace('\n\n',"\n")
            else:
                tmpBytes = bytes(" ",'utf-8')
                line = line.replace(bytes(",",'utf-8'), tmpBytes)
                line = line.replace(bytes(";",'utf-8'), tmpBytes)
                line = line.replace(bytes("\t",'utf-8'), tmpBytes)
                tmpBytes = bytes("\n",'utf-8')
                line = line.replace(bytes("\r","utf-8"), tmpBytes)
                line = line.replace(bytes('"',"utf-8"), bytes("", "utf-8"))
                line = line.replace(bytes('\n\n',"utf-8"), tmpBytes)
            while(len(line)):
                values = line.split()
                if len(values):
                    try:
                        reals = [float(x) for x in values]
                        ncols = len(reals)
                        if ncol0 < 0:ncol0 = ncols
                        if ncols == ncol0:
                            outdata.append(reals)
                            nlines += 1
                    except:
                        if len(line) > 1:
                            if sys.version < '3.0':
                                self.header.append(line.replace("\n",""))
                            else:
                                self.header.append(safe_str(line.replace(\
                                                    bytes("\n",'utf-8'),\
                                                    bytes("", 'utf-8'))))
                else:
                    if len(line) > 1:
                        if sys.version < '3.0':
                            self.header.append(line.replace("\n",""))
                        else:
                            self.header.append(safe_str(line.replace(bytes("\n",'utf-8'),
                                                            bytes("", 'utf-8'))))
                line = f.readline()
                if sys.version < '3.0':
                    line = line.replace(",","  ")
                    line = line.replace(";","  ")
                    line = line.replace("\t","  ")
                    line = line.replace("\r","\n")
                    line = line.replace('"',"")
                    line = line.replace('\n\n',"\n")
                else:
                    tmpBytes = bytes(" ",'utf-8')
                    line = line.replace(bytes(",",'utf-8'), tmpBytes)
                    line = line.replace(bytes(";",'utf-8'), tmpBytes)
                    line = line.replace(bytes("\t",'utf-8'), tmpBytes)
                    tmpBytes = bytes("\n",'utf-8')
                    line = line.replace(bytes("\r","utf-8"), tmpBytes)
                    line = line.replace(bytes('"',"utf-8"), bytes("", "utf-8"))
                    line = line.replace(bytes('\n\n',"utf-8"), tmpBytes)
            f.close()
            self.data = numpy.resize(numpy.array(outdata).astype(numpy.float),(nlines,ncol0))
        if self.amptek:
            self.scandata=[myscandata(self.data,'MCA','1.1',
                                      scanheader=self.header)]
        elif self.qxas:
            self.scandata=[myscandata(self.data,'MCA','1.1',
                                      scanheader=self.header,
                                      qxas=self._qxasHeader)]
        else:
            labels = None
            if len(self.header) == 1:
                if len(self.header[0]) > 0:
                    labels = self.header[0].split("  ")
                    if len(labels) != ncol0:
                        labels = None
            self.scandata=[myscandata(self.data,'SCAN','1.1',
                                      labels=labels,
                                      fileheader=self.header),
                           myscandata(self.data,'MCA','2.1',
                                      fileheader=self.header)]

    def list(self):
        if self.amptek or self.qxas or self.dta:
            return "1:1"
        else:
            return "1:2"

    def __getitem__(self,item):
        return self.scandata[item]

    def __len__(self):
        return self.scanno()

    def select(self,i):
        n = i.split(".")
        return self.__getitem__(int(n[0]) - 1)

    def scanno(self):
        if self.amptek or self.qxas:
            return 1
        else:
            return 2

class myscandata(object):
    def __init__(self, data, scantype=None, identification=None,
                 scanheader=None, qxas=None, labels=None, fileheader=None):
        if identification is None:
            identification='1.1'
        if scantype is None:
            scantype='SCAN'
        self.qxas = qxas
        self.scanheader = scanheader
        if fileheader is None:
            fileheader = []
        self._fileheader = fileheader
        #print shape(data)
        (rows, cols) = numpy.shape(data)
        if scantype == 'SCAN':
            self.__data = numpy.zeros((rows, cols +1 ), numpy.float)
            self.__data[:,0] = numpy.arange(rows) * 1.0
            self.__data[:,1:] = data * 1
            self.__cols = cols + 1
            self.labels = ['Point']
            if labels is None:
                for i in range(cols):
                    self.labels.append('Column %d'  % i)
            else:
                for label in labels:
                    self.labels.append('%s' % label)
        else:
            self.__data = data
            self.__cols = cols
            self.labels = []
        self.scantype = scantype
        self.rows = rows
        if scanheader is None:
            labels = '#L '
            for label in self.labels:
                labels += '  '+label
            if self.scantype == 'SCAN':
                self.scanheader = ['#S1 Unknown command',
                              '#N %d' % len(self.labels),
                              labels]
            else:
                self.scanheader = ['#S1 Unknown command']

        n = identification.split(".")
        self.__number = int(n[0])
        self.__order  = int(n[1])

    def alllabels(self):
        if self.scantype == 'SCAN':
            return self.labels
        else:
            return []

    def allmotorpos(self):
        return []

    def cols(self):
        return self.__cols

    def command(self):
        if DEBUG:
            print("command called")
        if self.qxas is not None:
            if 'S' in self.qxas:
                text = self.qxas['S']
        elif self.scanheader is not None:
            if len(self.scanheader):
                text = self.scanheader[0]
        return text

    def data(self):
        return numpy.transpose(self.__data)

    def datacol(self, col):
        # it is awful that starts at one ...
        if col <= 0:
            raise ValueError("Specfile column numberig starts at 1")
        return self.__data[:, col - 1]

    def dataline(self, line):
        # it is awful that starts at one ...
        if line <= 0:
            raise ValueError("Specfile line numberig starts at 1")
        return self.__data[line - 1,:]


    def date(self):
        text = 'sometime'
        if self.qxas is not None:
            if 'D' in self.qxas:
                return self.qxas['D']
        elif self.scanheader is not None:
            for line in self.scanheader:
                if 'START_TIME' in line:
                    text = "%s" % line
                    break
        return text

    def fileheader(self, key=''):
        # key is there for compatibility
        if DEBUG:
            print("file header called")
        return self._fileheader

    def header(self, key):
        if self.qxas is not None:
            if key in self.qxas:
                return self.qxas[key]
            elif key == "" or key == " ":
                return self.scanheader
        if key == 'S':
            return self.scanheader[0]
        elif key == 'N':
            return self.scanheader[-2]
        elif key == 'L':
            return self.scanheader[-1]
        elif key == '@CALIB':
            output = []
            if self.scanheader is None:
                return output
            if self.scanheader[0][0:2] == '<<':
                #amptek
                try:
                    amptekCalibrationLines = []
                    amptekInCalibrationLines = False
                    for line in self.scanheader:
                        if '<<CALIBRATION>>' in line:
                            amptekInCalibrationLines = True
                            continue
                        if line.startswith('<<'):
                            amptekInCalibrationLines = False
                            continue
                        if amptekInCalibrationLines and\
                           ('LABEL' not in line):
                            amptekCalibrationLines.append(line)
                    n = len(amptekCalibrationLines)
                    if n == 0 :
                        return output
                    if n == 1:
                        #one point calibration
                        x0,y0 = 0.0, 0.0
                        values = amptekCalibrationLines[0].split()
                        x1,y1 = map(float,values)
                        gain = (y1-y0)/(x1-x0)
                        zero = y0 - gain * x0
                    elif n == 2:
                        #two point calibration
                        values = amptekCalibrationLines[0].split()
                        x0,y0 = map(float,values)
                        values = amptekCalibrationLines[1].split()
                        x1,y1 = map(float,values)
                        gain = (y1-y0)/(x1-x0)
                        zero = y0 - gain * x0
                    else:
                        x = numpy.zeros((n,), numpy.float)
                        y = numpy.zeros((n,), numpy.float)
                        for i in range(n):
                            values = amptekCalibrationLines[i].split()
                            x[i], y[i] = map(float,values)
                        Sxy = numpy.dot(x, y.T)
                        Sxx = numpy.dot(x, x.T)
                        Sx  = x.sum()
                        Sy  = y.sum()
                        d = n * Sxx - Sx * Sx
                        zero = (Sxx * Sy - Sx * Sxy)/d
                        gain = (n * Sxy - Sx * Sy)/d
                    output = ['#@CALIB  %g  %g  0' % (zero, gain)]
                except:
                    pass
            return output
        elif key == "" or key == " ":
            return self.scanheader
        else:
            return []

    def order(self):
        return self.__order

    def number(self):
        return self.__number

    def lines(self):
        if self.scantype == 'SCAN':
            return self.rows
        else:
            return 0

    def nbmca(self):
        if self.scantype == 'SCAN':
            return 0
        else:
            return self.__cols

    def mca(self,number):
        if number <= 0:
            raise ValueError("Specfile mca numberig starts at 1")
        return self.__data[:,number-1]

class BufferedFile(object):
    def __init__(self, filename):
        f = open(filename, 'rb')
        self.__buffer = f.read()
        f.close()
        if sys.version < '3.0':
            self.__buffer=self.__buffer.replace("\r", "\n")
            self.__buffer=self.__buffer.replace("\n\n", "\n")
            self.__buffer = self.__buffer.split("\n")
        else:
            tmp = bytes("\n", 'utf-8')
            self.__buffer=self.__buffer.replace(bytes("\r", 'utf-8'), tmp)
            self.__buffer=self.__buffer.replace(bytes("\n\n", 'utf-8'), tmp)
            self.__buffer = self.__buffer.split(tmp)
        self.__currentLine = 0

    if sys.version < '3.0':
        def readline(self):
            if self.__currentLine >= len(self.__buffer):
                return ""
            line = self.__buffer[self.__currentLine] + "\n"
            self.__currentLine += 1
            return line
    else:
        def readline(self):
            if self.__currentLine >= len(self.__buffer):
                return bytes("", 'utf-8')
            line = self.__buffer[self.__currentLine] + bytes("\n", 'utf-8')
            self.__currentLine += 1
            return line

    def close(self):
        self.__currentLine = 0
        return

if __name__ == "__main__":
    filename = sys.argv[1]
    print(filename)
    sf=Specfile(filename)
    sf.list()
    print(sf[0].alllabels())
    print(dir(sf[0]))
