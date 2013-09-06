#!/usr/bin/env python
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
import numpy
from PyMca import EdfFile

DEBUG = 0

class PyMcaBatchBuildOutput(object):
    def __init__(self, inputdir=None, outputdir=None):
        self.inputDir  = inputdir
        self.outputDir = outputdir

    def buildOutput(self, inputdir=None, outputdir=None, delete=None):
        if inputdir is None:inputdir = self.inputDir
        if inputdir is None:inputdir = os.getcwd()
        if outputdir is None: outputdir = self.outputDir
        if outputdir is None: outputdir = inputdir
        if delete is None:
            if outputdir == inputdir:
                delete = True
        if DEBUG:
            print("delete option = ", delete)
        allfiles = os.listdir(inputdir)
        partialedflist = []
        partialdatlist = []
        partialconlist = []
        for filename in allfiles:
            if filename.endswith('000000_partial.edf'):partialedflist.append(filename)
            elif filename.endswith('000000_partial.dat'):partialdatlist.append(filename)
            elif filename.endswith('000000_partial_concentrations.txt'):partialconlist.append(filename)

        #IMAGES
        edfoutlist = []
        for filename in partialedflist:
            if DEBUG:
                print("Dealing with filename %s" % filename)
            edflist = self.getIndexedFileList(os.path.join(inputdir, filename))
            i = 0
            for edfname in edflist:
                edf    = EdfFile.EdfFile(edfname, access='rb', fastedf = 0)
                nImages = edf.GetNumImages()
                #get always the last image
                data0   = edf.GetData(nImages-1)
                data0[data0<0] = 0
                if i == 0:
                    header = edf.GetHeader(0)
                    data = data0.copy()
                else:
                    data += data0
                del edf
                i += 1
            edfname  = filename.replace('_000000_partial.edf',".edf")
            edfoutname = os.path.join(outputdir, edfname)
            if DEBUG:
                print("Dealing with output filename %s" % edfoutname)
            if os.path.exists(edfoutname):
                if DEBUG:
                    print("Output file already exists, trying to delete it")
                os.remove(edfoutname)
            edfout   = EdfFile.EdfFile(edfoutname, access="wb")
            edfout.WriteImage (header , data, Append=0)
            del edfout
            edfoutlist.append(edfoutname)
            if delete:
                for filename in edflist:
                    try:
                        os.remove(filename)
                    except:
                        print("Cannot delete file %s" % filename)

        #DAT IMAGES
        datoutlist = []
        for filename in partialdatlist:
            edflist = self.getIndexedFileList(os.path.join(inputdir, filename))
            first = True
            for edfname in edflist:
                f = open(edfname)
                lines = f.readlines()
                f.close()
                j = 1
                while (not len( lines[-j].replace("\n",""))):
                       j += 1
                if first:
                    first = False
                    labels = lines[0].replace("\n","").split("  ")
                    nlabels = len(labels)
                    nrows = len(lines) - j

                    data      = numpy.zeros((nrows, nlabels), numpy.double)
                    inputdata = numpy.zeros((nrows, nlabels), numpy.double)
                chisqIndex = labels.index('chisq')
                for i in range(nrows):
                    inputdata[i, :] = [float(x) for x in lines[i+1].split()]
                    if inputdata[i, chisqIndex] < 0.0:
                        inputdata[i, chisqIndex] = 0.0
                data += inputdata
            outfilename = os.path.join(outputdir, filename.replace("_000000_partial",""))
            if os.path.exists(outfilename):
                os.remove(outfilename)
            outfile=open(outfilename,'w+')
            outfile.write('%s' % lines[0])
            line=""
            for row in range(nrows):
                #line = "%d" % inputdata[row, 0]
                for col in range(nlabels):
                    if   col == 0 : line += "%d" % inputdata[row, col]
                    elif   col == 1 : line += "  %d" % inputdata[row, col]
                    else: line += "  %g" % data[row, col]
                line += "\n"
                outfile.write("%s" % line)
                line =""
            outfile.write("\n") 
            outfile.close()
            datoutlist.append(outfilename)
            if delete:
                for filename in edflist:
                    os.remove(filename)


        #CONCENTRATIONS
        outconlist = []
        for filename in partialconlist:
            edflist = self.getIndexedFileList(os.path.join(inputdir, filename))
            i = 0
            for edfname in edflist:
                edf    = open(edfname, 'rb')
                if i == 0:
                    outfilename = os.path.join(outputdir, filename.replace("_000000_partial",""))
                    if os.path.exists(outfilename):
                        os.remove(outfilename)
                    outfile = open(outfilename,'wb')
                lines = edf.readlines()
                for line in lines:
                    outfile.write(line)
                edf.close()
                i += 1
            outfile.close()
            outconlist.append(outfilename)
            if delete:
                for filename in edflist:
                    os.remove(filename)
        return edfoutlist, datoutlist, outconlist
        
    def getIndexedFileList(self, filename, begin=None,end=None, skip = None, fileindex=0):
        name = os.path.basename(filename)
        n = len(name)
        i = 1
        numbers = ['0', '1', '2', '3', '4', '5',
                   '6', '7', '8','9']
        while (i <= n):
            c = name[n-i:n-i+1]
            if c in ['0', '1', '2',
                                '3', '4', '5',
                                '6', '7', '8',
                                '9']:
                break
            i += 1
        suffix = name[n-i+1:]
        if len(name) == len(suffix):
            #just one file, one should use standard widget
            #and not this one.
            self.loadFileList(filename, fileindex=fileindex)
        else:
            nchain = []
            while (i<=n):
                c = name[n-i:n-i+1]
                if c not in ['0', '1', '2',
                                    '3', '4', '5',
                                    '6', '7', '8',
                                    '9']:
                    break
                else:
                    nchain.append(c)
                i += 1
            number = ""
            nchain.reverse()
            for c in nchain:
                number += c
            fformat = "%" + "0%dd" % len(number)
            if (len(number) + len(suffix)) == len(name):
                prefix = ""
            else:
                prefix = name[0:n-i+1]
            prefix = os.path.join(os.path.dirname(filename),prefix)
            if not os.path.exists(prefix + number + suffix):
                print("Internal error in EDFStack")
                print("file should exist: %s " % (prefix + number + suffix))
                return
            i = 0
            if begin is None:
                begin = 0
                testname = prefix+fformat % begin+suffix
                while not os.path.exists(prefix+ fformat % begin+suffix):
                    begin += 1
                    testname = prefix+fformat % begin+suffix
                    if len(testname) > len(filename):break
                i = begin
            else:
                i = begin
            if not os.path.exists(prefix+fformat % i+suffix):
                raise ValueError("Invalid start index file = %s" % \
                      (prefix+fformat % i+suffix))
            f = prefix+fformat % i+suffix
            filelist = []
            while os.path.exists(f):
                filelist.append(f)
                i += 1
                if end is not None:
                    if i > end:
                        break
                f = prefix+fformat % i+suffix
            return filelist


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage:")
        print("python PyMcaBatchBuildOutput.py directory")
        sys.exit(0)
    directory = sys.argv[1]
    w = PyMcaBatchBuildOutput(directory)
    w.buildOutput()
    """
    allfiles = os.listdir(directory)
    edflist = []
    datlist = []
    for filename in allfiles:
        if filename.endswith('000000_partial.edf'):edflist.append(filename)
        elif filename.endswith('000000_partial.dat'):datlist.append(filename)
    for filename in edflist:
        print w.getIndexedFileList(os.path.join(directory, filename))
    """
