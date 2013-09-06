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
import unittest
import sys
import os
import gc
import tempfile

class testSpecfile(unittest.TestCase):
    def setUp(self):
        """
        import the module
        """
        try:
            from PyMca import specfile
            self.specfileClass = specfile
        except:
            self.specfileClass = None
        if self.specfileClass is not None:
            text  = "#F \n"
            text += "\n"
            text += "#S 10  Undefined command 0\n"
            text += "#N 3\n"
            text += "#L First label  Second label  Third label\n"
            text += "10  100  1000\n"
            text += "20  400  8000\n"
            text += "30  900  270000\n"
            text += "\n"
            text += "#S 20  Undefined command 1\n"
            text += "#N 3\n"
            text += "#L First  Second  Third\n"
            text += "1.3  1  1\n"
            text += "2.5  4  8\n"
            text += "3.7  9  27\n"
            text += "\n"
            tmpFile = tempfile.mkstemp(text=False)
            if sys.version < '3.0':
                os.write(tmpFile[0], text)
            else:
                os.write(tmpFile[0], bytes(text, 'utf-8'))
            os.close(tmpFile[0])
            self.fname = tmpFile[1]

    def tearDown(self):
        """clean up any possible files"""
        # make sure the file handle is free
        self._sf = None
        self._scan = None
        # this should free the handle
        gc.collect()
        if self.specfileClass is not None:
            if os.path.exists(self.fname):
                os.remove(self.fname)

    def testSpecfileImport(self):
        #"""Test successful import"""
        self.assertTrue(self.specfileClass is not None,
                        'Unsuccessful PyMca.specfile import')

    def testSpecfileReading(self):
        #"""Test specfile readout"""
        self.testSpecfileImport()
        self._sf = self.specfileClass.Specfile(self.fname)
        # test the number of found scans
        self.assertEqual(len(self._sf), 2,
                         'Expected to read 2 scans, read %s' %\
                         len(self._sf))
        self.assertEqual(self._sf.scanno(), 2,
                         'Expected to read 2 scans, got %s' %\
                         self._sf.scanno())
        # test scan iteration selection method
        self._scan = self._sf[1]
        labels = self._scan.alllabels()
        expectedLabels = ['First', 'Second', 'Third']
        self.assertEqual(len(labels), 3,
                         'Expected to read 3 scans, got %s' % len(labels))
        for i in range(3):
            self.assertEqual(labels[i], expectedLabels[i],
                    'Read "%s" instead of "%s"' %\
                     (labels[i], expectedLabels[i]))
        # test scan number selection method
        self._scan = self._sf.select('20.1')
        labels = self._scan.alllabels()
        sf = None
        expectedLabels = ['First', 'Second', 'Third']
        self.assertEqual(len(labels), 3,
                         'Expected to read 3 labels, got %s' % len(labels))
        for i in range(3):
            self.assertEqual(labels[i], expectedLabels[i],
                'Read "%s" instead of "%s"' %\
                (labels[i], expectedLabels[i]))
        gc.collect()

    def testSpecfileReadingCompatibleWithUserLocale(self):
        #"""Test specfile compatible with C locale"""
        self.testSpecfileImport()
        self._sf = self.specfileClass.Specfile(self.fname)
        self._scan = self._sf[1]
        datacol = self._scan.datacol(1)
        data = self._scan.data()
        self._sf = None
        self.assertEqual(datacol[0], 1.3,
                    'Read %f instead of %f' %\
                    (datacol[0], 1.3))
        self.assertEqual(datacol[1], 2.5,
                    'Read %f instead of %f' %\
                    (datacol[1], 2.5))
        self.assertEqual(datacol[2], 3.7,
                    'Read %f instead of %f' %\
                    (datacol[2], 3.7))
        self.assertEqual(datacol[1], data[0][1],
                    'Read %f instead of %f' %\
                    (datacol[1], data[0][1]))
        gc.collect()

def getSuite(auto=True):
    testSuite = unittest.TestSuite()
    if auto:
        testSuite.addTest(\
            unittest.TestLoader().loadTestsFromTestCase(testSpecfile))
    else:
        # use a predefined order
        testSuite.addTest(testSpecfile("testSpecfileImport"))
        testSuite.addTest(testSpecfile("testSpecfileReading"))
        testSuite.addTest(\
            testSpecfile("testSpecfileReadingCompatibleWithUserLocale"))
    return testSuite

def test(auto=False):
    unittest.TextTestRunner(verbosity=2).run(getSuite(auto=auto))

if __name__ == '__main__':
    test()
