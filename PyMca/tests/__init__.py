import os
if os.path.exists('PyMca'):
    if os.path.exists('setup.py'):
        if os.path.exists('py2app_setup.py'):
            txt ='Tests cannnot be imported from top source directory'
            raise ImportError(txt)
from PyMca.tests.TestAll import main as testAll
from PyMca.tests.ConfigDictTest import test as testConfigDict
from PyMca.tests.EdfFileTest import test as testEdfFile
from PyMca.tests.ElementsTest import test as testElements
from PyMca.tests.GefitTest import test as testGefit
from PyMca.tests.PCAToolsTest import test as testPCATools
from PyMca.tests.SpecfileTest import test as testSpecfile
from PyMca.tests.specfilewrapperTest import test as testSpecfilewrapper
