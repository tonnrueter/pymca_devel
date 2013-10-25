from PyMca import PyMcaQt as qt
from PyMca import PyMca_Icons as Icons
from PyMca import FitParam
import numpy

DEBUG = 0
class XASBatchDialog(FitParam.FitParamDialog):
    def __init__(self, parent=None):
        FitParam.FitParamDialog.__init__(self, 
                                         parent=parent,
                                         name="XASBatchWidget")
        # Remove, but does not delete the pages!
        self.fitparam.mainTab.clear() 
        self.fitparam.mainTab.addTab(
                    self.fitparam.tabAttenuators,
                    'ATTENUATORS')
        self.fitparam.mainTab.addTab(
                    self.fitparam.tabMul,
                    'MATRIX')
                    
def openDialog():
    app= qt.QApplication([])
    app.connect(app, qt.SIGNAL("lastWindowClosed()"), app.quit)
    wid=XASBatchDialog()
    ret = wid.exec_()
    if ret == qt.QDialog.Accepted:
        npar = wid.getParameters()
        print(npar)
        del wid
    app.quit()

if __name__=="__main__":
    openDialog()
