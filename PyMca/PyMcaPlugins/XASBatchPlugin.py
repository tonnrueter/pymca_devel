try:
    from PyMca import Plugin1DBase
except ImportError:
    from . import Plugin1DBase

try:
    from PyMca import XASBatch
except ImportError:
    print("XMCDWindow importing from somewhere else")
    import XASBatch

class XASBatchPlugin(Plugin1DBase.Plugin1DBase):
    def __init__(self, plotWindow, **kw):
        Plugin1DBase.Plugin1DBase.__init__(self,  plotWindow,  **kw)
        self.methodDict = {
            'showDialog': [self.showDialog, 'Displays the XASBatchDialog', None]
        }
        self.widget = None
        
    #Methods to be implemented by the plugin
    def getMethods(self, plottype=None):
        """
        A list with the NAMES  associated to the callable methods
        that are applicable to the specified plot.

        Plot type can be "SCAN", "MCA", None, ...        
        """
        names = list(self.methodDict.keys())
        names.sort()
        return names

    def getMethodToolTip(self, name):
        """
        Returns the help associated to the particular method name or None.
        """
        return self.methodDict[name][1]

    def getMethodPixmap(self, name):
        """
        Returns the pixmap associated to the particular method name or None.
        """
        return self.methodDict[name][2]

    def applyMethod(self, name):
        """
        The plugin is asked to apply the method associated to name.
        """
        self.methodDict[name][0]()
        return

    def showDialog(self):
        if self.widget == None:
            self.widget = XASBatch.XASBatchDialog(None)
        self.widget.show()
        self.widget.raise_()

        
MENU_TEXT = "XAS Something"
def getPlugin1DInstance(plotWindow, **kw):
    ob = XASBatchPlugin(plotWindow)
    return ob

def main():
    from PyMca import PyMcaQt as qt
    from PyMca import Plot1D
    app = qt.QApplication([])
    plot = Plot1D.Plot1D()
    plugin = getPlugin1DInstance(plot)
    plugin.applyMethod(plugin.getMethods()[0])
    print app.exec_()
    
if __name__ == '__main__':
    main()

