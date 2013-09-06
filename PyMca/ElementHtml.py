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
from PyMca import Elements

class ElementHtml(object):
    def __init__(self,element=None):
        self.element = None
        
    def gethtml(self,element=None):
        if element is None:element = self.element
        if element is None:return ""
        ele = element
        #text="<center><b><font color=red size=5>Summary</font></b></center>"
        text=""
        if ele not in Elements.Element.keys(): 
            text+="<br><b><font color=blue size=4>Unknown Element</font></b>"
            return text
        symbol = Elements.getsymbol(Elements.getz(ele))
        omegak = Elements.Element[ele]['omegak']
        omegal = [  Elements.Element[ele]['omegal1'],
                    Elements.Element[ele]['omegal2'],
                    Elements.Element[ele]['omegal3'] ]
        omegam = [  Elements.Element[ele]['omegam1'],
                    Elements.Element[ele]['omegam2'],
                    Elements.Element[ele]['omegam3'],
                    Elements.Element[ele]['omegam4'],
                    Elements.Element[ele]['omegam5'] ]
        
        #text+="<center>
        text+="<br><b><font color=blue size=4>Element Info</font></b>"
        #text+="</center>"
        if 0:
            text+="<br><b><font size=3>Name = %s</font></b>"        % Elements.Element[ele]['name']
            text+="<br><b><font size=3>Symbol = %s</font></b>"      % symbol
            text+="<br><b><font size=3>At. Number = %d</font></b>"  % Elements.Element[ele]['Z']
            text+="<br><b><font size=3>At. Weight = %.5f</font></b>"  % Elements.Element[ele]['mass']
            text+="<br><b><font size=3>Density = %.5f</font></b>"     % Elements.Element[ele]['density']
        else:
            hcolor = 'white'
            finalcolor = 'white'
            

            text+="<nobr><table>"
            #symbol
            text+="<tr>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>Symbol</font></b>"
            text+="</td>"
            text+='<td align="center" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>=</font></b>"
            text+="</td>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>%s </font></b>"  % symbol
            text+="</td>"            
            #Z
            text+="<tr>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>At. Number</font></b>"
            text+="</td>"
            text+='<td align="center" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>=</font></b>"
            text+="</td>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>%d </font></b>"  % Elements.Element[ele]['Z']
            text+="</td>"                        
            #name
            text+="<tr>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>Name</font></b>"
            text+="</td>"
            text+='<td align="center" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>=</font></b>"
            text+="</td>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            name = Elements.Element[ele]['name'][0].upper()+Elements.Element[ele]['name'][1:]
            text+="<b><font size=3>%s </font></b>"  % name
            text+="</td>"
            #mass
            text+="<tr>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>At. Weight</font></b>"
            text+="</td>"
            text+='<td align="center" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>=</font></b>"
            text+="</td>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>%.5f </font></b>"  % Elements.Element[ele]['mass']
            text+="</td>"                        
            #density
            text+="<tr>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>Density</font></b>"
            text+="</td>"
            text+='<td align="center" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>=</font></b>"
            text+="</td>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>%.5f g/cm3</font></b>"  % Elements.Element[ele]['density']
            text+="</td>"
            text+="</tr>"
            text+="</table>"


        # Shell propierties
        hcolor = 'white'
        finalcolor = 'white'
        if Elements.Element[ele]['Z'] > 2:
            text+="<br><b><font color=blue size=4>Fluorescence Yields</font></b>"
            text+="<nobr><table><tr>"
            text+='<td align="left" bgcolor="%s"><b>' % hcolor
            text+='Shell'
            text+="</b></td>"
            text+='<td align="right" bgcolor="%s"><b>' % hcolor
            text+='Yield'
            text+="</b></td>"
            text+="</tr>"
            text+="<tr>"
            text+='<td align="left" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>%s </font></b>"  % "K"
            text+="</td>"
            text+='<td align="right" bgcolor="%s">' % finalcolor
            text+="<b><font size=3>%.3e </font></b>"  % omegak
            text+="</td>"
        for i in range(len(omegal)):
            if omegal[i] > 0.0:
                    text+="<tr>"
                    text+='<td align="left" bgcolor="%s">' % finalcolor
                    text+="<b><font size=3>L%d </font></b>"  % (i+1)
                    text+="</td>"
                    text+='<td align="right" bgcolor="%s">' % finalcolor
                    text+="<b><font size=3>%.3e </font></b>"  % omegal[i]
                    text+="</td>"
        for i in range(len(omegam)):
            if omegam[i] > 0.0:
                    text+="<tr>"
                    text+='<td align="left" bgcolor="%s">' % finalcolor
                    text+="<b><font size=3>M%d </font></b>"  % (i+1)
                    text+="</td>"
                    text+='<td align="right" bgcolor="%s">' % finalcolor
                    text+="<b><font size=3>%.3e </font></b>"  % omegam[i]
                    text+="</td>"
        text+="</tr>"
        text+="</table>"

        hcolor = 'white'
        finalcolor = 'white'
        f  = ['f12','f13','f23']
        ck = []
        doit = 0
        for item in f:
            value = Elements.Element[ele]['CosterKronig']['L'][item]
            if value > 0:doit=1
            ck.append(value)
        if doit:
            text+="<br><b><font color=blue size=4>L-Shell Coster-Kronig</font></b>"
            text+="<nobr><table><tr>"
            for item in f:
                text+='<td align="left" bgcolor="%s"><b>' % hcolor
                text+=item
                text+="</b></td>"
            text+="</tr>"
            text+="<tr>"
            for i in range(len(f)):
                text+='<td align="left" bgcolor="%s">' % finalcolor
                text+="<b><font size=3>%.3f </font></b>"  % ck[i]
            text+="</td>"
            text+="</tr>"
            text+="</table>"
            #M shell
            fs = [[ 'f12', 'f13', 'f14', 'f15'],
                          ['f23', 'f24', 'f25'],
                                 ['f34', 'f35'],
                                        ['f45']]
            doit = 0
            for f in fs:
                for item in f:
                    value = Elements.Element[ele]['CosterKronig']['M'][item]
                    if value > 0:doit=1
            if doit:
                text+="<br><b><font color=blue size=4>M-Shell Coster-Kronig</font></b>"
                text+="<nobr><table>"
                for f in fs:
                    text+="<tr>"
                    for item in f:
                        text+='<td align="left" bgcolor="%s"><b>' % hcolor
                        text+=item
                        text+="</b></td>"
                    text+="</tr>"
                    text+="<tr>"
                    for item in f:
                        text+='<td align="left" bgcolor="%s">' % finalcolor
                        text+="<b><font size=3>%.3f </font></b>"  % Elements.Element[ele]['CosterKronig']['M'][item]
                    text+="</td>"
                    text+="</tr>"
                text+="</table>"
        
        hcolor = 'white'
        finalcolor = 'white'
        for rays in Elements.Element[ele]['rays']:
            if rays == "Ka xrays":continue
            if rays == "Kb xrays":continue
            #text+="<center>"
            text+="<br><b><font color=blue size=4>%s Emission Energies</font></b>" % rays[0:-1]
            #text+="</center>"
            if 0:
                for transition in Elements.Element[ele][rays]:
                    text+="<br><b><font size=3>%s energy = %.5f  rate = %.5f</font></b>"  % (transition,Elements.Element[ele][transition]['energy'],
                                                                            Elements.Element[ele][transition]['rate'])

            else:
                text+="<nobr><table><tr>"
                text+='<td align="left" bgcolor="%s"><b>' % hcolor
                text+='Line'
                text+="</b></td>"
                text+='<td align="right" bgcolor="%s"><b>' % hcolor
                text+='Energy (keV)'
                text+="</b></td>"
                text+='<td align="right" bgcolor="%s"><b>' % hcolor
                text+='Rate'
                text+="</b></td>"
                text+="</tr>"
                for transition in Elements.Element[ele][rays]:
                    transitiontext = transition.replace('*','')
                    text+="<tr>"
                    text+='<td align="left" bgcolor="%s">' % finalcolor
                    text+="<b><font size=3>%s </font></b>"  % transitiontext
                    text+="</td>"
                    text+='<td align="right" bgcolor="%s">' % finalcolor
                    text+="<b><font size=3>%.5f</font></b>"  % Elements.Element[ele][transition]['energy']
                    text+="</td>"
                    text+='<td align="right" bgcolor="%s">' % finalcolor
                    text+="<b><font size=3>%.5f </font></b>"  % Elements.Element[ele][transition]['rate']
                    text+="</td>"
                text+="</tr>"
                text+="</table>"

        hcolor = 'white'
        finalcolor = 'white'
        #text+="<center>"
        text+="<br><b><font color=blue size=4>%s Binding Energies</font></b>" % "Electron"
        #text+="</center>"
        text+="<nobr><table><tr>"
        text+='<td align="left" bgcolor="%s"><b>' % hcolor
        text+='Shell'
        text+="</b></td>"
        text+='<td align="right" bgcolor="%s"><b>' % hcolor
        text+='Energy (keV)'
        text+="</b></td>"
        text+="</tr>"
        for shell in Elements.ElementShells:
            if Elements.Element[ele]['binding'][shell] > 0.0:
                text+="<tr>"
                text+='<td align="left" bgcolor="%s">' % finalcolor
                text+="<b><font size=3>%s </font></b>"  % shell
                text+="</td>"
                text+='<td align="right" bgcolor="%s">' % finalcolor
                text+="<b><font size=3>%.5f </font></b>"  % Elements.Element[ele]['binding'][shell]
                text+="</td>"
        text+="</tr>"
        text+="</table>"
        return text
        
        
if __name__ == "__main__":
    #
    # WARNING: This test method uses PyQt and therefore its license is GPL
    # If you use the methods of this module and not this test portion, you
    # can safely use it as LPGL for HTML generation.
    #
    import sys
    from PyMca import PyMcaQt as qt
    app  = qt.QApplication(sys.argv)
    if len(sys.argv) > 1:
        ele = sys.argv[1]
    else:
        ele = "Fe"
    w= qt.QWidget()
    l=qt.QVBoxLayout(w)
    html = ElementHtml()
    if qt.qVersion() < '4.0.0':
        text = qt.QTextEdit(w)
        text.setText(html.gethtml(ele))
        text.setReadOnly(1)
        app.setMainWidget(w)
    else:
        text = qt.QTextEdit(w)
        text.insertHtml(html.gethtml(ele))
        text.setReadOnly(1)
    l.addWidget(text)
    w.show()
    if qt.qVersion() < '4.0.0':
        app.exec_loop()
    else:
        app.exec_()
