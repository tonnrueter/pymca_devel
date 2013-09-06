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
import time
from PyMca import PyMcaLogo

class HtmlIndex(object):
    def __init__(self, htmldir):
        if htmldir is None:htmldir = "/tmp/HTML"
        self.htmldir = htmldir

    def isHtml(self, x):
        if len(x) < 5: return 0
        if x[-5:] == ".html":return 1

    def isHtmlDir(self, x):
        if len(x) < 7: return 0
        if x[-7:] == "HTMLDIR":return 1

    def getHeader(self,addlink=None):
        link= [ ['http://www.esrf.fr', 'ESRF home'],
                ['http://www.esrf.fr/computing/bliss/', 'BLISS home'],
          ]
        if addlink is not None:
            for item in addlink:
                link.append(item)
        text =""
        text+= "<HTML>"
        text+= "<HEAD>"
        text+= "<TITLE>PyMCA : Advanced Fit Results</TITLE>"
        text+= "</HEAD>"
        text+= "<BODY TEXT=#000000 BGCOLOR=#FFFFFF ALINK=#ff6600 LINK=#0000cc VLINK=#0000cc marginwidth=10 marginheight=10  topmargin=10 leftmargin=10>"
        text+= "<CENTER>"
        text+= "<TABLE WIDTH=100%% border=0 Height=70>"
        text+= "  <TR>"
        text+= "    <TD><Font Size=5 Color=#0000cc>"
        text+= "        <b>PyMCA : Advanced Fit Results</b></Font>"
        text+= "    </td>"
        text+= "    <td rowspan=2 ALIGN=RIGHT VALIGN=bottom>"
        text+= "        <a HREF=""http://www.esrf.fr/"">"
        logofile = self.htmldir + "/" + "PyMcaLogo.png"
        if not os.path.exists(logofile):
            try:
                import qt
                pixmap = qt.QPixmap(PyMcaLogo.PyMcaLogo)
                pixmap.save(logofile,"PNG")
            except:
                pass
        text+= "        <img SRC=%s ALT=""ESRF home"" WIDTH=55 HEIGHT=68 BORDER=0></a>" % "PyMcaLogo.png"
        text+= "    </td>"
        text+= "  </tr>"
        text+= "  <tr>"
        text+= "     <td width=100%%  VALIGN=bottom>"
        text+= "        <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0 WIDTH=100%%>"
        text+= "          <TR>"
        text+= "            <TD WIDTH=100%% BGCOLOR=#ee22aa HEIGHT=17  ALIGN=LEFT VALIGN=middle>"
        text+= "            <FONT color=#000000>&nbsp;"
        for name in link:
            text+= "|&nbsp;&nbsp;<A STYLE=""color: #FFFFFF"" HREF=""%s"">%s</a>&nbsp;&nbsp;"%(tuple(name))
        text+= "            </FONT>"
        text+= "            </TD>"
        text+= "          </TR>"
        text+= "        </TABLE>"
        text+= "     </td>"
        text+= "  </tr>"
        text+= "  <tr>"
        text+= "     <td colspan=2 height=5><spacer type=block height=10 width=0>"
        text+= "     </td>"
        text+= "  </tr>"
        text+= "</table>"
        text+= "</center>"
        return text

    def getFooter(self):
        now = time.time()
        text =""
        text+= "<center>"
        text+= "<table width=100%% border=0 cellspacing=0 cellpadding=0>"
        text+= "    <tr><td colspan=2 height=10><spacer type=block height=10 width=0></td></tr>"
        text+= "    <tr><td colspan=2 bgcolor=#cc0066 height=5><spacer type=block height=5 width=0></td></tr>"
        text+= "    <tr><td colspan=2 height=5><spacer type=block height=5 width=0></td></tr>"
        text+= "    <TR>"
        text+= "        <TD><FONT size=1 >created:  %s</font></TD>" % time.ctime(now)
        text+= "        <TD ALIGN=RIGHT><FONT size=1 >last modified: %s by" % time.ctime(now)
        #text+= "        <A STYLE=""color: #0000cc"" HREF=""mailto:papillon@esrf.fr"">papillon@esrf.fr</A></FONT></TD>"
        if sys.platform == 'win32':
            try:
                user = os.environ['USERNAME']
                text+= "        <A STYLE=""color: #0000cc"">%s</A></FONT></TD>" % user
            except:
                text +="</FONT></TD>"
        else:
            try:
                user = os.getlogin()
                text+= "        <A STYLE=""color: #0000cc"">%s</A></FONT></TD>" % user
            except:
                text +="</FONT></TD>"
        text+= "    </TR>"
        text+= "</TABLE>"
        text+= "</center>"
        text+= "</BODY>"
        text+= "</HTML>"
        return text

    def getBody(self, htmldir=None):
        if htmldir is None:htmldir = self.htmldir
        dirlist  = filter(self.isHtmlDir, os.listdir(htmldir))
        filelist = []
        for directory in dirlist:
            fulldir = os.path.join(self.htmldir,directory)
            filelist += filter(self.isHtml, os.listdir(fulldir))


        #I have a list of directories and of indexes
        for directory in dirlist:
            fulldir = os.path.join(htmldir,directory)
            index   = os.path.join(fulldir,"index.html")
            if os.path.exists(index):
                try:
                    os.remove(index)
                except:
                    print("cannot delete file %s" % index)
                    continue

    def _getHtmlFileList(self, directory):
        return filter(self.isHtml, os.listdir(directory))

    def _getHtmlDirList(self, directory):
        return filter(self.isHtmlDir, os.listdir(directory))

    def buildIndex(self, directory = None):
        if directory is None: directory = self.htmldir
        index = os.path.join(directory, "index.html")
        if os.path.exists(index):
            try:
                os.remove(index)
            except:
                print("cannot delete file %s" % index)
                return
        filelist = self._getHtmlFileList(directory)
        text = ""
        text += self.getHeader()
        for file in filelist:
            text +="<a href=""%s"">%s</a><BR>" % (file, file.split(".html")[0])
        text += self.getFooter()
        file=open(index,'wb')
        file.write(text)
        file.close()

    def buildRecursiveIndex(self, directory = None): 
        if directory is None: directory = self.htmldir
        index = os.path.join(directory, "index.html")
        if os.path.exists(index):
            try:
                os.remove(index)
            except:
                print("cannot delete file %s" % index)
                return
        directorylist = self._getHtmlDirList(directory)
        text = ""
        text += self.getHeader()
        for file in directorylist:
            fulldir = os.path.join(directory,file)
            self.buildIndex(fulldir)
            fileroot = file.split('_HTMLDIR')[0]
            link     = "./"+file+"/index.html"
            text +="<a href=""%s"">%s</a><BR>" % (link, fileroot)
        text += self.getFooter()
        file=open(index,'wb')
        file.write(text)
        file.close()
        
            
if __name__ == "__main__":
    if len(sys.argv) > 1:
        a = HtmlIndex(sys.argv[1])
    else:
        print("Trying /tmp/HTML as input directory")
        a = HtmlIndex('/tmp/HTML')
    a.buildRecursiveIndex()


