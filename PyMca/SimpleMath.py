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
import numpy
try:
    from PyMca import SGModule
except ImportError:
    print("SimpleMath importing SGModule directly")
    import SGModule
    
class SimpleMath(object):
    def derivate(self,xdata,ydata, xlimits=None):
        x=numpy.array(xdata, copy=False, dtype=numpy.float)
        y=numpy.array(ydata, copy=False, dtype=numpy.float)
        if xlimits is not None:
            i1=numpy.nonzero((xdata>=xlimits[0])&\
                               (xdata<=xlimits[1]))[0]
            x=numpy.take(x,i1)
            y=numpy.take(y,i1)
        i1 = numpy.argsort(x)
        x=numpy.take(x,i1)
        y=numpy.take(y,i1)  
        deltax=x[1:] - x[:-1]
        i1=numpy.nonzero(abs(deltax)>0.0000001)[0]
        x=numpy.take(x, i1)
        y=numpy.take(y, i1)
        minDelta = deltax.min()
        xInter = numpy.arange(x[0]-minDelta,x[-1]+minDelta,minDelta)
        yInter = numpy.interp(xInter, x, y, left=y[0], right=y[-1])
        if len(yInter) > 499:
            npoints = 5
        else:
            npoints = 3
        degree = 1
        order = 1
        coeff = SGModule.calc_coeff(npoints, degree, order)
        N = int(numpy.size(coeff-1)/2)
        yInterPrime = numpy.convolve(yInter, coeff, mode='valid')/minDelta
        i1 = numpy.nonzero((x>=xInter[N+1]) & (x <= xInter[-N]))[0]
        x = numpy.take(x, i1)
        result = numpy.interp(x, xInter[(N+1):-N],
                              yInterPrime[1:],
                              left=yInterPrime[1],
                              right=yInterPrime[-1])
        return x, result

    def average(self,xdata0,ydata0):
        #check if all the x axis are identical (no interpolation needed)
        allthesamex=1
        x0=xdata0[0]
        for xaxis in xdata0:
            if len(x0) == len(xaxis):
                if numpy.alltrue(x0==xaxis):
                    pass
                else:
                    allthesamex=0
                    break
            else:
                allthesamex=0
                break

        if allthesamex:
            xdata=[]
            ydata=[]
            i=0
            for x0 in xdata0:
                x=numpy.array(x0)
                xdata.append(x)
                ydata.append(numpy.array(ydata0[i]))
                i=i+1
                
            finalx=numpy.array(x0)
            finalx=xdata0[0]
            finaly=numpy.zeros(finalx.shape ,numpy.float)
            i = 0
            for x0 in xdata0:
                finaly += ydata[i]
                i=i+1
        else:
            #sort the data
            xdata=[]
            ydata=[]
            i=0
            for x0 in xdata0:
                x=numpy.array(x0)
                i1=numpy.argsort(x)
                xdata.append(numpy.take(x,i1))
                ydata.append(numpy.take(numpy.array(ydata0[i]),i1))
                i=i+1         
            
            #get the max and the min x axis
            xmin=xdata[0][0]
            xmax=xdata[0][-1]
            for x in xdata:
                if xmin < x[0]:
                    xmin=x[0]
                if xmax > x[-1]:
                    xmax=x[-1]
            #take the data in between
            x=[]
            y=[]
            i=0
            minimumLength = len(xdata[0])
            for x0 in xdata:
                i1=numpy.nonzero((x0>=xmin) & (x0<=xmax))[0]
                x.append(numpy.take(x0,i1))
                y.append(numpy.take(numpy.array(ydata[i]),i1))
                if len(x0) < minimumLength:
                    minimumLength = len(x0)
                i=i+1

            if minimumLength < 2:
                raise ValueError("Not enough points to take a meaningfull average")
            #take as x axis the first
            finalx=x[0]
            for i in range(len(x)):
                if x[i][0] > finalx[0]:
                    finalx = x[i] 
            finaly=numpy.zeros(finalx.shape, numpy.float)
            j=-1
            allthesamex=0
            for p in range(len(finalx)):
              point=finalx[p] 
              i=0
              j=j+1
              try:            
                for x0 in x:
                    if allthesamex:
                        finaly[p]+=y[i][p]
                    else:
                        i1=max(numpy.nonzero(x0<=point)[0])
                        i2=min(numpy.nonzero(x0>=point)[0])
                        if i1 >= i2:
                            #take the point as it is
                            finaly[p]+=y[i][i1]
                        else:
                            #interpolation
                            A=(x0[i2]-point)/(x0[i2]-x0[i1])
                            B=1.-A
                            finaly[p]+=A*y[i][i1]+B*y[i][i2]
                    i=i+1
              except:
                break
        if allthesamex:
              finalx=finalx[0:]
              finaly=finaly[0:]/len(xdata0)      
        else:
              finalx=finalx[0:j]
              finaly=finaly[0:j]/len(xdata0)
     
        return finalx,finaly

    def smooth(self, *var, **kw):
        """
        smooth(self,*vars,**kw)
        Usage: self.smooth(y)
               self.smooth(y=y)
               self.smooth()
        """
        if 'y' in kw:
            ydata=kw['y']
        elif len(var) > 0:
            ydata=var[0]
        else:
            ydata=self.y
        f=[0.25,0.5,0.25]
        result=numpy.array(ydata, copy=False, dtype=numpy.float)
        if len(result) > 1:
            result[1:-1]=numpy.convolve(result,f,mode=0)
            result[0]=0.5*(result[0]+result[1])
            result[-1]=0.5*(result[-1]+result[-2])
        return result
        
if __name__ == "__main__":
    x = numpy.arange(100.)*0.25
    y = x*x + 2 * x
    a = SimpleMath()
    #print(a.average(x,y))
    xplot, yprime = a.derivate(x, y)
    print("Found:")
    for i in range(0,10):
        print("x = %f  y'= %f expected = %f" % (xplot[i], yprime[i], 2*xplot[i]+2))

