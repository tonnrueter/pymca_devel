#/*##########################################################################
# Copyright (C) 2004-2013 European Synchrotron Radiation Facility
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
__revision__ = "$Revision: 1.26 $"
__author__ = "V.A. Sole - ESRF Software Group"
import sys
import copy
import numpy
from PyMca import Elements
from PyMca.XRFMC import XRFMCHelper

class ConcentrationsConversion(object):
    def getConcentrationsAsHtml(self, concentrations=None):
        text = ""
        if concentrations is None:
            return text

        result = concentrations
        #the header
        if 'mmolar' in result:
            mmolarflaglist = [False, True]
        else:
            mmolarflaglist = [False]

        for mmolarflag in mmolarflaglist:
            text += "\n"
            text += "<H2><a NAME=""%s""></a><FONT color=#009999>" %\
                'Concentrations'
            if mmolarflag:
                text += "%s:" % 'mM Concentrations'
            else:
                text += "%s:" % 'Concentrations'
            text += "</FONT></H2>"
            text += "<br>"
            labels = ['Element', 'Group', 'Fit Area', 'Sigma Area']
            if mmolarflag:
                labels += 'mM concentration'
            else:
                labels += 'Mass fraction'

            #the table
            if 'layerlist' in result:
                if type(result['layerlist']) != type([]):
                    result['layerlist'] = [result['layerlist']]
                for label in result['layerlist']:
                    labels += [label]
            lemmon = ("#%x%x%x" % (255, 250, 205)).upper()
            white = '#FFFFFF'
            hcolor = ("#%x%x%x" % (230, 240, 249)).upper()
            text += "<CENTER>"
            text += "<nobr>"
            text += '<table width="80%" border="0" cellspacing="1" cellpadding="1" >'
            text += "<tr>"
            for l in range(len(labels)):
                if l < 2:
                    text += '<td align="left" bgcolor=%s><b>%s</b></td>' %\
                        (hcolor, labels[l])
                elif l == 2:
                    text += '<td align="center" bgcolor=%s><b>%s</b></td>' %\
                        (hcolor, labels[l])
                else:
                    text += '<td align="right" bgcolor=%s><b>%s</b></td>' %\
                        (hcolor, labels[l])
            text += "</tr>"
            line = 0
            for group in result['groups']:
                text += ("<tr>")
                element, group0 = group.split()
                fitarea    = "%.6e" % result['fitarea'][group]
                sigmaarea  = "%.2e" % result['sigmaarea'][group]
                area       = "%.6e" % result['area'][group]
                if mmolarflag:
                    fraction   = "%.4g" % result['mmolar'][group]
                else:
                    fraction   = "%.4g" % result['mass fraction'][group]
                if 'Expected Area' in labels:
                    fields = [element, group0, fitarea, sigmaarea, area, fraction]
                else:
                    fields = [element, group0, fitarea, sigmaarea, fraction]
                if 'layerlist' in result:
                    for layer in result['layerlist']:
                        if result[layer]['mass fraction'][group] < 0.0:
                            fraction   = "Unknown"
                        else:
                            if mmolarflag:
                                fraction = "%.4g" % result[layer]['mmolar'][group]
                            else:
                                fraction = "%.4g" % result[layer]['mass fraction'][group]
                        fields += [fraction]
                if line % 2:
                    color = lemmon
                else:
                    color = white
                i = 0
                for field in fields:
                    if (i<2):
                        text += '<td align="left"  bgcolor=%s>%s</td>' % (color, field)
                    else:
                        text += '<td align="right" bgcolor=%s>%s</td>' % (color, field)
                    i += 1
                text += '</tr>'
                line += 1
            text += ("</table>")
            text += ("</nobr>")
            text += "</CENTER>"
        return text

    def getConcentrationsAsAscii(self, concentrations=None):
        text = ""
        if concentrations is None:
            return text
        result = concentrations
        #the table
        if 'mmolar' in result:
            mmolarflaglist = [False, True]
        else:
            mmolarflaglist = [False]
        for mmolarflag in mmolarflaglist:
            labels = ['Element', 'Group', 'Fit_Area', 'Sigma_Area']
            if mmolarflag:
                labels += 'mM_Concentration'
            else:
                labels += 'Mass_fraction'
            if 'layerlist' in result:
                if type(result['layerlist']) != type([]):
                    result['layerlist'] = [result['layerlist']]
                for label in result['layerlist']:
                    labels += [label.replace(' ', '')]
            for l in labels:
                text += "%s  " % l
            text += ("\n")
            for group in result['groups']:
                element, group0 = group.split()
                fitarea = "%.6e" % result['fitarea'][group]
                sigmaarea = "%.2e" % result['sigmaarea'][group]
                area = "%.6e" % result['area'][group]
                if mmolarflag:
                    fraction = "%.4g" % result['mmolar'][group]
                else:
                    fraction = "%.4g" % result['mass fraction'][group]
                if 'Expected Area' in labels:
                    fields = [element, group0, fitarea, sigmaarea, area,
                              fraction]
                else:
                    fields = [element, group0, fitarea, sigmaarea, fraction]
                if 'layerlist' in result:
                    for layer in result['layerlist']:
                        if result[layer]['mass fraction'][group] < 0.0:
                            fraction   = "Unknown"
                        else:
                            if mmolarflag:
                                fraction = "%.4g" %\
                                    result[layer]['mmolar'][group]
                            else:
                                fraction = "%.4g" %\
                                    result[layer]['mass fraction'][group]
                        fields += [fraction]
                for field in fields:
                    text += '%s  ' % (field)
                text += '\n'
        return text


class ConcentrationsTool(object):
    def __init__(self, config=None, fitresult=None):
        self.config = {}
        self.config['usematrix'] = 0
        self.config['useattenuators'] = 1
        self.config['usemultilayersecondary'] = 0
        self.config['usexrfmc'] = 0
        self.config['flux'] = 1.0E10
        self.config['time'] = 1.0
        self.config['area'] = 30.0
        self.config['distance'] = 10.0
        self.config['reference'] = "Auto"
        self.config['mmolarflag'] = 0
        if config is not None:
            self.configure(config)
        self.fitresult = fitresult

    def configure(self, ddict=None):
        if ddict is None:
            ddict = {}
        for key in ddict:
            if key in self.config.keys():
                self.config[key] = ddict[key]
        return copy.deepcopy(self.config)

    def processFitResult(self, config=None, fitresult=None,
                         elementsfrommatrix=False, fluorates=None,
                         addinfo=False):
        # I should check if fit was successful ...
        if fitresult is None:
            fitresult = self.fitresult
        else:
            self.fitresult = fitresult
        if config is None:
            config = self.config
        else:
            self.config = config
        if 'usemultilayersecondary' not in self.config:
            self.config['usemultilayersecondary'] = 0
        if 'usexrfmc' not in self.config:
            self.config['usexrfmc'] = 0
        secondary = self.config['usemultilayersecondary']
        xrfmcSecondary = self.config['usexrfmc']
        if secondary and xrfmcSecondary:
            txt = "Only one of built-in secondary and Monte Carlo correction can be used"
            raise ValueError(txt) 
        # get attenuators and matrix from fit
        attenuators = []
        beamfilters = []
        funnyfilters = []
        matrix = None
        detectoratt = None
        multilayer = None
        for attenuator in fitresult['result']['config']['attenuators'].keys():
            if not fitresult['result']['config']['attenuators'][attenuator][0]:
                continue
            if attenuator.upper() == "MATRIX":
                matrix = fitresult['result']['config']['attenuators'][attenuator][1:4]
                alphain  = fitresult['result']['config']['attenuators'][attenuator][4]
                alphaout = fitresult['result']['config']['attenuators'][attenuator][5]
            elif attenuator.upper()[:-1] == "BEAMFILTER":
                beamfilters.append(fitresult['result']['config']['attenuators']\
                                                                [attenuator][1:])
            elif attenuator.upper() == "DETECTOR":
                detectoratt = fitresult['result']['config']['attenuators'][attenuator][1:]
            else:
                if len(fitresult['result']['config']['attenuators'][attenuator][1:]) == 4:
                   fitresult['result']['config']['attenuators'][attenuator].append(1.0)
                if abs(fitresult['result']['config']['attenuators'][attenuator][4]-1.0) > 1.0e-10:
                    #funny attenuator
                    funnyfilters.append(fitresult['result']['config']['attenuators']\
                                                                [attenuator][1:])
                else:
                    attenuators.append(fitresult['result']['config']['attenuators']\
                                                                [attenuator][1:])
        if matrix is None:
            raise ValueError("Invalid or undefined sample matrix")

        if matrix[0].upper() == "MULTILAYER":
            layerlist = fitresult['result']['config']['multilayer'].keys()
            layerlist.sort()
            for layer in layerlist:
                if fitresult['result']['config']['multilayer'][layer][0]:
                    if multilayer is None:
                        multilayer = []
                    multilayer.append(fitresult['result']['config']['multilayer'][layer][1:])
                    if not Elements.isValidMaterial(multilayer[-1][0]):
                        raise ValueError("Material %s is not defined" % multilayer[-1][0])

        else:
            layerlist = ["Layer0"]
            multilayer = [matrix]
            if not Elements.isValidMaterial(matrix[0]):
                raise ValueError("Material %s is not defined" % matrix[0])
        if xrfmcSecondary and (len(layerlist) > 1):
            txt = "Multilayer Monte Carlo correction not implemented yet"
            raise ValueError(txt)
        energyList = fitresult['result']['config']['fit']['energy']
        if energyList is None:
            raise ValueError("Invalid energy")
        if type(energyList) != type([]):
            energyList    = [energyList]
            flagList   = [1]
            weightList = [1.0]
        else:
            flagList   = fitresult['result']['config']['fit']['energyflag']
            weightList = fitresult['result']['config']['fit']['energyweight']
        finalEnergy = []
        finalWeight = []
        finalFlag = []
        for idx in range(len(energyList)):
            if flagList[idx]:
                energy = energyList[idx]
                if energy is None:
                    raise ValueError(\
                          "Energy %d isn't a valid energy" % idx)
                if energy <= 0.001:
                    raise ValueError(\
                          "Energy %d with value %f isn't a valid energy" %\
                          (idx, energy))
                if weightList[idx] is None:
                    raise ValueError(\
                          "Weight %d isn't a valid weight" % idx)
                if weightList[idx] < 0.0:
                    raise ValueError(\
                          "Weight %d with value %f isn't a valid weight" %\
                          (idx, weightList[idx]))
                finalEnergy.append(energy)
                finalWeight.append(weightList[idx])
                finalFlag.append(1)
        totalWeight = sum(weightList)
        if totalWeight == 0.0:
            raise ValueError("Sum of energy weights is 0.0")
        weightList = [x / totalWeight for x in finalWeight]
        energyList = finalEnergy
        flagList   = finalFlag

        # get elements list from fit, not from matrix
        groupsList = fitresult['result']['groups'] * 1
        if type(groupsList) != type([]):
            groupsList = [groupsList]

        todelete = []
        for i in range(len(groupsList)):
            ele = groupsList[i].split()[0]
            if len(ele) > 2:
                todelete.append(i)
        if len(todelete):
            todelete.reverse()
            for i in todelete:
                del groupsList[i]

        elements = []
        newelements = []
        for group in groupsList:
            splitted = group.split()
            ele = splitted[0]
            newelements.append([Elements.getz(splitted[0]),
                                splitted[0], splitted[1]])
            if len(elements):
                if elements[-1] != ele:
                    elements.append(ele)
            else:
                elements.append(ele)
        newelements.sort()
        elements.sort()
        if not config['useattenuators']:
            attenuators  = None
            funnyfilters = None
        #import time
        #t0=time.time()
        if elementsfrommatrix:
            newelementsList = []
            for ilayer in range(len(multilayer)):
                pseudomatrix = multilayer[ilayer]
                eleDict = Elements.getMaterialMassFractions([pseudomatrix[0]], [1.0])
                if eleDict == {}:
                    raise ValueError(\
                        "Invalid layer material %s" % pseudomatrix[0])
                keys = eleDict.keys()
                for ele in keys:
                    for group in newelements:
                        if ele == group[1]:
                            if not group in newelementsList:
                                newelementsList.append(group)
            newelementsList.sort()
            fluo0 = Elements.getMultilayerFluorescence(multilayer,
                         energyList,
                         layerList=None,
                         weightList=weightList,
                         flagList=weightList,
                         fulloutput=1,
                         beamfilters=beamfilters * 1,
                         attenuators=attenuators * 1,
                         elementsList=newelementsList * 1,
                         alphain=alphain,
                         alphaout=alphaout,
                         cascade=True,
                         detector=detectoratt,
                         funnyfilters=funnyfilters * 1,
                         forcepresent=0,
                         secondary=secondary)
            fluototal = fluo0[0]
            fluolist = fluo0[1:]
        else:
            if matrix[0].upper() != "MULTILAYER":
                multilayer = [matrix * 1]
            if fluorates is None:
                fluo0 = Elements.getMultilayerFluorescence(multilayer,
                             energyList,
                             layerList=None,
                             weightList=weightList,
                             flagList=flagList,
                             fulloutput=1,
                             beamfilters=beamfilters * 1,
                             attenuators=attenuators * 1,
                             elementsList=newelements * 1,
                             alphain=alphain,
                             alphaout=alphaout,
                             cascade=True,
                             detector=detectoratt,
                             funnyfilters=funnyfilters * 1,
                             forcepresent=1,
                             secondary=secondary)
            else:
                fluo0 = fluorates
            fluototal = fluo0[0]
            fluolist = fluo0[1:]
        #I'll need total fluo element by element at some point
        #print "getMatrixFluorescence elapsed = ",time.time()-t0
        if config['usematrix']:
            present = []
            referenceLayerDict = {}
            materialComposition = []
            for ilayer in range(len(multilayer)):
                pseudomatrix = multilayer[ilayer] * 1
                #get elemental composition from matrix
                materialComposition.append(Elements.getMaterialMassFractions([pseudomatrix[0]], [1.0]))
                keys = materialComposition[-1].keys()
                materialElements = [[Elements.getz(x), x] for x in keys]
                materialElements.sort()
                for z, key in materialElements:
                    for ele in elements:
                        if key == ele:
                            present.append(key)
                            if not (ele in referenceLayerDict):
                                referenceLayerDict[ele] = []
                            referenceLayerDict[ele].append(ilayer)
            if len(present) == 0:
                text = "Matrix must contain at least one fitted element\n"
                text += "in order to estimate flux and efficiency from it."
                raise ValueError(text)
            referenceElement = config['reference'].replace(' ', "")
            if len(referenceElement) and (referenceElement.upper() != 'AUTO'):
                if Elements.isValidFormula(referenceElement):
                    if len(referenceElement) == 2:
                        referenceElement = referenceElement.upper()[0] +\
                                           referenceElement.lower()[1]
                    elif len(referenceElement) == 1:
                        referenceElement = referenceElement.upper()[0]
                    if not (referenceElement in elements):
                        text = "Element %s not among fitted elements" % referenceElement
                        raise ValueError(text)
                    elif not (referenceElement in present):
                        text = "Element %s not among matrix elements" % referenceElement
                        raise ValueError(text)
                    referenceLayers  = referenceLayerDict[referenceElement]
                else:
                    text = "Element %s not a valid element" % referenceElement
                    raise ValueError(text)
            elif len(present) == 1:
                referenceElement = present[0]
                referenceLayers  = referenceLayerDict[referenceElement]
            else:
                # how to choose? Best fitted, largest fit area or
                # greater concentration?  or better to give a weight to
                # the different shells, energies , ...?
                referenceElement = present[0]
                fom = self._figureOfMerit(present[0],fluototal,fitresult)
                for key in present:
                    #if materialComposition[key] > materialComposition[referenceElement]:
                    #    referenceElement = key
                    newfom = self._figureOfMerit(key,fluototal,fitresult)
                    if newfom > fom:
                        fom = newfom
                        referenceElement = key
                referenceLayers  = referenceLayerDict[referenceElement]

            referenceTransitions = None
            for group in groupsList:
                item = group.split()
                element = item[0]
                if element == referenceElement:
                    transitions = item[1] + " xrays"
                    if referenceTransitions is None:
                        referenceTransitions = transitions
                        referenceGroup = group
                    elif (referenceTransitions[0] == transitions[0]) and\
                         (referenceTransitions[0] == 'L'):
                        # this prevents selecting L1 and selects L3 although
                        # given the appropriate area, L2 can be a safer choice.
                        referenceGroup = group
                        referenceTransitions = transitions
                elif referenceTransitions is not None:
                    break
            theoretical = 0.0
            for ilayer in referenceLayers:
                if elementsfrommatrix:
                    theoretical += fluolist[ilayer][referenceElement]['rates'][referenceTransitions] * \
                                   fluolist[ilayer][referenceElement]['mass fraction']
                else:
                    theoretical  += materialComposition[ilayer][referenceElement] * \
                                    fluolist[ilayer][referenceElement]['rates'][referenceTransitions]
            if theoretical <= 0.0:
                raise ValueError(\
                    "Theoretical rate is almost 0.0 Impossible to determine flux")
            else:
                if (config['distance'] > 0.0) and (config['area'] > 0.0):
                    #solidangle = config['area']/(4.0 * numpy.pi * pow(config['distance'],2))
                    radius2 = config['area']/numpy.pi
                    solidangle = 0.5 * (1.0 -  (config['distance']/numpy.sqrt(pow(config['distance'],2)+ radius2)))
                else:
                    solidangle = 1.0
                flux = fitresult['result'][referenceGroup]['fitarea'] / (theoretical * solidangle)
        else:
            referenceElement = None
            referenceTransitions = None
            #solidangle = config['area']/(4.0 * numpy.pi * pow(config['distance'],2))
            radius2 = config['area']/numpy.pi
            solidangle = 0.5 * (1.0 -  (config['distance']/numpy.sqrt(pow(config['distance'],2)+ radius2)))
            flux       = config['flux'] * config['time']

        #print "OBTAINED FLUX * SOLID ANGLE= ",flux * solidangle
        #print "flux * time = ",flux
        #print "actual solid angle = ",0.5 * (1.0 -  (config['distance']/sqrt(pow(config['distance'],2)+ config['area']/pi)))
        #print "solid angle factor= ",solidangle
        #ele  = 'Pb'
        #rays = "L xrays"
        #print "theoretical = ",fluototal[ele]['rates'][rays]
        #print "expected    = ",flux * solidangle * fluototal[ele]['rates'][rays]

        #for ilayer in range(len(multilayer)):
        #    print "ilayer = ",ilayer, "theoretical = ",fluolist[ilayer][ele]['rates'][rays]
        #    print "ilayer = ",ilayer, "expected = ",flux * solidangle * fluolist[ilayer][ele]['rates'][rays]
        ddict = {}
        ddict['groups'] = groupsList
        ddict['elements'] = elements
        ddict['mass fraction'] = {}
        if 'mmolarflag' in config:
            if config['mmolarflag']:
                ddict['mmolar'] = {}
        else:
            config['mmolarflag'] = 0
        ddict['area'] = {}
        ddict['fitarea'] = {}
        ddict['sigmaarea'] = {}
        fluo = fluototal

        for group in groupsList:
            item = group.split()
            element = item[0]
            transitions = item[1] + " xrays"
            if element in fluo.keys():
                if transitions in fluo[element]:
                    #this SHOULD be with concentration one
                    theoretical = fluo[element]['rates'][transitions] * 1.0
                    expected = theoretical * flux * solidangle
                    concentration = fitresult['result'][group]['fitarea']/expected
                else:
                    theoretical   = 0.0
                    concentration = 0.0
            else:
                theoretical   = 0.0
                concentration = 0.0
            #ddict['area'][group]    = theoretical * flux * solidangle * concentration
            ddict['fitarea'][group] = fitresult['result'][group]['fitarea']
            ddict['sigmaarea'][group] = fitresult['result'][group]['sigmaarea']
            if elementsfrommatrix:
                if element in fluo.keys():
                    ddict['mass fraction'][group] = 1.0 * fluo[element]['mass fraction']
                else:
                    ddict['mass fraction'][group] = 0.0
                ddict['area'][group]    = theoretical * flux * solidangle *\
                                         ddict['mass fraction'][group]
            else:
                ddict['mass fraction'][group] = concentration
                ddict['area'][group]    = theoretical * flux * solidangle
            if config['mmolarflag']:
                #mM = (mass_fraction * density)/atomic_weight
                ddict['mmolar'] [group]= 1000000. *\
                                        (multilayer[0][1] * ddict['mass fraction'][group])/Elements.Element[element]['mass']

        #I have the globals/average values now I calculate layer per layer
        #if necessary
        ddict['layerlist'] = []
        if matrix[0].upper() == "MULTILAYER":
            ilayer = 0
            for layer in layerlist:
                if fitresult['result']['config']['multilayer'][layer][0]:
                    ddict['layerlist'].append(layer)
                    ddict[layer] = {}
                    dict2 = ddict[layer]
                    dict2['groups'] = groupsList
                    dict2['elements'] = elements
                    dict2['mass fraction'] = {}
                    if config['mmolarflag']:
                        dict2['mmolar'] = {}
                    dict2['area'] = {}
                    dict2['fitarea'] = {}
                    fluo = fluolist[ilayer]
                    for group in groupsList:
                        item = group.split()
                        element = item[0]
                        transitions = item[1] + " xrays"
                        if element in fluo.keys():
                            if transitions in fluo[element]:
                                theoretical = fluo[element]['rates'][transitions] * 1
                                expected = theoretical * flux * solidangle
                                if expected > 0.0:
                                    concentration = fitresult['result'][group]['fitarea']/expected
                                else:
                                    concentration = -1
                            else:
                                theoretical   = 0.0
                                concentration = 0.0
                        else:
                            theoretical   = 0.0
                            concentration = 0.0
                        dict2['fitarea'][group] = 1 * fitresult['result'][group]['fitarea']
                        if elementsfrommatrix:
                            if element in fluo.keys():
                                dict2['mass fraction'][group] = 1 * fluo[element]['mass fraction']
                            else:
                                dict2['mass fraction'][group] = 0.0
                            #I calculate matrix in optimized form,
                            #so I have to multiply by the mass fraction
                            dict2['area'][group] = theoretical * flux * solidangle *\
                                dict2['mass fraction'][group]
                        else:
                            dict2['mass fraction'][group] = concentration
                            dict2['area'][group] = theoretical * flux * solidangle
                        if config['mmolarflag']:
                            #mM = (mass_fraction * density)/atomic_weight
                            dict2['mmolar'][group] = 1000000. *\
                                (multilayer[ilayer][1] * dict2['mass fraction'][group]) /\
                                Elements.Element[element]['mass']
                        #if group == "Pb L":
                        #    print "layer", ilayer,'area ', dict2['area'][group]
                        #    print "layer", ilayer,'mass fraction =', dict2['mass fraction'][group]
                    ilayer += 1
            if elementsfrommatrix:
                for group in groupsList:
                    ddict['area'][group] = 0.0
                    for layer in ddict['layerlist']:
                        if group in ddict[layer]['area'].keys():
                            ddict['area'][group] += ddict[layer]['area'][group]
        if xrfmcSecondary and (not elementsfrommatrix):
            xrfmcCorrections = None
            if 'xrfmc' in fitresult:
                xrfmcCorrections = fitresult['xrfmc'].get('corrections', None)
            if xrfmcCorrections is None:
                if 'xrfmc' in fitresult['result']:
                    xrfmcCorrections = fitresult['result']['xrfmc'].get('corrections', None)
            if xrfmcCorrections is None:
                # try to see if they were in the configuration
                if 'xrfmc' in fitresult['result']['config']:
                    xrfmcCorrections = fitresult['result']['config']['xrfmc'].get('corrections',
                                                                                  None)
            if xrfmcCorrections is None:
                # calculate the corrections
                xrfmcCorrections = XRFMCHelper.getXRFMCCorrectionFactors(fitresult['result']['config'])

            
            if referenceElement is not None:
                referenceLines = referenceTransitions.split()[0]
                referenceCorrection = xrfmcCorrections[referenceElement][referenceLines]\
                                            ['correction_factor'][-1]
                xrfmcCorrections[referenceElement][referenceLines]\
                                            ['correction_factor'][-1] = 1.0
                for group in groupsList:
                    item = group.split()
                    element = item[0]
                    lines = item[1]
                    if element in xrfmcCorrections:
                        if element != referenceElement:
                            if lines != referenceLines:
                                xrfmcCorrections[element][lines]['correction_factor'][-1] *= referenceCorrection
            # now we have to apply the corrections
            for group in groupsList:
                item = group.split()
                element = item[0]
                lines = item[1]
                if element in xrfmcCorrections:
                    correction = xrfmcCorrections[element][item[1]]['correction_factor'][-1]
                    ddict['mass fraction'][group] /= correction
        if addinfo:
            addInfo = {}
            addInfo['ReferenceElement'] = referenceElement
            addInfo['ReferenceTransitions'] = referenceTransitions
            addInfo['SolidAngle'] = solidangle
            if config['time'] > 0.0:
                addInfo['Time'] = config['time']
            else:
                addInfo['Time'] = 1.0
            addInfo['Flux'] = flux / addInfo['Time']       
            addInfo['I0'] = flux
            addInfo['DetectorDistance'] = config['distance']
            addInfo['DetectorArea'] = config['area']
            return ddict , addInfo
        else:
            return ddict

    def _figureOfMerit(self, element, fluo, fitresult):
        weight = 0.0
        for transitions in fluo[element]['rates'].keys():
            if fluo[element]['rates'][transitions] > 0.0:
                if   (transitions[0] == "K") and (Elements.getz(element) > 18):
                    factor = 2.0
                elif (transitions[0] == "L") and (Elements.getz(element) > 54):
                    factor = 1.5
                else:
                    factor = 1.0
                group = element + " " + transitions.split()[0]
                if group in fitresult['result']['groups']:
                    fitarea = fitresult['result'][group]['fitarea']
                    weightHelp = fitarea * fluo[element]['rates'][transitions] * factor * \
                                fluo[element]['mass fraction']
                    if weightHelp > weight:
                        weight = weightHelp
        return weight

def main():
    import sys
    import getopt

    from PyMca import ConfigDict

    if len(sys.argv) > 1:
        options = ''
        longoptions = ['flux=', 'time=', 'area=', 'distance=',
                       'attenuators=', 'usematrix=']
        tool = ConcentrationsTool()
        opts, args = getopt.getopt(
                        sys.argv[1:],
                        options,
                        longoptions)
        config = tool.configure()
        for opt, arg in opts:
            if opt in ('--flux'):
                config['flux'] = float(arg)
            elif opt in ('--area'):
                config['area'] = float(arg)
            elif opt in ('--time'):
                config['time'] = float(arg)
            elif opt in ('--distance'):
                config['distance'] = float(arg)
            elif opt in ('--attenuators'):
                config['useattenuators'] = int(float(arg))
            elif opt in ('--usematrix'):
                config['usematrix'] = int(float(arg))
        tool.configure(config)
        filelist = args
        for filename in filelist:
            d = ConfigDict.ConfigDict()
            d.read(filename)
            for material in d['result']['config']['materials'].keys():
                Elements.Material[material] =\
                    copy.deepcopy(d['result']['config']['materials'][material])
            print(tool.processFitResult(fitresult=d, elementsfrommatrix=True))
    else:
        print("Usage:")
        print("ConcentrationsTool [--flux=xxxx --area=xxxx] fitresultfile")

if __name__ == "__main__":
    main()
