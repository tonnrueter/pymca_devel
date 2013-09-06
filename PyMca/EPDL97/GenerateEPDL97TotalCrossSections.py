__doc__= "Generate specfile from EPL97 total cross sections in keV and barn" 
import os
import sys
import EPDL97Parser as EPDLParser
Elements = EPDLParser.Elements
AVOGADRO_NUMBER = EPDLParser.AVOGADRO_NUMBER
import numpy
log = numpy.log
exp = numpy.exp
getTotalCoherentCrossSection = EPDLParser.getTotalCoherentCrossSection
getTotalIncoherentCrossSection = EPDLParser.getTotalIncoherentCrossSection
getTotalPhotoelectricCrossSection = EPDLParser.getTotalPhotoelectricCrossSection
getTotalPairCrossSection = EPDLParser.getTotalPairCrossSection
getTotalTripletCrossSection = EPDLParser.getTotalTripletCrossSection

if len(sys.argv) < 3:
    print("Usage:")
    print("python EPDLGenerateTotalCrossSections SPEC_output_filename barns_flag")
    sys.exit(0)

def getHeader(filename):
    text  = '#F %s\n' % filename
    text += '#U00 This file is a direct conversion to specfile format of \n'
    text += '#U01 the original EPDL97 total cross sections contained in the\n'
    text += '#U02 EPDL97.DAT from the library.\n'
    text += '#U03 EPDL97 itself can be found at:\n'
    text += '#U04           http://www-nds.iaea.org/epdl97/libsall.htm\n'
    text += '\n'
    return text

fname = sys.argv[1]
if os.path.exists(fname):
    os.remove(fname)

if int(sys.argv[2]):
    BARNS = True
else:
    BARNS = False
print("BARNS = %s" % BARNS)
outfile = open(fname, 'wb')
outfile.write(getHeader(fname))

for i in range(1, 101):
    print("i = %d element = %s" % (i, Elements[i-1]))
    #coherent
    energy_cohe, value_cohe, mode_cohe = getTotalCoherentCrossSection(i,
                                                            getmode=True)
    #incoherent
    energy_incohe, value_incohe, mode_incohe = getTotalIncoherentCrossSection(i,
                                                            getmode=True)

    #photoelectric
    energy_photo, value_photo, mode_photo = getTotalPhotoelectricCrossSection(i,
                                                            getmode=True)

    #check to see the energies:
    #for j in range(10):
    #    print energy_cohe[j], energy_incohe[j], energy_photo[j]


    #to select an appropriate energy grid as close as possible to the original
    #while keeping in mind the PyMca goals, I use the coherent energy grid till
    #the non-zero first value of the photoelectric cross section. At that point,
    #I use the photoelectric energy grid.
    energy = numpy.concatenate((energy_cohe[energy_cohe<energy_photo[0]],
                               energy_photo))

    #now perform a log-log interpolation when needed
    #lin-lin interpolation:
    #
    #              y0 (x1-x) + y1 (x-x0)
    #        y = -------------------------
    #                     x1 - x0
    #
    #log-log interpolation:
    #
    #                  log(y0) * log(x1/x) + log(y1) * log(x/x0)
    #        log(y) = ------------------------------------------
    #                                  log (x1/x0)
    #
    cohe    = numpy.zeros(len(energy), numpy.float)
    incohe  = numpy.zeros(len(energy), numpy.float)
    photo   = numpy.zeros(len(energy), numpy.float)
    total   = numpy.zeros(len(energy), numpy.float)    

    #coherent needs to interpolate
    indices = numpy.nonzero(energy_cohe<energy_photo[0]) 
    cohe[indices]  = value_cohe[indices]
    for n in range(len(indices),len(energy)):
        x = energy[n]
        j1 = len(indices)
        while energy_cohe[j1] < x:
            j1 += 1
        j0 = j1 - 1
        x0 = energy_cohe[j0]
        x1 = energy_cohe[j1]
        y0 = value_cohe[j0]
        y1 = value_cohe[j1]
        cohe[n] = exp((log(y0) * log(x1/x) + log(y1) * log(x/x0))/log(x1/x0))

    #compton needs to interpolate everything
    for n in range(len(energy)):
        x = energy[n]
        j1 = 0
        while energy_incohe[j1] < x:
            j1 += 1
        j0 = j1 - 1
        x0 = energy_incohe[j0]
        x1 = energy_incohe[j1]
        y0 = value_incohe[j0]
        y1 = value_incohe[j1]
        incohe[n] = exp((log(y0) * log(x1/x) + log(y1) * log(x/x0))/log(x1/x0))

    #photoelectric does not need to interpolate anything
    photo[energy>=energy_photo[0]] = value_photo[:]


    #convert to keV and cut at 500 keV
    energy *= 1000.
    indices = numpy.nonzero(energy<=500.)
    energy = energy[indices]
    photo  = photo[indices]
    cohe   = cohe[indices]
    incohe = incohe[indices]

    #I cut at 500 keV, I do not need to take the pair production
    total = photo + cohe + incohe

    #now I am ready to write a Specfile
    ele = Elements[i-1]
    text  = '#S %d %s\n' % (i, ele)
    text += '#N 5\n'
    labels = '#L PhotonEnergy[keV]'
    labels += '  Rayleigh(coherent)[barn/atom]' 
    labels += '  Compton(incoherent)[barn/atom]'
    labels += '  CoherentPlusIncoherent[barn/atom]'
    labels += '  Photoelectric[barn/atom]'
    labels += '  TotalCrossSection[barn/atom]\n'
    if not BARNS:
        labels = labels.replace("barn/atom", "cm2/g")
        factor = (1.0E-24*AVOGADRO_NUMBER)/EPDLParser.getAtomicWeights()[i-1]
    else:
        factor = 1.0
    text += labels
    if 0:
        fformat = "%g %g %g %g %g %g\n"
    else:
        fformat = "%.6E %.6E %.6E %.6E %.6E %.6E\n"
    outfile.write(text)
    for n in range(len(energy)):
        line = fformat % (energy[n],
                          cohe[n] * factor,
                          incohe[n] * factor,
                          (cohe[n]+incohe[n]) * factor,
                          photo[n] * factor,
                          total[n] * factor)
        outfile.write(line)
outfile.write('\n')
outfile.close()
