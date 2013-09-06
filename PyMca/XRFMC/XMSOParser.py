import sys
import os
import xml.etree.ElementTree as ElementTree

DEBUG = 0

def getXMSOFileFluorescenceInformation(xmsoFile):
    f = ElementTree.parse(xmsoFile)
    ddict = {}
    root = f.getroot()
    transitions = ['K', 'Ka', 'Kb', 'L', 'L1', 'L2', 'L3', 'M']
    for i in root.iter('fluorescence_line_counts'):
        if DEBUG:
            print(i.attrib)
            for key in ['symbol', 'total_counts']:
                print(key, '= ', i.get(key))
        element = i.get('symbol') 
        ddict[element] = {}
        #ddict[element]['z'] = i.get('atomic_number')
        for key in transitions:
            ddict[element][key] = { 'total':0.0,
                                    'counts': [],
                                    'correction_factor':[]}
        for a in i.iter('fluorescence_line'):
            if DEBUG:
                print(a.attrib)
                for key in ['type', 'total_counts']:
                    print(key, '= ', a.get(key))
            line = a.get('type')
            ddict[element][line] = {}
            #ddict[element][line]['total'] = float(a.get('total_counts'))
            ddict[element][line]['counts'] = []
            ddict[element][line]['total']=0
            transitionsAffected = []
            for key in transitions:
                if line.startswith(key):
                    transitionsAffected.append(key)
                elif line.startswith('KL') and (key == 'Ka'):
                    transitionsAffected.append(key)
                elif line.startswith('K') and (key == 'Kb'):
                    if not line.startswith('KL'):
                        transitionsAffected.append(key)
            cumulator = 0
            for b in a.iter('counts'):
                if DEBUG:
                    print(b.attrib)
                value = float(b.text)
                ddict[element][line]['counts'].append(value)
                cumulator += value
            ddict[element][line]['total'] = cumulator
            single = ddict[element][line]['counts'][0]
            multiple = 0.0
            ddict[element][line]['correction_factor'] = []
            excitationCounter = 0
            for value in ddict[element][line]['counts']:
                multiple += value
                ddict[element][line]['correction_factor'].append(\
                    multiple/single)
                for key in transitionsAffected:
                    nValues = len(ddict[element][line]['counts']) 
                    while(len(ddict[element][key]['counts']) < nValues):
                        ddict[element][key]['counts'].append(0.0)
                    ddict[element][key]['counts'][excitationCounter] += value
                excitationCounter += 1
        ddict[element][key]['correction_factor'] = []
        for key in transitions:
            multiple = 0.0
            if len(ddict[element][key]['counts']) == 0:
                nValues = len(ddict[element][line]['counts'])
                ddict[element][key]['counts'] = [0.0] * nValues
                ddict[element][key]['correction_factor'] = [1.0] * nValues
            else:
                single = ddict[element][key]['counts'][0]
                for value in ddict[element][key]['counts']:
                    multiple += value
                    ddict[element][key]['correction_factor'].append(\
                        multiple/single)
            ddict[element][key]['total'] = multiple
    return ddict

def test(xmsoFile='t.xmso'):
    ddict = getXMSOFileFluorescenceInformation(xmsoFile)
    for element in ddict:
        for line in ddict[element]:
            if line == "z":
                #atomic number
                continue
            if 1 or DEBUG or line in ['K', 'Ka', 'Kb', 'L', 'L1', 'L2', 'L3', 'M']:
                correction1 = ddict[element][line]['correction_factor'][1]
                correctionn = ddict[element][line]['correction_factor'][-1]
                print("Element %s Line %s Correction 2 = %f Correction n = %f" %\
                            (element, line,correction1, correctionn))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        if os.path.exists('t.xmso'):
            test()
        else:
            print("Usage:")
            print("python XMSOParser.py xmso_file") 
            sys.exit(0)
    else:
        test(sys.argv[1])
