#/*##########################################################################
# Copyright (C) 2004 - 2013 European Synchrotron Radiation Facility
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
__author__ = "E. Papillon & V.A. Sole - ESRF Software Group"
import sys
if sys.version < '3.0':
    import ConfigParser
else:
    import configparser as ConfigParser
try:
    import numpy
    USE_NUMPY = True
except ImportError:
    # do not use numpy, use lists
    USE_NUMPY = False


class ConfigDict(dict):
    def __init__(self, defaultdict=None, initdict=None, filelist=None):
        if defaultdict is None:
            defaultdict = {}
        dict.__init__(self, defaultdict)
        self.default = defaultdict
        self.filelist = []

        if initdict is not None:
            self.update(initdict)
        if filelist is not None:
            self.read(filelist)

    def reset(self):
        """ Revert to default values
        """
        self.clear()
        self.update(self.default)

    def clear(self):
        """ Clear dictionnary
        """
        dict.clear(self)
        self.filelist = []

    def _check(self):
        pass

    def __tolist(self, mylist):
        if mylist is None:
            return None
        if not isinstance(mylist, list):
            return [mylist]
        else:
            return mylist

    def getfiles(self):
        return self.filelist

    def getlastfile(self):
        return self.filelist[len(self.filelist) - 1]

    def __convert(self, option):
        return option

    def read(self, filelist, sections=None):
        """
        read the input filename into the internal dictionary
        """
        filelist = self.__tolist(filelist)
        sections = self.__tolist(sections)
        cfg = ConfigParser.ConfigParser()
        cfg.optionxform = self.__convert
        cfg.read(filelist)
        self.__read(cfg, sections)

        for ffile in filelist:
            self.filelist.append([ffile, sections])
        self._check()

    def __read(self, cfg, sections=None):
        cfgsect = cfg.sections()

        if sections is None:
            readsect = cfgsect
        else:
            readsect = [sect for sect in cfgsect if sect in sections]

        for sect in readsect:
            ddict = self
            for subsectw in sect.split('.'):
                subsect = subsectw.replace("_|_", ".")
                if not (subsect in ddict):
                    ddict[subsect] = {}
                ddict = ddict[subsect]
            for opt in cfg.options(sect):
                ddict[opt] = self.__parse_data(cfg.get(sect, opt))

    def __parse_data(self, data):
        if len(data):
            if data.find(',') == -1:
                # it is not a list
                if USE_NUMPY and (data[0] == '[') and (data[-1] == ']'):
                    # this looks as an array
                    try:
                        return numpy.array([float(x) for x in data[1:-1].split()])
                    except ValueError:
                        try:
                            if (data[2] == '[') and (data[-3] == ']'):
                                nrows = len(data[3:-3].split('] ['))
                                indata = data[3:-3].replace('] [', ' ')
                                indata = numpy.array([float(x) for x in
                                                      indata.split()])
                                indata.shape = nrows, -1
                                return indata
                        except ValueError:
                            pass
        dataline = [line for line in data.splitlines()]
        if len(dataline) == 1:
            return self.__parse_line(dataline[0])
        else:
            return [self.__parse_line(line) for line in dataline]

    def __parse_line(self, line):
        if line.find(',') != -1:
            if line.endswith(','):
                if ',' in line[:-1]:
                    return [self.__parse_string(sstr.strip())
                            for sstr in line[:-1].split(',')]
                else:
                    return [self.__parse_string(line[:-1].strip())]
            else:
                return [self.__parse_string(sstr.strip())
                        for sstr in line.split(',')]
        else:
            return self.__parse_string(line.strip())

    def __parse_string(self, sstr):
        try:
            return int(sstr)
        except ValueError:
            try:
                return float(sstr)
            except ValueError:
                return sstr

    def tostring(self, sections=None):
        import StringIO
        tmp = StringIO.StringIO()
        sections = self.__tolist(sections)
        self.__write(tmp, self, sections)
        return tmp.getvalue()

    def write(self, filename, sections=None):
        """
        Write the current dictionary to the given filename
        """
        sections = self.__tolist(sections)
        fp = open(filename, "w")
        self.__write(fp, self, sections)
        fp.close()

    def __write(self, fp, ddict, sections=None, secthead=None):
        dictkey = []
        listkey = []
        valkey = []
        for key in ddict.keys():
            if isinstance(ddict[key], list):
                listkey.append(key)
            elif hasattr(ddict[key], 'keys'):
                dictkey.append(key)
            else:
                valkey.append(key)

        for key in valkey:
            if USE_NUMPY:
                if isinstance(ddict[key], numpy.ndarray):
                    fp.write('%s =' % key + ' [ ' +
                             ' '.join([str(val) for val in ddict[key]]) +
                             ' ]\n')
                    continue
            fp.write('%s = %s\n' % (key, ddict[key]))

        for key in listkey:
            fp.write('%s = ' % key)
            llist = []
            sep = ', '
            for item in ddict[key]:
                if isinstance(item, list):
                    if len(item) == 1:
                        llist.append('%s,' % item[0])
                    else:
                        llist.append(', '.join([str(val) for val in item]))
                    sep = '\n\t'
                else:
                    llist.append(str(item))
            fp.write('%s\n' % (sep.join(llist)))
        if 0:
            # this optimization method does not pass the tests.
            # disable it for the time being.
            if sections is not None:
                dictkey= [ key for key in dictkey if key in sections ]
        for key in dictkey:
            if secthead is None:
                newsecthead = key.replace(".", "_|_")
            else:
                newsecthead = '%s.%s' % (secthead, key.replace(".", "_|_"))
            #print "newsecthead = ", newsecthead
            fp.write('\n[%s]\n' % newsecthead)
            self.__write(fp, ddict[key], key, newsecthead)


def prtdict(ddict, lvl=0):
    for key in ddict.keys():
        if hasattr(ddict[key], 'keys'):
            print('\t' * lvl),
            print('+', key)
            prtdict(ddict[key], lvl + 1)
        else:
            print('\t' * lvl),
            print('-', key, '=', ddict[key])


def main():
    if len(sys.argv) > 1:
        config = ConfigDict(filelist=sys.argv[1:])
        prtdict(config)
    else:
        print("USAGE: %s <filelist>" % sys.argv[0])


if __name__ == '__main__':
    main()
