"""
IntelParser class

Implementation of DataParser for Intel building data set
http://cres.usc.edu/radishrepository/view-one.php?name=intel_lab
"""

from DataParser import DataParser
import re

class IntelParser(DataParser):
    def parse(self, file_in):
        fin = open(file_in, 'r')
        for l in fin:
            line = l.rstrip()
            if line.startswith('POS'):
                _, s, us, f1, f2, f3 = line.split(' ')
                # TODO: why is the angle 90 - f3?
                # this is how the rec2script.c file handles it
                x,y,o = float(f1), float(f2), (90.0 - float(f3))
                self.posData.append((x,y,o))
            elif line.startswith('LASER-RANGE'):
                header, rangeDataStr = line.split(':')
                _, s, us, lasernr, numValues, arange = (re.sub(' +', ' ', header.strip())).split()
                rangeDataStr = re.sub(' +', ' ', rangeDataStr.strip())
                self.rangeData.append([float(x) for x in rangeDataStr.split(' ')])
