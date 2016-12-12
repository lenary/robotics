#!/usr/bin/env python

import sys
import os
import math

from CarmenParser import CarmenParser
from SDF import SDF

# check for DISPLAY environment variable and use matplotlib 'Agg' if not present
# this is a check for Windows Linux Subsystem GUI issues
if os.environ.get('DISPLAY') is None:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

def plotScan(scans, name):
    """
    Plot laser range data scan taken with robot at specified position.

    This plots position data stored with the laser scan data

    Commented out lines are for plotting the laser scan data points
    """

    res = 0.01

    plt.clf()

    global_sdf = SDF.asGlobal(initextents=(200.0, 200.0), res=res)

    for scan in scans:
        local_sdf = SDF.fromRangeLine(scan, res=res, threshold=10)
        global_sdf.addLocalSDF(local_sdf)

    print "Scanning Complete: %s" % name

    global_sdf.plotSurfaces(name)
    global_sdf.plot(name)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        fpin = sys.argv[1]
        print("Logfile: %s" % fpin)

        parser = CarmenParser()
        parser.parse(fpin)

        plotScan(parser.rangeData, fpin.rsplit("/",1)[0])
