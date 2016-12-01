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

def plotPos(pos):
    """
    Plot robot position from odometry data
    """
    x = [p[0] for p in pos]
    y = [p[1] for p in pos]
    plt.clf()
    plt.plot(x, y, 'k.')
    plt.axis('equal')
    plt.xlabel('x pos')
    plt.ylabel('y pos')
    plt.savefig('pos.png')


def plotScan(scans):
    """
    Plot laser range data scan taken with robot at specified position.

    This plots position data stored with the laser scan data

    Commented out lines are for plotting the laser scan data points
    """

    plt.clf()

    global_sdf = SDF.asGlobal(initextents=(200.0, 200.0), res=0.1)

    for scan in scans:
        local_sdf = SDF.fromRangeLine(scan, res=0.1, threshold=10)
        global_sdf.addLocalSDF(local_sdf)

        # plot robot location and heading
        # (x,y,theta,_,lrData) = scan
        # thetaLen = 0.5
        # brgs = [brg for (brg,_) in lrData]
        # maxBrg = max(brgs)
        # minBrg = min(brgs)
        # plt.plot([x],[y], 'bo')
        # plt.plot([x, x+(2*thetaLen*math.cos(theta))],
        #          [y, y+(2*thetaLen*math.sin(theta))],
        #          'b-')
        # plt.plot([x, x+(thetaLen*math.cos(theta + maxBrg))],
        #          [y, y+(thetaLen*math.sin(theta + maxBrg))],
        #          'r-')
        # plt.plot([x, x+(thetaLen*math.cos(theta + minBrg))],
        #          [y, y+(thetaLen*math.sin(theta + minBrg))],
        #          'g-')


    global_sdf.plot(plt)

    # plt.axis('equal')
    # plt.xlabel('x pos')
    # plt.ylabel('y pos')
    # plt.show()
    # plt.savefig('scan.png')

if __name__ == '__main__':
    if len(sys.argv) == 2:
        fpin = sys.argv[1]
        print(fpin)

        parser = CarmenParser()
        parser.parse(fpin)

        plotPos(parser.posData)

        plotScan(parser.rangeData)
