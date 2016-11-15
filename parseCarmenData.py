#!/usr/bin/env python

import sys
import os
import math

from CarmenParser import CarmenParser

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

    for i in range(len(scans)):
        scan = scans[i]

        thetaLen = 10
        x,y,theta = scan[0:3]
        dx = []
        dy = []
        lrData = scan[4]
        stepSize = 0.5 if (len(lrData) == 360 or len(lrData) == 361) else 1.0
        for i in range(len(lrData)):
            dx.append(x + (lrData[i] * math.cos(math.radians(theta + 90.0 - stepSize*float(i)))))
            dy.append(y + (lrData[i] * math.sin(math.radians(theta + 90.0 - stepSize*float(i)))))

        # plot robot location and heading
        #plt.plot([x],[y], 'ro')
        #plt.plot([x, x+(thetaLen*math.cos(math.radians(theta)))],[y, y+(thetaLen*math.sin(math.radians(theta)))], 'r-')

        # plot a smaller robot location marker than above (and no heading)
        plt.plot([x],[y], 'r.')

        # plot laser scan data points
        #plt.plot(dx, dy, 'b.')

    plt.xlabel('x pos')
    plt.ylabel('y pos')
    plt.savefig('scan.png')

if __name__ == '__main__':
    if len(sys.argv) == 2:
        fpin = sys.argv[1]
        print(fpin)

        parser = CarmenParser()
        parser.parse(fpin)

        plotPos(parser.posData)

        plotScan(parser.rangeData)
