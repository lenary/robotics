#!/usr/bin/env python

import sys
import os
import math

from IntelParser import IntelParser

# check for DISPLAY environment variable and use matplotlib 'Agg' if not present
# this is a check for Windows Linux Subsystem GUI issues
if os.environ.get('DISPLAY') is None:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

def plotPos(pos):
    x = [p[0] for p in pos]
    y = [p[1] for p in pos]
    plt.clf()
    plt.plot(x,y)
    plt.xlabel('x pos')
    plt.ylabel('y pos')
    plt.savefig('pos.png')


def plotScan(pos, scanData):
    """
    Plot laser range data scan taken with robot at specified position.
    Also plot a robot bearing line, and robot position

    Assumes scanData is 180 points, in 180 degree arc (0 deg to 179 deg)
    ** really should be 181 points, but data in file is only 180

    This only plots a single scan at a single position
    """
    thetaLen = 500
    x,y,theta = pos
    dx = []
    dy = []
    for i in range(len(scanData)):
        dx.append(x + (scanData[i] * math.cos(math.radians(theta + 90.0 - float(i)))))
        dy.append(y + (scanData[i] * math.sin(math.radians(theta + 90.0 - float(i)))))
    plt.clf()
    plt.plot([x],[y], 'ro')
    plt.plot([x, x+(thetaLen*math.cos(math.radians(theta)))],[y, y+(thetaLen*math.sin(math.radians(theta)))], 'r-')
    plt.plot(dx, dy, 'b.')
    plt.xlabel('x pos')
    plt.ylabel('y pos')
    plt.savefig('scans.png')

if __name__ == '__main__':
    if len(sys.argv) == 2:
        fpin = sys.argv[1]
        print(fpin)

        ip = IntelParser()
        ip.parse(fpin)

        plotPos(ip.posData)

        print('length of range data: {0}'.format(len(ip.rangeData[0])))
        print('x,y,theta = {0}, {1}, {2}'.format(ip.posData[0][0], ip.posData[0][1], ip.posData[0][2]))
        plotScan(ip.posData[0], ip.rangeData[0])
