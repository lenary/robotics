#!/usr/bin/env python

import sys
import os
import math

import numpy as np

from CarmenParser import CarmenParser
from ICP import ICP
from SDF import SDF

# check for DISPLAY environment variable and use matplotlib 'Agg' if not present
# this is a check for Windows Linux Subsystem GUI issues
if os.environ.get('DISPLAY') is None:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

def plotPoints(points, c='k'):
    px,py = zip(*points)
    plt.plot(list(px),list(py),'{0}.'.format(c))

    plt.axis('equal')
    plt.xlabel('x pos')
    plt.ylabel('y pos')
    plt.show()
    plt.savefig('testICP.png')

def plotSingleScan(scan, c='k', rngMax=10.0):
    # plot robot location and heading
    (x,y,theta,_,lrData) = scan

    print '(x,y,theta,deg) = ({0},{1},{2},{3})'.format(x,y,theta,np.rad2deg(theta))

    lrData = [(brg,rng) for (brg,rng) in lrData if rng < rngMax]

    thetaLen = 0.5
    brgs = [brg for (brg,_) in lrData]
    maxBrg = max(brgs)
    minBrg = min(brgs)
    plt.plot([x],[y], '{0}o'.format(c))
    plt.plot([x, x+(2*thetaLen*math.cos(theta))],
             [y, y+(2*thetaLen*math.sin(theta))],
             '{0}-'.format(c))
    """
    plt.plot([x, x+(thetaLen*math.cos(theta + maxBrg))],
             [y, y+(thetaLen*math.sin(theta + maxBrg))],
             'r-')
    plt.plot([x, x+(thetaLen*math.cos(theta + minBrg))],
             [y, y+(thetaLen*math.sin(theta + minBrg))],
             'g-')
    """

    points = [(x + rng*math.cos(theta+brg), y + rng*math.sin(theta+brg)) for (brg,rng) in lrData]

    plotPoints(points, c)

    DATA = np.zeros((2,len(points)))
    for i in range(len(points)):
        DATA[0,i] = points[i][0]
        DATA[1,i] = points[i][1]
    return DATA

def singleScanData(scan, rngMax=10.0):
    # plot robot location and heading
    (x,y,theta,_,lrData) = scan

    print '(x,y,theta,deg) = ({0},{1},{2},{3})'.format(x,y,theta,np.rad2deg(theta))

    lrData = [(brg,rng) for (brg,rng) in lrData if rng < rngMax]
    thetaLen = 0.5
    brgs = [brg for (brg,_) in lrData]
    maxBrg = max(brgs)
    minBrg = min(brgs)
    points = [(x + rng*math.cos(theta+brg), y + rng*math.sin(theta+brg)) for (brg,rng) in lrData]
    DATA = np.zeros((2,len(points)))
    for i in range(len(points)):
        DATA[0,i] = points[i][0]
        DATA[1,i] = points[i][1]
    return DATA

if __name__ == '__main__':
    if len(sys.argv) == 2:
        fpin = sys.argv[1]
        print(fpin)

        parser = CarmenParser()
        parser.parse(fpin)

        A = plotSingleScan(parser.rangeData[0])

        # ICP settings
        N = 10
        e = 1e-5
        k = 3

        for i in range(1,10):
            B = singleScanData(parser.rangeData[i])
            R,T = ICP.icp(A,B,N,e,k)
            Bp = (R.dot(B).T + T).T
            BpPoints = [(p[0],p[1]) for p in Bp.T]
            plotPoints(BpPoints)

        """
        A = plotSingleScan(parser.rangeData[12])
        B = plotSingleScan(parser.rangeData[13], 'r')

        # ICP settings
        N = 10
        e = 1e-5
        k = 3
        
        R,T = ICP.icp(A, B, N, e, k)

        print R
        print T

        #cA,cB,E,iA,iB = ICP.correspondingPoints(A,B,k,ICP.pwSSE)
        #ApPoints = [(p[0],p[1]) for p in cB.T]
        #plotPoints(ApPoints,'m')

        Ap = (R.dot(B).T + T).T
        ApPoints = [(p[0],p[1]) for p in Ap.T]
        #plotPoints(ApPoints,'m')
        """
        
