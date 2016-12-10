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

def plotPoints(points, c='k', fname='out.png'):
    px,py = zip(*points)
    plt.plot(list(px),list(py),'{0}.'.format(c))

    plt.axis('equal')
    plt.xlabel('x pos')
    plt.ylabel('y pos')
    plt.show()
    plt.savefig('{0}'.format(fname))

def plotRobot(x,y,theta,c='g',thetaLen=0.5):
    maxBrg = np.pi/2
    minBrg = -np.pi/2
    plt.plot([x],[y], '{0}o'.format(c))
    plt.plot([x, x+(2*thetaLen*math.cos(theta))],
             [y, y+(2*thetaLen*math.sin(theta))],
             '{0}-'.format(c))
    plt.plot([x, x+(thetaLen*math.cos(theta + maxBrg))],
             [y, y+(thetaLen*math.sin(theta + maxBrg))],
             'r-')
    plt.plot([x, x+(thetaLen*math.cos(theta + minBrg))],
             [y, y+(thetaLen*math.sin(theta + minBrg))],
             'g-')
    plt.show()

def plotSingleScan(scan, c='k', rngMax=10.0):
    # plot robot location and heading
    (x,y,theta,_,lrData) = scan

    #print '(x,y,theta,deg) = ({0},{1},{2},{3})'.format(x,y,theta,np.rad2deg(theta))

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

    plotPoints(points, c, 'laserICP.png')

    DATA = np.zeros((2,len(points)))
    for i in range(len(points)):
        DATA[0,i] = points[i][0]
        DATA[1,i] = points[i][1]
    return DATA

def singleScanData(scan, offset, rngMax=10.0):
    # plot robot location and heading
    (x,y,theta,_,lrData) = scan

    #print '(x,y,theta,deg) = ({0},{1},{2},{3})'.format(x,y,theta,np.rad2deg(theta))

    lrData = [(brg,rng) for (brg,rng) in lrData if rng < rngMax]
    brgs = [brg for (brg,_) in lrData]
    maxBrg = max(brgs)
    minBrg = min(brgs)
    points = [(x + rng*math.cos(theta+brg) - offset[0], y + rng*math.sin(theta+brg) - offset[1]) for (brg,rng) in lrData]
    DATA = np.zeros((2,len(points)))
    for i in range(len(points)):
        DATA[0,i] = points[i][0]
        DATA[1,i] = points[i][1]
    return x-offset[0],y-offset[1],DATA

def correctScanData(scans, offset):
    # scans[i] = (x,y,theta,ts,readings)
    # corData[i] = (x,y)
    corData = []

    # get first scan
    ax,ay,A = singleScanData(scans[0], offset, 20.0)
    corData.append((ax,ay))

    # ICP settings
    N = 10
    e = 1e-5
    k = 3
    T = np.zeros(2)

    R = np.eye(2)
    T = np.zeros(2)

    for i in range(1,len(scans)):
        # get original offset adjusted position and scan points
        bx,by,B = singleScanData(scans[i], offset, 20.0)
        # run ICP to find R and T
        dR,dT = ICP.laserDataIcp(A,B,N,e,k)

        R = dR.dot(R)
        T = T + dT

        pos = np.array([bx,by])
        pos = (R.dot(pos.T).T + T).T
        corData.append((pos[0], pos[1]))
        A = B

    return corData

def runICP(scans,a,b,fname):
    # ICP settings
    N = 30
    e = 1e-5
    k = 2.5
    maxRng = 20.0
    R,T,A,B,errors,loc = ICP.laserDataIcp(scans[a],scans[b], N, e, k, ICP.pwSSE, maxRng)
    plt.clf()
    plotRobot(0,0,loc[0],'k')
    plotRobot(loc[1],loc[2],loc[3],'b')
    plotPoints(A.T, 'k', fname)
    plotPoints(B.T, 'b', fname)
    C = (R.dot(B).T + T).T
    plotPoints(C.T, 'm', fname)
    return R,T,loc

def plotPos(posData, fname):
    """
    Plot robot position from odometry data
    """
    x = [p[0] for p in posData]
    y = [p[1] for p in posData]
    plt.clf()
    plt.plot(x, y, 'k.')
    plt.axis('equal')
    plt.xlabel('x pos')
    plt.ylabel('y pos')
    plt.savefig('{0}'.format(fname))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print '{0} rawfile corfile'.format(sys.argv[0])
    else:
        frawIn = sys.argv[1]
        fcorIn = sys.argv[2]
        print(frawIn)
        print(fcorIn)

        rawParser = CarmenParser()
        rawParser.parse(frawIn)

        corParser = CarmenParser()
        corParser.parse(fcorIn)

        # plot the corrected position data
        rawParser.plotPos('pos-raw.png')
        corParser.plotPos('pos-cor.png')

        # attempt to plot path data by correcting the raw data using ICP to align consecutive scans
        #correctedRawScanPosData = correctScanData(rawParser.rangeData, rawParser.posOffset)
        #plotPos(correctedRawScanPosData, 'pos-rawcor.png')

        print 'len laser raw {0}'.format(len(rawParser.rangeData))
        plotPos(rawParser.rangeData, 'pos-laser-raw.png')
        print 'len laser cor {0}'.format(len(corParser.rangeData))
        plotPos(corParser.rangeData, 'pos-laser-cor.png')

        """
        thetaLim = np.pi/4
        for i in range(1,len(rawParser.rangeData)):
            if abs(rawParser.rangeData[i][2] - rawParser.rangeData[i-1][2]) > thetaLim:
                print '{0}, {1}'.format(rawParser.rangeData[i][2],rawParser.rangeData[i-1][2])
        """

        outdir = os.path.join(os.getcwd(),'scans')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        pos = [(0,0)]
        for i in range(1,len(rawParser.rangeData)):
            # loc = (atheta, dx, dy, btheta)
            R,T,loc = runICP(rawParser.rangeData,i-1,i,os.path.abspath(os.path.join(outdir,'laserICP-{0}.png'.format(str(i).zfill(6)))))
            dxDy = np.array([loc[1],loc[2]])
            corDxDy = (R.dot(dxDy).T + T).T
            pos.append((pos[-1][0] + corDxDy[0], pos[-1][1] + corDxDy[1]))

        plotPos(pos, 'correctedPos.png');
