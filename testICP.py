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

def runICP(scans,a,b,fname):
    # ICP settings
    N = 30
    e = 1e-5
    k = 2.5
    maxRng = 20.0
    # loc = (atheta, dx, dy, btheta)
    R,T,A,B,errors,loc = ICP.laserDataIcp(scans[a],scans[b], N, e, k, ICP.pwSSE, maxRng)
    theta = ICP.recoverTheta(R)
    plt.clf()
    plt.title('(Rtheta, Tdx, Tdy) = ({0:.5f}, {1:.2f}, {2:.2f}), (dx,dy,dtheta) = ({3:.2f},{4:.2f},{5:.2f})'.format(theta,T[0],T[1],loc[1],loc[2],ICP.normalizeAngle(loc[3]-loc[0])))
    plotRobot(0,0,loc[0],'k')
    plotRobot(loc[1],loc[2],loc[3],'b')
    plotPoints(A.T, 'k', fname)
    plotPoints(B.T, 'b', fname)
    C = (R.dot(B).T + T).T
    plotPoints(C.T, 'm', fname)
    return R,T,loc,A,B,C

def plotPos(posData, fname, c='r'):
    """
    Plot robot position from odometry data
    """
    x = [p[0] for p in posData]
    y = [p[1] for p in posData]
    plt.clf()
    plt.plot(x, y, '{0}.'.format(c))
    plt.axis('equal')
    plt.xlabel('x pos')
    plt.ylabel('y pos')
    plt.savefig('{0}'.format(fname))

#def plotRotatedPos

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print '{0} rawfile corfile'.format(sys.argv[0])
    else:
        frawIn = sys.argv[1]
        fcorIn = sys.argv[2]
        print(frawIn)
        print(fcorIn)

        # parse the raw laser data
        rawParser = CarmenParser()
        rawParser.parse(frawIn)

        # parse the corrected laser data
        corParser = CarmenParser()
        corParser.parse(fcorIn)

        # plot the corrected position data
        rawParser.plotPos('pos-raw.png')
        corParser.plotPos('pos-cor.png')

        print 'len laser raw {0}'.format(len(rawParser.rangeData))
        plotPos(rawParser.rangeData, 'pos-laser-raw.png')
        print 'len laser cor {0}'.format(len(corParser.rangeData))
        plotPos(corParser.rangeData, 'pos-laser-cor.png')

        outdir = os.path.join(os.getcwd(),'scans')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        pos = [(0,0)]
        totalTheta = 0.0
        for i in range(1,200):#len(rawParser.rangeData)):
            # loc = (atheta, dx, dy, btheta)
            R,T,loc,A,B,C = runICP(rawParser.rangeData,i-1,i,os.path.abspath(os.path.join(outdir,'laserICP-{0}.png'.format(str(i).zfill(6)))))

            # at each step, plot A rotated by totalTheta, which is the accumulated rotational offset
            Rp = np.array([[np.cos(totalTheta), -1.0*np.sin(totalTheta)],[np.sin(totalTheta), np.cos(totalTheta)]])
            Ap = Rp.dot(A)
            plotPos(Ap.T,os.path.abspath(os.path.join(outdir,'laserICP-{0}-cor.png'.format(str(i).zfill(6)))))

            # rotation angle from ICP (i.e., correction to rotation)
            rTheta = ICP.recoverTheta(R)
            # rotation reported by robot between consecutive scan poses
            dTheta = ICP.normalizeAngle(loc[3]-loc[0])
            # translation reported by robot between consecutie scan poses
            dx = loc[1]
            dy = loc[2]

            totalTheta = totalTheta + rTheta
            print '{0}: {1:.5f}, {2:.2f}, {3:.5f}'.format(i, totalTheta, np.rad2deg(totalTheta), rTheta)

            dxDy = np.array([loc[1],loc[2]])
            Rp = np.array([[np.cos(totalTheta), -1.0*np.sin(totalTheta)],[np.sin(totalTheta), np.cos(totalTheta)]])
            corDxDy = (Rp.dot(dxDy).T + T).T
            pos.append((pos[-1][0] + corDxDy[0], pos[-1][1] + corDxDy[1]))




        plotPos(pos, 'correctedPos.png');
