#!/usr/bin/env python

import sys
import os
import math

import numpy as np

from CarmenParser import CarmenParser
from SDF import SDF
from ICP import ICP

# check for DISPLAY environment variable and use matplotlib 'Agg' if not present
# this is a check for Windows Linux Subsystem GUI issues
if os.environ.get('DISPLAY') is None:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

def runICP(scana, scanb):
    """
    Run ICP on two scans, given by indices a and b, saving a plot
    of the two scans, plus the "corrected" scan to fname

    Rreturns the rotation and translation matrices, R and T;
    the location vector giving heading of scan a, (dx,dy) to robot in
    scan b, and heading of scan b; and the point sets A, B, and C, where
    A has the robot at the origin, robot for B is at (dx,dy), and C is
    the correct scan b
    """

    # ICP settings
    N = 30
    e = 1e-5
    k = 2.5
    maxRng = 10.0
    # loc = (atheta, dx, dy, btheta)
    R,T,A,B,errors,loc = ICP.laserDataIcp(scana, scanb, N, e, k, ICP.pwSSE, maxRng)
    theta = ICP.recoverTheta(R)
    # compute "corrected" scanB
    C = (R.dot(B).T + T).T

    return R,T,loc,A,B,C


def plotScan(scans, name):
    """
    Plot laser range data scan taken with robot at specified position.

    This plots position data stored with the laser scan data

    Commented out lines are for plotting the laser scan data points
    """

    res = 0.10

    global_sdf = SDF.asGlobal(initextents=(200.0, 200.0), res=res)

    # Run ICP for each pair of consecutive frames
    pos = [(0,0)]
    totalTheta = 0.0

    prevScan = None

    sys.stdout.write("Scanning %s..." % name)
    sys.stdout.flush()

    c = 0
    for scan in scans:
        if prevScan:
            R,T,loc,A,B,C = runICP(scan, prevScan)

            # at each step, plot A rotated by totalTheta, which is the accumulated rotational offset
            Rp = np.array([[np.cos(totalTheta), -1.0*np.sin(totalTheta)],[np.sin(totalTheta), np.cos(totalTheta)]])
            Ap = Rp.dot(A)

            # rotation angle from ICP (i.e., correction to rotation)
            rTheta = ICP.recoverTheta(R)

            # accumulate rotation
            totalTheta = totalTheta + rTheta

            # compute "corrected" position
            dxDy = np.array([loc[1],loc[2]])
            Rp = np.array([[np.cos(totalTheta), -1.0*np.sin(totalTheta)],[np.sin(totalTheta), np.cos(totalTheta)]])
            corDxDy = (Rp.dot(dxDy).T + T).T

            currPos = pos[-1]

            newLaserData = []
            for p in Ap.T:
                x,y = p
                newBrg = np.arctan2(y,x) - loc[0]
                newRng = np.hypot(y,x)
                newLaserData.append((newBrg,newRng))

            corrScan = (currPos[0], currPos[1], loc[0], None, newLaserData)

            lastPos = (pos[-1][0] + corDxDy[0], pos[-1][1] + corDxDy[1])
            pos.append(lastPos)


            local_sdf = SDF.fromRangeLine(corrScan, res=res, threshold=10)
            global_sdf.addLocalSDF(local_sdf)

            sys.stdout.write(".")
            sys.stdout.flush()


        prevScan = scan

        c += 1
        if (c % 10) == 0:
            sys.stdout.write(" %d " % c)
        if (c > 50):
            break

    sys.stdout.write("\n")
    print "Scanning Complete: %s" % name

    global_sdf.plotSurfaces(name, prefix="corr-")
    global_sdf.plot(name, prefix="corr-")

if __name__ == '__main__':
    if len(sys.argv) == 2:
        fpin = sys.argv[1]

        parser = CarmenParser()
        parser.parse(fpin)

        plotScan(parser.rangeData, fpin.rsplit("/",1)[0])
