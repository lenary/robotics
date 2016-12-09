import numpy as np
import math
import os

# check for DISPLAY environment variable and use matplotlib 'Agg' if not present
# this is a check for Windows Linux Subsystem GUI issues
if os.environ.get('DISPLAY') is None:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
from operator import itemgetter

"""
Iterative Closest Point algorithm
"""
class ICP(object):
    def __init__(self):
        pass
    """
    Sources for algorithm and implementation ideas:
    * http://connor-johnson.com/2014/06/06/an-iterative-closest-point-algorithm/
    * https://engineering.purdue.edu/kak/distICP/ICP-2.0.html
    """
    @staticmethod
    def pwSSE(p, a):
        """
        Point-wise smallest squared error.
        Squared Distance from p to closest point in a
        p is a single point [x,y]
        a is an 2 x n array of points (column vector per point)

        This function serves as the default error/distance metric for ICP
        """
        difference = (a.T - p).T
        # x and y columns
        xcol = np.ravel(difference[0,:])
        ycol = np.ravel(difference[1,:])
        # sum of the squared differences between p
        sqr_difference = xcol**2.0 + ycol**2.0
        # nearest squared distance
        distance = np.min(sqr_difference)
        # index of nearest point to p in a
        nearest_pt = np.argmin(sqr_difference)
        return distance, nearest_pt

    @staticmethod
    def correspondingPoints(a,b,sf=3.0,distfunc=pwSSE):
        """
        For each point in a, find the point in b that is closest, according
        to the provided distance function 'distfunc'.

        Assumes a and b have shapes (k,m) and (k,n), where k is the dimensionality
        of the points, and m and n are the number of points in a and b, respectively

        Returns tuple (A, B, E, iA, iB) where A and B are shape (k,w), and E is shape(w)
        A is the set of selected points from a
        B is the set of corresponding points from b
        E[i] is the distance/error value from A[i] to B[i]
        iA are the indices of points selected from a
        iB are the indices of corresponding points from b

        Pairs A[i], B[i] are only accepted if they are not outliers if
        E[i] <= T*stdv(E)

        sf is a scale factor for outlier thresholding. 

        """
        # list of (error, index of element in b)
        PE = []
        i = 0
        for i in range(a.shape[1]):
            # compute error from a[:,i] to nearest point in b
            # also return ib, the index in b of nearest neighbor to a[:,i]
            e, ib = distfunc(a[:,i], b)
            PE.append((e,ib))

        # compute outlier rejection threshold
        stdv = np.std([pe[0] for pe in PE])
        RT = sf*stdv

        bToa = {}
        for i in range(len(PE)):
            if (PE[i][0] <= RT):
                data = (a[:,i],b[:,PE[i][1]],PE[i][0],i,PE[i][1])
                if not PE[i][1] in bToa:
                    bToa[PE[i][1]] = [data]
                else:
                    bToa[PE[i][1]].append(data)

        # at this point, each point in a is matched to a point in b
        # but, the same point in b may be present in B multiple times
        # (i.e., it may be the closest point to multiple points from a)
        data = []
        for k,v in bToa.iteritems():
            data.append(min(v,key=itemgetter(2)))

        
        A,B,E,iA,iB = zip(*data)
        return np.array(A).T, np.array(B).T, np.array(E), iA, iB

    @staticmethod
    def computeRotation(A,B):
        """
        Given corresponding point sets A and B, compute the Rotation and Translation
        matrices R and T, respectively, using SVD, with initial translation guess T

        A are the model points
        B are the data points
        """
        R = A.dot(B.T).dot(np.linalg.inv(B.dot(B.T)))
        U,S,VT = np.linalg.svd(R)
        R = U.dot(VT)
        return R

    @staticmethod
    def computeOriginRotation(A,B):
        """
        Compute the rotation matrix R that aligns points in A to points in B,
        where A and B are first centered on the origin.

        A and B are corresponding point sets, but are not required to be centered
        around the origin, and will be recentered to the origin before computing R

        A are the model points
        B are the data points
        """
        # center the point sets on the origin
        muA = np.mean(A,1)
        muB = np.mean(B,1)
        Ap = (A.T - muA).T
        Bp = (B.T - muB).T        

        R = Ap.dot(Bp.T).dot(np.linalg.inv(Bp.dot(Bp.T)))
        U,S,VT = np.linalg.svd(R)
        R = U.dot(VT)
        return R

    @staticmethod
    def icp(A, B, N=5, e=1e-5, k=3, distfunc=pwSSE):
        """
        Returns the Rotation and Translation matrices, R and T, that transform points in B
        to points in A.
        
        Algorithm runs until error between iterations is <= e, or N iterations occur.

        Assumptions:
        1. A and B have shape (k,n) and (k,m), respectively, where k is the dimensionality of the
        points in A and B
        2. A are the model points, B are the data points
        3. The Rotation and Translation between A and B are small, especially Rotation, for some
        definition of small
        """

        assert(A.shape[0] == B.shape[0])

        ### ICP Pre-Processing
        # clean up the point sets
        A = ICP.removeOutliers(A)
        B = ICP.removeOutliers(B)

        # create point sets Ap and Bp that are centered around (0,0) by subtracting
        # centers of mass muA and muB from A and B, respectively
        muA = np.mean(A,1)
        muB = np.mean(B,1)
        Ap = (A.T - muA).T
        Bp = (B.T - muB).T

        ### Initialize variables for iteration
        # rotation and translation matrices
        R=np.eye(A.shape[0])
        T=np.zeros(A.shape[0])
        # list of the errors at each step
        errors = list()
        # current error
        old_error = np.inf

        ### ICP iteration

        # create a copy of Bp that is used and updated each iteration
        # Bi should get closer to Ap every iteration
        Bi = np.copy(Bp)

        # compute a guess at T based on difference of centers of mass
        # under the assumption of a small rotational angle theta, this should be okay
        Ti = np.copy(T)

        # iterate until error converges, or up to N iterations
        for i in range(N):
            # find the corresponding points between A and Bp
            cBi, cAp, E, iB, iA = ICP.correspondingPoints(Bi,Ap,k,distfunc)

            # compute the error
            err = np.sum(E)
            #print err
            # check for convergence
            if abs(old_error - err) > e:
                old_error = err
            else:
                break

            # compute dR, the incremental Rotation matrix between
            # corresponding point sets pBi and pAp
            dR = ICP.computeRotation(cAp,cBi)

            # update R using dR
            R = dR.dot(R)
            # update guess for T using the corresponding point sets
            AA = np.array([A[:,i] for i in iA]).T
            BB = np.array([B[:,i] for i in iB]).T
            T = np.mean(AA,1) - np.mean(R.dot(BB),1)

            Ti = np.mean(Ap,1) - np.mean(R.dot(Bi),1)
            Bi = (R.dot(Bp).T + Ti).T

        return R,T
        

    @staticmethod
    def convertScanToPoints(scan, maxRng=10.0, robotCentric=False):
        """
        Convert scan data, produced by the CarmenParser, to a set of points (x,y)

        Points greater than maxRange in distance are removed

        Setting robotCentric true gives points in (x,y) coordinates with the robot at the origin
        Setting to false provides the points in the global coordinate system as given by the
        (x,y) position of the robot from the scan data

        Returns the a list of points (x,y) and the (x,y,theta) position of the robot from scan data
        """
        # extract the scan fields
        (x,y,theta,_,lrData) = scan

        # remove points further than maxRng from robot
        lrData = [(brg,rng) for (brg,rng) in lrData if rng < maxRng]

        # compute (x,y) points relative to robot
        points = [(rng*math.cos(theta+brg), rng*math.sin(theta+brg)) for (brg,rng) in lrData]

        # adjust points to global coordinates if requested
        if not robotCentric:
            points = [(x + p[0], y + p[1]) for p in points]

        # convert to a numpy array
        points = np.array([list(zip(*points)[0]),list(zip(*points)[1])])

        return x,y,theta,points


    @staticmethod
    def normalizeAngle(theta):
        """
        Given an angle in rads, convert angle to [0,2*pi]
        """
        while theta < -np.pi:
            theta = theta + 2*np.pi
        while theta > 2*np.pi:
            theta = theta - 2*np.pi
        return theta

    @staticmethod
    def laserDataIcp(scanA, scanB, N=5, e=1e-5, k=3, distfunc=pwSSE, maxRng=10.0, estTheta=0, estT=(0,0)):
        """
        ICP specialized for robot laser range scanner data

        Returns the Rotation and Translation matrices, R and T, that transform points in B
        to points in A.
        
        estTheta is a guess at initial rotation offset between A and B
        estT is a guess at initial translation offset between A and B

        Algorithm runs until error between iterations is <= e, or N iterations occur.
        Assumptions:
        1. A and B have shape (k,n) and (k,m), respectively, where k is the dimensionality of the
        points in A and B
        2. A are the model points, B are the data points
        3. The Rotation and Translation between A and B are small, especially Rotation, for some
        definition of small
        """

        ### ICP Pre-Processing
        # get point sets from scans
        # points returned assuming the robot is located at the origin
        ax,ay,atheta,A = ICP.convertScanToPoints(scanA, maxRng, True)
        atheta = ICP.normalizeAngle(atheta)
        # points from B are located with the robot at its perceived location
        bx,by,btheta,B = ICP.convertScanToPoints(scanB, maxRng, False)
        btheta = ICP.normalizeAngle(btheta)
        # adjust points in B by the robot location of scanA
        B = (B.T - np.array([ax,ay])).T

        # compute delta position between robot from scanA to scanB
        # also normalizes dtheta to [0,2*pi]
        dx,dy,dtheta = (bx-ax),(by-ay),ICP.normalizeAngle(btheta-atheta)

        ## !! If there is no noise, the two scans will be aligned (as they were in the global frame
        ## of reference), but with the robot at scanA located at the origin, and the robot at
        ## scanB located at (dx,dy).
        ## If there is noise the two point sets will be off by some R and T
        ## running ICP should find these matrices.

        ### Initialize variables for iteration
        # rotation and translation matrices
        # ensure same dimensionality of points in A and B
        # note: A and B are not required to contain the same number of points
        assert(A.shape[0] == B.shape[0])
        # initial guess of R and T if provided
        R = np.array([[np.cos(estTheta), -1.0*np.sin(estTheta)],[np.sin(estTheta), np.cos(estTheta)]])
        T = np.array([estT[0], estT[1]])
        
        # list of the errors at each step
        errors = list()
        # current error
        old_error = np.inf

        # current B points at iteration i
        Bi = np.copy(B)
        # cumulative T matrix at iteration i
        Ti = np.copy(T)

        ### ICP iteration

        # iterate until error converges, or up to N iterations
        for i in range(N):
            # First, apply current R and T to B, then find the corresponding points between
            # A and B, to generate a new point set B'
            Bp = (R.dot(B).T + T).T

            # find the corresponding points between A and Bp
            cBp, cA, E, iB, iA = ICP.correspondingPoints(Bp,A,k,ICP.pwSSE)

            # compute the error
            err = np.sum(E)
            #print err
            # check for convergence
            if abs(old_error - err) > e:
                old_error = err
            else:
                break

            # compute dR, the incremental Rotation matrix between
            # corresponding point sets pBi and pAp
            dR = ICP.computeRotation(cA,cBp)

            # update R using dR
            R = dR.dot(R)

            ### TODO: how do we compute T?
            # the means of the two point sets will not align at the end because the
            # two point sets have no perfect correspondence

            # update guess for T using the corresponding point sets
            T = np.mean(cA,1) - np.mean(R.dot(cBp),1)

        return R,T,A,B,errors,(atheta,dx,dy,btheta)

    @staticmethod
    def recoverRT(r,t):
        """
        Given R and T matrices that transform points from B to A,
        return R' and T' matrices that transfrom points from A to B
        """
        return r.T, -t.dot(np.linalg.inv(r.T))


def plotSets(A,B,C,R,T):
    plt.clf()
    plt.plot([x[0] for x in A.T], [x[1] for x in A.T], 'k.')
    plt.plot([x[0] for x in B.T], [x[1] for x in B.T], 'b.')
    plt.plot([x[0] for x in C.T], [x[1] for x in C.T], 'r+')

    muA = np.mean(A,1)
    muB = np.mean(B,1)
    muC = np.mean(C,1)

    plt.plot(muA[0],muA[1],'ko')
    plt.plot(muB[0],muB[1],'bo')
    plt.plot(muC[0],muC[1],'ro')

    plt.axis('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('points.png')

if __name__ == '__main__':

    N = 10
    e = 1e-5
    k = 3

    # some random points
    A = np.random.rand(2,10)
    theta = (np.random.rand(1)*np.pi/4)[0]
    T = np.random.rand(1,2)#np.array([2.5, 5.0])
    # build rotation matrix
    R = np.array([[np.cos(theta), -1.0*np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    # compute rotated and translated points (perfect transformation)
    B = (R.dot(A).T + T).T

    print np.mean(A,1)
    print np.mean(B,1)

    r,t = ICP.icp(A,B,N,e,k)
    print '\nActual R:'
    print R
    print 'Estimated R:'
    print r

    print '\nActual T:'
    print T
    print 'Estimated T:'
    print t

    Rp,Tp = ICP.recoverRT(r,t)
    print '\nRecovered R:'
    print Rp
    print 'Recovered T:'
    print Tp

    # B = (R.dot(A).T + T).T
    # A' = rotation and translation from B that should transform the points to A
    Ap = (r.dot(B).T + t).T
    plotSets(A, B, Ap, r, T)

    if (np.sum(T-Tp) < e and np.sum(R-Rp) < e):
        print 'pass'
    else:
        print 'fail'
