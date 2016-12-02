import numpy as np
import math
import os

# check for DISPLAY environment variable and use matplotlib 'Agg' if not present
# this is a check for Windows Linux Subsystem GUI issues
if os.environ.get('DISPLAY') is None:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt


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
    def correspondingPoints(a,b,distfunc):
        """
        for each point in a, find the point in b that is closest, according
        to the provided distance function 'distfunc'

        assumes a and b have shapes (k,m) and (k,n), where k is the dimensionality
        of the points, and m and n are the number of points in a and b, respectively

        returns tuple (P,E) where P is shape (k,m) and E is shape (m)
        P is the set of corresponding points, E is the distance for each point

        """
        P = np.zeros(a.shape)
        E = np.zeros(a.shape[1])
        i = 0
        for i in range(a.shape[1]):
            e, ip = distfunc(a[:,i], b)
            E[i] = e
            P[:,i] = b[:,ip]

        # TODO: trim points that have error > median(E) (or too large error)
        # TODO: if two points in a match to same point in b, only use one of the points
        #       from a, choosing the one with smaller distance

        return P,E

    @staticmethod
    def computeRT(A,B,T):
        """
        Given corresponding point sets A and B, compute the Rotation and Translation
        matrices R and T, respectively, using SVD, with initial translation guess T

        A are the model points
        B are the data points
        """
        R = ((A.T - T).T).dot(B.T).dot(np.linalg.inv(B.dot(B.T)))
        U,S,VT = np.linalg.svd(R)
        R = U.dot(VT)
        T = np.mean(A,1) - np.mean(R.dot(B),1)
        return R, T

    @staticmethod
    def removeOutliers(P):
        """
        Remove outlier points from dataset P
        """
        # TODO: for now , just return P (i.e., assume no outliers)
        return P

    @staticmethod
    def icp(A, B, N=5, e=1e-5):
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
        # list of the errors at each step
        errors = list()

        # Rotation and Translation matrices
        R = np.eye(A.shape[0])
        T = np.zeros(A.shape[0])
        # current error
        old_error = np.inf

        # create a copy of Bp that is used and updated each iteration
        # Bi should get closer to Ap every iteration
        Bi = np.copy(Bp)

        # compute a guess at T based on difference of centers of mass
        # under the assumption of a small rotational angle theta, this should be okay
        Ti = muB - muA

        # iterate until error converges, or up to N iterations
        for i in range(N):
            # find the corresponding points in Ap for every point in Bi, using pwSSE as distance
            # Bi is the current set of data points (with incremental R and T applied each step)
            # P is the set of points from Ap that are closest to each point in Bi
            P,E = ICP.correspondingPoints(Bi,Ap,ICP.pwSSE)

            # compute the error
            err = np.sum(E)
            print err
            # check for convergence
            if abs(old_error - err) > e:
                old_error = err
            else:
                break

            # compute dR and dT, the incremental Rotation and Translation matrices between
            # corresponding point sets Bi and P
            dR,dT = ICP.computeRT(P,Bi,Ti)

            # update R using dR
            R = dR.dot(R)

            # update guess for T
            T = np.mean(A,1) - np.mean(R.dot(B),1)

            # update incremental T
            Ti = np.mean(Ap,1) - np.mean(R.dot(Bi),1)
            # update Bi
            Bi = (R.dot(Bp).T + Ti).T

        return R,T
        

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
    
    # some random points
    A = np.random.rand(2,10)
    theta = (np.random.rand(1)*np.pi/4)[0]
    T = np.array([2.5, 5.0])
    # build rotation matrix
    R = np.array([[np.cos(theta), -1.0*np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    # compute rotated and translated points (perfect transformation)
    B = (R.dot(A).T + T).T

    print np.mean(A,1)
    print np.mean(B,1)

    r,t = ICP.icp(A,B,N,e)
    print '\nActual R:'
    print R
    print 'Estimated R:'
    print r

    print '\nActual T:'
    print T
    print 'Estimated T:'
    print t
    
    # B = (R.dot(A).T + T).T
    # A' = rotation and translation from B that should transform the points to A
    Ap = (r.dot(B).T + t).T
    plotSets(A, B, Ap, r, T)
