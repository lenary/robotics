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
    def correspondingPoints(a,b,T,distfunc):
        """
        For each point in a, find the point in b that is closest, according
        to the provided distance function 'distfunc'.

        Assumes a and b have shapes (k,m) and (k,n), where k is the dimensionality
        of the points, and m and n are the number of points in a and b, respectively

        Returns tuple (A, B, E) where A and B are shape (k,w), and E is shape(w)
        A is the set of selected points from a
        B is the set of corresponding points from b
        E[i] is the distance/error value from A[i] to B[i]

        Pairs A[i], B[i] are only accepted if they are not outliers if
        E[i] <= T*stdv(E)

        T is a scale factor for outlier thresholding. 

        """

        # list of (error, index of element in a, index of element in b)
        PE = []
        i = 0
        for i in range(a.shape[1]):
            # compute error from a[:,i] to nearest point in b
            # also return ib, the index in b of nearest neighbor to a[:,i]
            e, ib = distfunc(a[:,i], b)
            PE.append((e,ib))

        # compute outlier rejection threshold
        stdv = np.std([pe[0] for pe in PE])
        RT = T*stdv

        A = []
        B = []
        E = []
        for i in range(len(PE)):
            if (PE[i][0] <= RT):
                A.append(a[:,i])
                B.append(b[:,PE[i][1]])
                E.append(PE[i][0])

        return np.array(A).T, np.array(B).T, np.array(E)

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
    def computeOriginRotation(A,B):
        """
        Compute the rotation matrix R that aligns points in A to points in B,
        where A and B are first centered on the origin.

        A and B are corresponding point sets, but are not required to be centered
        around the origin

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
    def icp(A, B, N=5, e=1e-5, k=3):
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

        #cBp, cAp, cE = ICP.correspondingPoints(Bp,Ap,k,ICP.pwSSE)
        #Rorig = ICP.computeOriginRotation(cAp,cBp)
        #print Rorig

        # iterate until error converges, or up to N iterations
        for i in range(N):
            # find the corresponding points in Ap for every point in Bi, using pwSSE as distance
            # Bi is the current set of data points (with incremental R and T applied each step)
            # P is the set of points from Ap that are closest to each point in Bi
            pBi, pAp, E = ICP.correspondingPoints(Bi,Ap,k,ICP.pwSSE)

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
            dR,dT = ICP.computeRT(pAp,pBi,Ti)

            # update R using dR
            R = dR.dot(R)

            # update guess for T
            T = np.mean(A,1) - np.mean(R.dot(B),1)

            # update incremental T
            Ti = np.mean(Ap,1) - np.mean(R.dot(Bi),1)
            # update Bi
            Bi = (R.dot(Bp).T + Ti).T

        return R,T
        
    @staticmethod
    def convertRT(r,t):
        """
        Given R and T matrices that transform points from A to B,
        return R' and T' matrices that transfrom points from B to A
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
    T = np.array([2.5, 5.0])
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

    # this gives t also
    # r is a rotation performed on points B
    # R is a rotation performed on points A
    # this is why t and T are different; applying the rotations to the two point sets
    # aligns the orientation of the set, but results in different translations to bring
    # the centers of mass together
    #print (np.mean(A,1) - np.mean(r.dot(B),1))
    
    # this gives us the original T translation vector
    #print -t.dot(np.linalg.inv(r.T))

    Rp,Tp = ICP.convertRT(r,t)
    print '\nRecovered R:'
    print Rp
    print 'Recovered T:'
    print Tp

    # B = (R.dot(A).T + T).T
    # A' = rotation and translation from B that should transform the points to A
    Ap = (r.dot(B).T + t).T
    plotSets(A, B, Ap, r, T)
