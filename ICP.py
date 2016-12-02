import numpy as np
import math

"""
Iterative Closest Point algorithm
"""
class ICP(object):
    def __init__(self):
        pass

    """

    """

    """
    Pieces of the implementation borrowed from:
    http://connor-johnson.com/2014/06/06/an-iterative-closest-point-algorithm/
    """
    @staticmethod
    def SSE(a, b):
        return np.sum(np.array(a-b)**2.0)

    @staticmethod
    def pwSSE(p, a):
        """
        Point-wise smallest squared error.
        Distance from p to closest point in a
        """
        difference = p - a
        # x and y columns
        xcol = np.ravel(difference[:,0])
        ycol = np.ravel(difference[:,1])
        # sum of the squared differences between p
        sqr_difference = xcol**2.0 + ycol**2.0
        # nearest squared distance
        distance = np.min(sqr_difference)
        # index of nearest point to p in a
        nearest_pt = np.argmin(sqr_difference)
        return distance, nearest_pt

    @staticmethod
    def NSSE(a, b):
        """
        Nearest point sum squared error
        """
        err = 0.0
        for p in a:
            d, np = ICP.pwSSE(p, b)
            err += d
        return err

    @staticmethod
    def icp(a, b, iterations=1, e=1e-5):
        """
        Returns the translation and rotation matrices that transform points a
        to points b, with an error threshold of e

        The goal is to minimize the NSSE between a and b until it is less than e
        or the full number of iterations has occurred

        a: n by 2 numpy.ndarray of (x,y) point locations
        b: n by 2 numpy.ndarray of (x,y) point locations

        Assumptions:
        a and b have same shape, namely (n, 2)
        """

        # list of the errors at each step
        errors = list()
        errors.append(ICP.NSSE(a,b))
        print errors[-1]

        if (a.shape != b.shape):
            return None

        n = a.shape[0]

        # compute center of mass for a and b
        mu_a = np.sum(a,0) / float(n)
        mu_b = np.sum(b,0) / float(n)

        # subtract center of mass from a and b
        ap = a - mu_a
        bp = b - mu_b

        # SVD
        # compute W
        W = np.zeros((2,2))
        for i in range(ap.shape[0]):
            W += np.outer(ap[i,:].T, bp[i,:])
        U,S,V = np.linalg.svd(W)

        # compute Rotation and translation
        # If perfect, a = Rb + t
        #R = U.dot(V.T)
        R = V.dot(U.T) # this is just a transpose of the above
        t = mu_b.dot(R) - mu_a # provides translation from a to b

        return R,t


if __name__ == '__main__':

    # create some test points
    A = np.array([[1.0,0.0],[2.0,3.5],[2.2,3.6],[6.6,2.1],[3.0,4.3]])
    # define a rotation theta and translation T
    theta = np.pi/24
    T = np.array([0.1, 0.1])
    # build rotation matrix
    R = np.array([[np.cos(theta), -1.0*np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    # compute rotated and translated points
    B = A.dot(R.T) + T
    print A
    print B

    r,t = ICP.icp(A,B)
    print '\nActual R:'
    print R
    print 'Estimated R:'
    print r

    print '\nActual T:'
    print T
    print 'Estimated T:'
    print t
