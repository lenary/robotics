import numpy as np
import math

"""
Ok, this is a signed distance function.

There will be a single global SDF, which is updated by merging the contents of a local SDF per set of laser readings.

Global SDFs have an offset of 0,0.

These are represented by a numpy 2d array underneath, with each entry being a pair of the average distance to a close surface, and the number of readings that the average was calculated from.

Each entry represents the distance from that point (or the centre of the res by res square) to the surface. To make merging far easier, it's important the points in a local SDF are aligned with the points in the global SDF. This means some slightly ugly munging of distances, etc.
"""
class SDF(object):
    @staticmethod
    def fromRangeLine(logline, res=0.1, threshold=10):
        (x, y, theta, ts, readings) = logline

        size = int(2*threshold/res)

        offset_x = closestResVal(x-threshold, res)
        offset_y = closestResVal(y-threshold, res)

        sdf = SDF(initsize=(size,size), res=res, offset=(offset_x, offset_y))
        sdf.timestamp = ts

        for brg,rng in readings:
            if rng < threshold:
                sdf.addLaserRangeReading(normalizeRadians(theta+brg), rng, origin=(x,y))

        return sdf

    def __init__(self, initsize=(50,50), res=0.1, offset=(0,0)):
        self.size = initsize
        self.res = res
        self.offset = offset

        self._dtype = np.dtype([('dist', np.float32), ('nreadings', np.int)])

        self._array = np.zeros(self.size, dtype=self._dtype)

    """

    """
    def addLaserRangeReading(self, brg, rng, origin=(0,0)):
        visited = set()

        (x,y) = origin

        reading_x = x + (rng * math.cos(brg))
        reading_y = y + (rng * math.sin(brg))

        current = 0.0
        end = 2.0*rng
        step = self.res/5.0
        while current <= end:
            current += step

            step_x = x + (current * math.cos(brg))
            step_y = y + (current * math.sin(brg))

            idx = self._convertCoordsToLocalIdx(step_x, step_y)
            if idx is None:
                continue
            elif idx in visited:
                continue

            dist_from_origin = math.hypot(step_x - x, step_y - y)
            dist = rng - dist_from_origin # roughly correct

            ave_dist, count = self._array[idx]
            new_ave = (ave_dist * count + dist)/(count+1)
            self._array[idx] = (new_ave, count+1)

    def _convertCoordsToLocalIdx(self, x, y):
        x = int(round((x-self.offset[0])/self.res))
        y = int(round((y-self.offset[1])/self.res))
        if (0 <= x < self.size[0]) and (0 <= y < self.size[1]):
            return (x,y)
        else:
            return None

    """
    This only works if the resolution is the same, and the offsets align based on the resolution.

    This should expand the current SDF to include all the local data, averaged into the global view.
    """
    def addLocalSDF(self, local_sdf):
        pass

    def plot(self, plt):
        plt.clf()
        plt.axis('equal')
        plt.imshow(self._array['dist'])
        # plt.imshow(self._array['nreadings'])
        plt.show()

def closestResVal(value, res):
    return res*round(value/res)

def normalizeRadians(radians):
    x = radians
    while x >= 2*math.pi:
        x -= 2*math.pi
    while x < 0:
        x += 2*math.pi
    return x
