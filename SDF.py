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
    def fromRangeLine(logline, res=0.1, threshold=10, truncate=1.0):
        (x, y, theta, ts, readings) = logline

        size = int(2*threshold/res)

        offset_x = closestResVal(x-threshold, res)
        offset_y = closestResVal(y-threshold, res)

        sdf = SDF(initsize=(size,size), res=res, offset=(offset_x, offset_y), truncate=truncate)
        sdf.timestamp = ts

        for brg,rng in readings:
            if rng < threshold:
                sdf.addLaserRangeReading(normalizeRadians(theta+brg), rng, origin=(x,y))

        return sdf

    def __init__(self, initsize=(50,50), res=0.1, offset=(0,0), truncate=False):
        self.size = initsize
        self.res = res
        self.offset = offset
        self.truncate = truncate

        self._dtype = np.dtype([('dist', np.float32), ('nreadings', np.int)])

        self._array = np.zeros(self.size, dtype=self._dtype)

    def is_global_sdf(self):
        return self.offset == (0,0) and self.truncate is False

    """

    """
    def addLaserRangeReading(self, brg, rng, beam_origin=(0,0)):
        visited = set()

        (x,y) = beam_origin

        # Now we begin marching out from beam_origin
        current = 0.0
        step = self.res/5.0
        while True:
            # Take a step
            current += step

            # Where are we in the global space?
            step_x = x + (current * math.cos(brg))
            step_y = y + (current * math.sin(brg))

            # Where is the closest square to where we are right now
            idx = self._convertCoordsToLocalIdx(step_x, step_y)
            if idx is None:
                # we're outside the limits of the course, stop
                # marching

                # this is the only thing that breaks from the while
                # loop if we're not truncating
                break
            elif idx in visited:
                # we've already visited this square, march onwards
                continue
            else:
                # Find out where the centre of this square is, and how
                # far it is from the origin
                sq_x, sq_y = self._convertLocalIdxToCoords(idx)
                sq_to_origin = math.hypot(sq_x - x, sq_y - y)

                # The distance to the surface is roughly the
                # difference between the range reading and how far the
                # current square is from the origin.
                #
                # This only works because we know the box is maximally
                # about self.res away from the line we're walking along
                dist = rng - sq_to_origin

                # We've done this square, we don't want to update it (e.g)
                # 3 times if 3 steps all land in the same square
                visited.add(idx)

                # Add the reading to the local SDF
                res = self._addSurfaceDistanceReading(idx, dist)

                # This is used in the TSDF to stop adding readings
                # once the distance is too negative, i.e. we've gone
                # far beyond the range reading.
                if res == "STOP":
                    break

    def _addSurfaceDistanceReading(self, idx, dist, count=1):
        # This just accumulates the averages
        ave_dist, old_count = self._array[idx]
        new_count = old_count + count
        new_ave_dist = (ave_dist * old_count + dist * count) / new_count

        if (not self.truncate) or (abs(dist) < self.truncate):
            # If we're not truncating, or the distance is less than
            # the truncated distance from the surface, add the reading
            self._array[idx] = (new_ave_dist, new_count)
        elif dist < -self.truncate:
            # Otherwise, we've gone further than the truncate distance
            # beyond the surface reading, so stop
            return "STOP"


    def _convertCoordsToLocalIdx(self, x, y):
        x = int(round((x-self.offset[0])/self.res))
        y = int(round((y-self.offset[1])/self.res))

        if (0 <= x < self.size[0]) and (0 <= y < self.size[1]):
            return (x,y)
        else:
            return None

    def _convertLocalIdxToCoords(self, idx):
        x = (idx[0] * self.res) + self.offset[0]
        y = (idx[1] * self.res) + self.offset[1]
        return (x,y)

    """
    This only works if the resolution is the same, and the offsets align based on the resolution.

    This should expand the current SDF to include all the local data, averaged into the global view.
    """
    def addLocalSDF(self, local_sdf):
        # TODO:
        # - Check resolutions agree
        # - Check sizes, expand local array if needed
        # - Iterate over local_sdf, adding all entries from it to the
        # current sdf. Use coordinate conversion methods above, and
        # _addSurfaceDistanceReading with all 3 arguments.
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
