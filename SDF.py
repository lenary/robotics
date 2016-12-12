import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.transforms as plt_tr

"""
Ok, this is a signed distance function.

There will be a single global SDF, which is updated by merging the contents of a local SDF per set of laser readings.

Global SDFs have an offset of 0,0.

These are represented by a numpy 2d array underneath, with each entry being a pair of the average distance to a close surface, and the number of readings that the average was calculated from.

Each entry represents the distance from that point (or the centre of the res by res square) to the surface. To make merging far easier, it's important the points in a local SDF are aligned with the points in the global SDF. This means some slightly ugly munging of distances, etc.
"""
class SDF(object):
    @staticmethod
    def asGlobal(initextents=(100.0, 100.0), res=0.1):
        size_y = int(initextents[0]/res)
        size_x = int(initextents[1]/res)

        offset_y = closestResVal(-initextents[0]/2, res)
        offset_x = closestResVal(-initextents[1]/2, res)

        return SDF(initsize=(size_y,size_x), res=res, offset=(offset_y,offset_x,), truncate=False)

    @staticmethod
    def fromRangeLine(logline, res=0.1, threshold=10.0, truncate=1.0):
        (x, y, theta, ts, readings) = logline

        size = int(2*threshold/res)

        offset_y = closestResVal(y-threshold, res)
        offset_x = closestResVal(x-threshold, res)

        sdf = SDF(initsize=(size,size), res=res, offset=(offset_y, offset_x), truncate=truncate)
        sdf.timestamp = ts

        for brg,rng in readings:
            if rng < threshold:
                sdf.addLaserRangeReading(theta+brg, rng, beam_origin=(y,x))

        return sdf

    def __init__(self, initsize=(50,50), res=0.1, offset=(0,0), truncate=False):
        self.size = initsize
        self.res = res
        self.offset = offset
        self.truncate = truncate

        self._dtype = np.dtype([('dist', np.float32), ('nreadings', np.int)])

        self._array = np.zeros(self.size, dtype=self._dtype)

        self._pltextents = None

    """

    """
    def addLaserRangeReading(self, brg, rng, beam_origin=(0,0)):
        visited = set()

        (y,x) = beam_origin

        # Now we begin marching out from beam_origin
        current = 0.0
        step = self.res/5.0
        while True:
            # Take a step
            current += step

            # Where are we in the global space?
            step_y = y + (current * math.sin(brg))
            step_x = x + (current * math.cos(brg))

            # Where is the closest square to where we are right now
            idx = self._convertCoordsToLocalIdx(step_y, step_x)
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
                sq_y, sq_x = self._convertLocalIdxToCoords(idx)
                sq_to_origin = math.hypot(sq_y - y, sq_x - x)

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
        if count == 0:
            return

        # This just accumulates the averages
        ave_dist, old_count = self._array[idx]
        new_count = old_count + count
        new_ave_dist = (ave_dist * old_count + dist * count) / new_count

        if (not self.truncate) or (abs(dist) < self.truncate):
            # If we're not truncating, or the distance is less than
            # the truncated distance from the surface, add the reading
            self._array[idx] = (new_ave_dist, new_count)
            self._extendPlotExtent(idx)
        elif dist < -self.truncate:
            # Otherwise, we've gone further than the truncate distance
            # beyond the surface reading, so stop
            return "STOP"

    def _extendPlotExtent(self, idx):
        y = idx[0]
        x = idx[1]
        if self._pltextents is None:
            self._pltextents = [x,x,y,y]
        else:
            [min_x,max_x,min_y,max_y] = self._pltextents
            self._pltextents = [min(min_x,x),max(max_x,x),min(min_y,y),max(max_y,y)]



    def _convertCoordsToLocalIdx(self, y, x):
        # Think of offset as the global coordinates (in metres) of the bottom left corner
        y = int(round((y-self.offset[0])/self.res))
        x = int(round((x-self.offset[1])/self.res))

        # Of course, only convert coords if we know they're in the right space
        if (0 <= y < self.size[0]) and (0 <= x < self.size[1]):
            return (y,x)
        else:
            return None

    def _convertLocalIdxToCoords(self, idx):
        y = (idx[0] * self.res) + self.offset[0]
        x = (idx[1] * self.res) + self.offset[1]
        return (y,x)

    def extents(self):
        bl = self._convertLocalIdxToCoords((0,0))
        tr = self._convertLocalIdxToCoords((self.size[0]-1, self.size[1]-1))
        return {"bl":bl, "tr":tr}

    """
    This only works if the resolution is the same, and the offsets align based on the resolution.

    This should expand the current SDF to include all the local data, averaged into the global view.
    """
    def addLocalSDF(self, local_sdf):
        if self.res != local_sdf.res:
            raise ValueError("Resolutions of Local SDF (%f) must match Global SDF (%f)" % (local_sdf.res, self.res))

        global_extents = self.extents()
        local_extents = self.extents()
        local_in_global = (global_extents["bl"][0] <= local_extents["bl"][0] <= global_extents["tr"][0]) and \
                          (global_extents["bl"][0] <= local_extents["tr"][0] <= global_extents["tr"][0]) and \
                          (global_extents["bl"][1] <= local_extents["bl"][1] <= global_extents["tr"][1]) and \
                          (global_extents["bl"][1] <= local_extents["tr"][1] <= global_extents["tr"][1])

        if not local_in_global:
            self._extendToCover(local_extents)


        it = np.nditer(local_sdf._array, flags=['multi_index'], op_flags=['readonly'])
        while not it.finished:
            current = it[0]
            local_idx = it.multi_index
            local_coords = local_sdf._convertLocalIdxToCoords(local_idx)
            global_idx = self._convertCoordsToLocalIdx(local_coords[0], local_coords[1])
            if global_idx is None:
                raise ValueError("something's fucky, local is within global but coords don't resolve")
            else:
                self._addSurfaceDistanceReading(global_idx, current['dist'], count=current['nreadings'])

            it.iternext()

    # TODO
    def _extendToCover(self, new_extents):
        print "Extending..."
        pass

    def plot(self, name, prefix=""):
        # This is odd, but gives us really nice plots.
        masked = np.ma.masked_where(self._array["nreadings"] == 0, self._array["dist"])

        fig = plt.figure(dpi=400, tight_layout=True)
        ax = fig.add_subplot(1, 1, 1, frame_on=False, xticks=[], yticks=[], autoscale_on=True)
        ax.imshow(masked, origin='lower', interpolation='nearest', aspect='equal')
        fig.savefig("%s/%ssdf-all.%d.png" % (name, prefix, int(self.res*100)), dpi=fig.dpi)


    def plotSurfaces(self, name, prefix=""):
        # so we first mask to only the places we have readings
        masked = np.ma.masked_where(self._array["nreadings"] == 0, self._array["dist"])

        # each square has built up an average distance to a surface
        # from all the readings that have intersected with it. We want
        # to map all the squares that are within 3*self.res of that
        # surface, or at least squares that think they're within that
        # distance.
        # We'll plot the number of readings, higher number = more sure.
        surfaces = np.ma.masked_where(np.absolute(masked) > 3*self.res, masked)

        fig = plt.figure(dpi=400, tight_layout=True)
        ax = fig.add_subplot(1, 1, 1, frame_on=False, xticks=[], yticks=[], autoscale_on=True)
        ax.imshow(surfaces, origin='lower', interpolation='nearest', aspect='equal')
        fig.savefig("%s/%ssdf-surfaces.%d.png" % (name, prefix, int(self.res*100)), dpi=fig.dpi)


def closestResVal(value, res):
    return res*round(value/res)
