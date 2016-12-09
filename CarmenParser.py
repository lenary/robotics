"""
CarmenParser class

Implementation of DataParser for Carmen Log Files

CARMEN Logfile:
file format is one message per line
message_name [message contents] ipc_timestamp ipc_hostname logger_timestamp
message formats defined: PARAM SYNC ODOM FLASER RLASER
PARAM param_name param_value
SYNC tagname
ODOM x y theta tv rv accel
FLASER num_readings [range_readings] x y theta odom_x odom_y odom_theta
RLASER num_readings [range_readings] x y theta odom_x odom_y odom_theta
TRUEPOS true_x true_y true_theta odom_x odom_y odom_theta
"""

from DataParser import DataParser
import re
import math
import os

# check for DISPLAY environment variable and use matplotlib 'Agg' if not present
# this is a check for Windows Linux Subsystem GUI issues
if os.environ.get('DISPLAY') is None:
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

class CarmenParser(DataParser):
    def parse(self, file_in):
        self.posData = []
        fin = open(file_in, 'r')
        for l in fin:
            line = l.strip().split(' ')

            if line[0].startswith('ODOM'):
                x,y,theta = float(line[1]), float(line[2]), float(line[3])
                ts = float(line[-1])
                self.posData.append((x,y,theta,ts))

            elif line[0].startswith('FLASER'):
                n = int(line[1])
                x,y,theta = float(line[-9]), float(line[-8]), float(line[-7])
                ox,oy,otheta = float(line[-6]), float(line[-5]), float(line[-4])
                ts = float(line[-1])

                if n == 180 or n == 360:
                    res = 180.0/n
                    readings = [(math.radians(i*res-90), float(rng)) for i,rng in enumerate(line[2:2+n])]
                elif n == 181 or n == 361:
                    res = 180.0/(n-1)
                    readings = [(math.radians(i*res-90), float(rng)) for i,rng in enumerate(line[2:2+n])]

                entry = (x,y,theta,ts,readings)
                self.rangeData.append(entry)

            elif line[0].startswith('RLASER'):
                # likely do same as FLASER (format is same), but store in separate data structure?
                pass

            elif line[0] == 'PARAM':
                # Parameters
                param_name = line[1]
                param_val = line[2]
                if param_val == "on":
                    param_val = True
                elif param_val == "off":
                    param_val = False
                else:
                    # pay no attention to the man behind the curtain.
                    try:
                        param_val = int(param_val)
                    except ValueError:
                        try:
                            param_val = float(param_val)
                        except ValueError:
                            pass

                self.params[param_name] = param_val

        if self.posData:
            self.posOffset = (self.posData[0][0], self.posData[0][1], self.posData[0][1])
        else:
            self.posOffset = (0.0, 0.0, 0.0)


    def plotPos(self, fname):
        """
        Plot robot position from odometry data

        Position is adjusted by any initial offset in x and y, but not by any angle offset
        """
        x = [p[0]-self.posOffset[0] for p in self.posData]
        y = [p[1]-self.posOffset[1] for p in self.posData]
        plt.clf()
        plt.plot(x, y, 'k.')
        plt.axis('equal')
        plt.xlabel('x pos')
        plt.ylabel('y pos')
        plt.savefig('{0}'.format(fname))

