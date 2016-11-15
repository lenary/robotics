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

class CarmenParser(DataParser):
    def parse(self, file_in):
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
                self.rangeData.append((x,y,theta,ts, [float(x) for x in line[2:2+n]]))

            elif line[0].startswith('RLASER'):
                # likely do same as FLASER (format is same), but store in separate data structure?
                pass
