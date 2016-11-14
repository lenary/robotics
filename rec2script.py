#!/usr/bin/env python

import sys
import os

# check for DISPLAY environment variable and use matplotlib 'Agg' if not present
# this is a check for Windows Linux Subsystem GUI issues
if os.environ.get('DISPLAY') is None:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Python implementation of rec2script.c

def write_laser_data(fout, laserval1, laserval2):
    fout.write('#LASER {0} {1}:'.format(len(laserval1), len(laserval2)))
    for l in laserval1:
        fout.write(' {0}'.format(i))
    for l in laserval2:
        fout.write(' {0}'.format(i))
    fout.write('\n')

def main(fpin, fpout):
    fin = open(fpin, 'r')
    fout = open(fpout, 'w')
    
    pos = []

    for l in fin:
        line = l.rstrip()
        if line.startswith('POS'):
            _, s, us, f1, f2, f3 = line.split(' ')
            x,y,o = float(f1), float(f2), 90.0 - float(f3)
            #print('#ROBOT {0} {1} {2}'.format(x,y,o))
            pos.append((x,y,o))
        elif line.startswith('LASER-RANGE'):
            pass
            #print(line)
        elif line.startswith('LASER'):
            pass
            #print(line)

    for p in pos:
        fout.write('#ROBOT {0:.6f} {1:.6f} {2:.6f}\n'.format(p[0],p[1],p[2]))

    return pos

def plotPos(pos):
    x = [p[0] for p in pos]
    y = [p[1] for p in pos]
    plt.plot(x,y)
    plt.xlabel('x pos')
    plt.ylabel('y pos')
    plt.savefig('pos.png')

if __name__ == '__main__':
    if len(sys.argv) == 3:
        fpin = sys.argv[1]
        fpout = sys.argv[2]
        print(fpin)
        print(fpout)

        pos = main(fpin, fpout)

        plotPos(pos)
