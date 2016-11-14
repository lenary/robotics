#!/usr/bin/env python

import sys
import os

from IntelParser import IntelParser

# check for DISPLAY environment variable and use matplotlib 'Agg' if not present
# this is a check for Windows Linux Subsystem GUI issues
if os.environ.get('DISPLAY') is None:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

def plotPos(pos):
    x = [p[0] for p in pos]
    y = [p[1] for p in pos]
    plt.plot(x,y)
    plt.xlabel('x pos')
    plt.ylabel('y pos')
    plt.savefig('pos.png')

if __name__ == '__main__':
    if len(sys.argv) == 2:
        fpin = sys.argv[1]
        print(fpin)

        ip = IntelParser()
        ip.parse(fpin)

        plotPos(ip.posData)
