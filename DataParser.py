"""
DataParser base class

All DataParsers should implement the parse method that reads a data log file
and populates the position and laser range data lists: posData, rangeData

Position Data: (x, y, theta, timestamp)
Range Data: (x, y, theta, timestamp, [laser range readings])
"""
class DataParser:
    def __init__(self):
        self.posData = []
        self.rangeData = []

    def parse(self, file_in):
        return NotImplemented
