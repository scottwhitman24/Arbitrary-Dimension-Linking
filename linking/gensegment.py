from hyperPlane import *


class genSegment:

    def __init__(self, start, end, dimension):
        self.start = start
        self.end = end

    def vector(self):
        return difference(self.start,self.end)

    def getStart(self):
        return self.start

    def setStart(self, start):
        self.start = start

    def getEnd(self):
        return self.end

    def setEnd(self, end):
        self.end = end