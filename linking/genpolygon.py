from gensegment import *


class Polygon:
    def __init__Helper1(self, segments):
        self.edges = segments

    def __init__Helper2(self, points):
        self.edges = []
        for i in range(len(points) - 1):
            self.edges.append(genSegment(points[i], points[i + 1], self.dimension))

        self.edges.append(genSegment(points[-1], points[0], self.dimension))

    def __init__(self, arg1, dimension):
        self.edges = None
        self.dimension = dimension
        if isinstance(arg1, list):
            if isinstance(arg1[0], genSegment):
                self.__init__Helper1(arg1)
            if isinstance(arg1[0], list):
                self.__init__Helper2(arg1)

    def getEdges(self):
        return self.edges

    def setEdges(self, edges):
        self.edges = edges

    def getPoints(self):
        points = []
        for i in range(len(self.edges)):
            points.append(self.edges[i].getStart())

        return points
