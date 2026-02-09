from polygon import *
from collections import *


class PointCurve:
    segments = []

    def __init__(self, points):
        self.segments = []
        for i in range(len(points) - 1):
            self.segments.append(Segment(points[i], points[i + 1]))

        self.segments.append(Segment(points[len(points) - 1], points[0]))

    def getTriangles(self):
        triangles = [Polygon([self.segments[0].getStart(), self.segments[1].getStart()
                                 , self.segments[2].getStart()])]
        for i in range(2, len(self.segments) - 1):
            previous = triangles[i - 2]
            current = Polygon([self.segments[0].getStart(), self.segments[i].getStart()
                                  , self.segments[i + 1].getStart()])
            if previous.alignmentWithNormal(current.plane.getNorm()) != 0:
                current.plane.setNorm(
                    product(current.plane.getNorm(), previous.alignmentWithNormal(current.plane.getNorm())))
            triangles.append(current)
        return triangles

    def signedLinkWith(self, that):
        output = 0
        for triangle in self.getTriangles():
            output += triangle.signedLinksWith(that)
        return output
