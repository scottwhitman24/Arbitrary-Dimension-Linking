from multipledispatch import *
from collections import *
from math import *
import numpy as np

tolerance = 0.000001


class HyperPlane:

    def __init__(self, points, dimension):
        self.norm = crossProduct([difference(points[i], points[0]) for i in range(1, len(points))], dimension)
        self.k = dotProduct(self.norm, points[0])
        self.dimension = dimension  # dimension of ambient space

    def vector(self):
        return difference(self.start, self.end)

    def getStart(self):
        return self.start

    def setStart(self, start):
        self.start = start

    def getEnd(self):
        return self.end

    def setEnd(self, end):
        self.end = end

    def getNorm(self):
        return self.norm

    def setNorm(self, norm):
        normsum = sum(self.norm)
        self.norm = norm
        if compare(self.k, 0) == 0:
            self.k = dotProduct(self.norm, [0 for i in range(self.dimension)])
            return
        self.k = dotProduct(self.norm, [self.k / normsum for i in range(self.dimension)])

    def getK(self):
        return self.k

    def setK(self, k):
        self.k = k

    @dispatch(list)
    def contains(self, point):
        # determines if the input point is contained in the plane
        # System.out.println("Contains:");
        # System.out.println(Main.dotProduct(point, this.getNorm()));
        # System.out.println(this.getK());
        return compare(dotProduct(point, self.getNorm()), self.getK()) == 0

    @dispatch(object)
    def contains(self, segment):
        # determines if the input segment is contained in the plane
        # print("Contains:\n"+this+"\n"+segment)
        orthogonality = compare(dotProduct(segment.vector(), self.getNorm()), 0)
        startContained = compare(dotProduct(segment.getStart(), self.getNorm()), self.getK())
        # print(orthogonality);
        # print(startContained);
        return orthogonality == 0 and startContained == 0

    def __str__(self):
        return "Norm=" + str(self.norm) + " k=" + str(self.k)

    def compareTo(self, point, plane):
        # Returns a 1 if the given point is strictly to the right of this segment
        # as determined by the orientation induced by the normal of the plane,
        # a -1 if the point is strictly to the left, a 2 if it is contained in the segment and
        # a -2 if it is in line with the segment but not contained by it*/
        # Main.print(Main.crossProduct(Main.difference(
        #        point,this.getStart()),this.vector()));
        # Main.print(plane.getNorm());
        if not plane.contains(point):
            print(str(self) + " is -3 with:")
            print(point)
            return -3
        toPoint = difference(point, self.getStart())  # the vector from start to point
        # if point is equal to start, return 2
        if compare(magnitude(toPoint), 0) == 0:
            # System.out.println(this+" is 2.1 with:");
            # Main.print(point);
            return 2

        if sameDirections(toPoint, self.vector()):
            # Main.print(Main.difference(point,this.getStart()));
            # Main.print(this.vector());
            # if toPoint is not larger than this.vector(), then point is inside this segment, and otherwise it isn't
            if compare(magnitude(toPoint), magnitude(self.vector())) <= 0:
                # System.out.println(this+" is 2.2 with:");
                # Main.print(point);
                return 2
            else:
                # System.out.println(this+" is -2.1 with:");
                # Main.print(point);
                return -2

        # if point-start is in the opposite direction of this.vector() return -2
        invertedToPoint = product(toPoint, -1)
        if sameDirections(invertedToPoint, self.vector()):
            # System.out.println(this+" is 2.2 with:");
            # Main.print(point);
            # Main.print(Main.difference(point,this.getStart()));
            # Main.print(this.vector());
            return -2
        # if (point-start)x(end-start) is in the same direction as the normal to the plane,
        # point is to the right of this segment with respect to this plane, so return 1
        print("cross")
        cross = crossProduct(toPoint, self.vector())
        # print("point is "+str(point))
        # print("topoint is "+str(toPoint))
        # print("segment is "+str(self))
        # print("vector is "+str(self.vector()))
        # print("cross is "+str(cross))
        # print("norm is "+str(plane.getNorm()))
        if sameDirections(cross, plane.getNorm()):
            # System.out.println(this+" is 1 with:");
            # Main.print(point);
            # System.out.println("toPoint");
            # Main.print(toPoint);
            # System.out.println("vector");
            # Main.print(this.vector());
            # System.out.println("cross");
            # Main.print(cross);
            # System.out.println("norm");
            # Main.print(plane.getNorm());
            return 1

        # if (point-start)x(end-start) is in the opposite direction as the normal to the plane,
        # then point is to the left of this segment with respect to this plane, so return 1
        invertedCross = product(cross, -1)
        if sameDirections(invertedCross, plane.getNorm()):
            # System.out.println(this+" is -1 with:");
            # Main.print(point);
            return -1

        # if we have gotten this far, either the segment is not in the plane or
        # # we have failed to determine which side the point is on, so pass a unique value
        # print(str(self)+" is 0 with:")
        # print(point)
        # print("toPoint")
        # print(toPoint)
        # print("vector")
        # print(self.vector())
        # print("cross")
        # print(cross)
        # print("invertedCross")
        # print(invertedCross)
        # print("norm")
        # print(plane.getNorm())
        # print("normalized cross")
        # print(normalize(cross))
        # print("normalized invertedCross")
        # print(normalize(invertedCross))
        # print("normalized norm")
        # print(normalize(plane.getNorm()))
        return 0

    def __str__(self):

        return ("(" + str(self.getStart().x) + "," + str(self.getStart().y) + "," + str(self.getStart().z) + ")->"
                + "(" + str(self.getEnd().x) + "," + str(self.getEnd().y) + "," + str(self.getEnd().z) + ")")


def basisElement(index, dimension):
    return [i == index - 1 for i in range(dimension)]


def crossProduct(vectors, dimension):
    cross = [0 for i in range(dimension)]
    varray = np.array(vectors)
    axis=1
    if varray.shape[0]>=varray.shape[1]:
        axis=0
    for i in range(dimension):
        cross[i] = ((-1) ** dimension) * vectors[0][i] * np.linalg.det(np.delete(varray, i, axis))
    return cross


def dotProduct(v, u):
    return sum([v[i] * u[i] for i in range(len(v))])


def magnitude(v):
    magnitudeSquared = dotProduct(v, v)
    return sqrt(magnitudeSquared)


def normalize(v):
    magnitudeSquared = dotProduct(v, v)
    if compare(magnitudeSquared, 0) == 0:
        raise ZeroDivisionError

    inverseMagnitude = sqrt(1 / magnitudeSquared)
    # System.out.println("\\"+magnitudeSquared+"/");
    # System.out.println("/"+inverseMagnitude+"\\");
    # System.out.println("v");
    # print(v);
    output = product(v, inverseMagnitude)
    # System.out.println("output");
    # print(output);
    return output


def difference(v, u):
    return [v[i] - u[i] for i in range(len(v))]


def same(v, u):
    return all([compare(v[i], u[i]) for i in range(len(v))])


def compare(n, m):
    if abs(n - m) <= tolerance:
        return 0
    if n > m:
        return 1
    if n < m:
        return -1


def sameDirections(v, u):
    return same(normalize(v), normalize(u))


def product(v, k):
    return [v[i] * k for i in range(len(v))]
