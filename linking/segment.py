from multipledispatch import *
from collections import *
from math import *

tolerance = 0.000001
Point = namedtuple('Point', ['x', 'y', 'z'])


class Segment:
    start = Point(0, 0, 0)
    end = Point(0, 0, 0)

    def __init__(self, start, end):
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
        cross = crossProduct(toPoint, self.vector())
        #print("point is "+str(point))
        #print("topoint is "+str(toPoint))
        #print("segment is "+str(self))
        #print("vector is "+str(self.vector()))
        #print("cross is "+str(cross))
        #print("norm is "+str(plane.getNorm()))
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


@dispatch(Segment, Segment)
def dotProduct(seg1, seg2):
    return dotProduct(seg1.vector(), seg2.vector())


@dispatch(Point, Point)
def crossProduct(v, u):
    return Point(v.y * u.z-v.z*u.y, v.z * u.x-v.x*u.z, v.x*u.y - v.y*u.x)


@dispatch(Segment, Segment)
def crossProduct(seg1, seg2):
    return crossProduct(seg1.vector(), seg2.vector())


@dispatch(Point, Point)
def dotProduct(v, u):
    return v.x * u.x+v.y * u.y+v.z * u.z


@dispatch(Point)
def magnitude(v):
    magnitudeSquared = dotProduct(v, v)
    return sqrt(magnitudeSquared)


@dispatch(Point)
def normalize(v):
    # print(v);
    magnitudeSquared = v.x ** 2 + v.y ** 2 + v.z ** 2
    # System.out.println("|"+magnitudeSquared+"|");
    if compare(magnitudeSquared, 0) == 0:
        raise ZeroDivisionError
        return v

    inverseMagnitude = sqrt(1 / magnitudeSquared)
    # System.out.println("\\"+magnitudeSquared+"/");
    # System.out.println("/"+inverseMagnitude+"\\");
    # System.out.println("v");
    # print(v);
    output = product(v, inverseMagnitude)
    # System.out.println("output");
    # print(output);
    return output

@dispatch(list)
def normalize(v):
    # print(v);
    magnitudeSquared = 0
    for element in v:
        magnitudeSquared+=element**2
    # System.out.println("|"+magnitudeSquared+"|");
    if compare(magnitudeSquared, 0) == 0:
        raise ZeroDivisionError
        return v

    inverseMagnitude = sqrt(1 / magnitudeSquared)
    # System.out.println("\\"+magnitudeSquared+"/");
    # System.out.println("/"+inverseMagnitude+"\\");
    # System.out.println("v");
    # print(v);
    output = product(v, inverseMagnitude)
    # System.out.println("output");
    # print(output);
    return output


@dispatch(Point, Point)
def difference(v, u):
    return Point(u.x - v.x, u.y - v.y, u.z - v.z)


@dispatch(Point, Point)
def same(v, u):
    # print(v1);
    # print(v2);
    return (compare(v.x, u.x) == 0 and compare(v.y, u.y) == 0
            and compare(v.z, u.z) == 0)


def compare(n, m):
    if abs(n - m) <= tolerance:
        return 0
    if n > m:
        return 1
    if n < m:
        return -1


@dispatch(Point, Point)
def sameDirections(v, u):
    return same(normalize(v), normalize(u))

@dispatch(Point,object)
def product(v, k):
    return Point(v.x * k, v.y * k, v.z * k)
@dispatch(list,object)
def product(v, k):
    return [v[i]*k for i in range(len(v))]