from mathPlane import *
from multipledispatch import *
from ursina import *
from collections import *


@dispatch(Segment, mathPlane)
def intersects(segment, plane):
    # Returns whether the given segment intersects this plane
    # if the segment isn't parallel to the plane or its start is contained in the plane
    if (compare(dotProduct(segment.vector(), plane.getNorm()), 0) != 0 or
            compare(dotProduct(segment.getStart(), plane.getNorm()), plane.getK()) == 0):
        # print(tComponents);
        # solution for t of the equation obtained by plugging in the parametrized form of the line the segment part of
        if compare(dotProduct(segment.vector(), plane.getNorm()), 0) == 0:
            #print(segment.vector())
            #print(plane.getNorm())
            #print("orth")
            return False
        t = (plane.getK() - dotProduct(plane.getNorm(), segment.getStart())) / (
            dotProduct(segment.vector(), plane.getNorm()))
        #print("t: " + str(t))
        #print(plane.getK())
        #print(dotProduct(plane.getNorm(), segment.getStart()))
        #print(dotProduct(segment.vector(), plane.getNorm()))
        #print(str(segment))
        #print(str(plane))
        if compare(t, 1) > 0:
            #print("greater")
            return False

        if compare(t, 0) < 0:
            #print("lesser")
            return False

        return True

    return False


@dispatch(Segment, mathPlane)
def intersection(segment, plane):
    # print(tComponents);
    # if the segment is parallel to the plane
    if compare(dotProduct(segment.vector(), plane.getNorm()), 0) == 0:
        # if the segments start is in the plane
        if compare(dotProduct(segment.getStart(), plane.getNorm()), plane.getK()) == 0:
            # System.out.println("!");
            return segment.getStart()

        else:
            print(",")
            return None

    # solution for t of the equation obtained by plugging in the parametrized form of the line the segment is part of
    t = (plane.getK() - dotProduct(plane.getNorm(), segment.getStart())) / (
        dotProduct(segment.vector(), plane.getNorm()))
    #print("-"+str(t)+"-")
    # print(sum(segment.getStart(), product(segment.vector(), t)));
    #print(plane)
    #print(str(segment)+"\n")
    #print(compare(t, 1))
    #print(compare(t, 0))
    # System.out.println(compare(t,BigDecimal.ZERO));
    if compare(t, 1) > 0:
        # System.out.println(".");
        return None

    if compare(t, 0) < 0:
        # System.out.println(".");
        return None
    #print(segment)
    #print(segment.vector())
    #print(t)
    return sum(segment.getStart(), product(segment.vector(), t))


class Polygon:
    def __init__Helper1(self, segments):
        self.edges = segments
        self.plane = mathPlane(self.edges[0], self.edges[1])
        self.valid = True

    def __init__Helper2(self, points):

        # Assumes the input points form a convex polygon lying within a plane
        self.edges = []
        for i in range(len(points) - 1):
            self.edges.append(Segment(points[i], points[i + 1]))

        self.edges.append(Segment(points[-1], points[0]))
        self.plane = mathPlane(self.edges[0], self.edges[1])

    def __init__(self, arg1):
        self.edges = None
        self.plane = None
        if isinstance(arg1, list):
            #print(str(arg1[0]))
            if isinstance(arg1[0], Segment):
                self.__init__Helper1(arg1)
            if isinstance(arg1[0], Point):
                self.__init__Helper2(arg1)

    def intersectedBy(self, segment):
        #print("function call")

        # determines if this polygon is intersected by the given segment
        # if the segment doesn't intersect the plane this polygon lies in,
        # then it clearly doesn't intersect this polygon
        if not intersects(segment, self.plane):
            #print("nointersection")
            return False

        # System.out.println(segment);
        # System.out.println(this.plane);
        # where the segment intersects the plane this polygon lies in
        intersect = intersection(segment, self.plane)
        # if something has gone wrong with finding the intersection, print what led to this
        if intersect == None:
            # System.out.println(segment);
            # System.out.println(this.getPlane());
            valid = False
            return False

        # System.out.println("["+intersection[0]);
        # System.out.println(intersection[1]);
        # System.out.println(intersection[2]+"]");
        all1 = True  # whether all edges checked so far returned 1 with compareTo
        all_1 = True  # whether all edges checked so far returned -1 with compareTo
        # for every edge of the polygon
        for edge in self.edges:
            val = edge.compareTo(intersect, self.plane)
            # if the point is in this edge, then it is in the entire polygon
            if val == 0:
                print("zero")
            if val == 2:
                # System.out.println("[0");
                # System.out.println(this.getPlane());
                # Main.print(intersection);
                # System.out.println(edge);
                # System.out.println("]");
                #print("2")
                return True

            # if the point is in line with this edge but not contained in it,
            # since the polygon is convex the point is not in the polygon
            if val == -2:
                #print(-2)
                return False

            if val == 1:
                # System.out.println("[1");
                # System.out.println(this.getPlane());
                # Main.print(intersection);
                # System.out.println("]");
                #print("1")
                all_1 = False

            if val == -1:
                # System.out.println("[-1");
                # System.out.println(this.getPlane());
                # Main.print(intersection);
                # System.out.println("]");
                #print("-1")
                all1 = False

            if not all1 and not all_1:
                return False

            # System.out.println("---");
            # System.out.println(val);
            # System.out.println(all1);
            # System.out.println(all_1);
            # System.out.println("---");

        # if(all1 || all_1){
        #    Main.print(intersection);
        # }
        #print(all1)
        #print(all_1)
        return all1 or all_1

    def linksWith(self, that):
        # Returns true if this polygon and the passed one have an edge that intersects the other
        thisIntersectsThat = False
        for edge in self.edges:
            if that.intersectedBy(edge):
                thisIntersectsThat = True
                # System.out.println(that.getPlane());
                # System.out.println(" with:");
                # System.out.println(edge);
                intersection(edge, that.getPlane())
                break

        thatIntersectsThis = False
        for edge in self.edges:
            if self.intersectedBy(edge):
                thatIntersectsThis = True
                # System.out.println(this.getPlane());
                # System.out.println(" with:");
                # System.out.println(edge);
                intersection(edge, self.getPlane())
                break

        return thisIntersectsThat and thatIntersectsThis

    def getEdges(self):
        return self.edges

    def setEdges(self, edges):
        self.edges = edges

    def getPlane(self):
        return self.plane

    def setPlane(self, plane):

        self.plane = plane

    def __str__(self):

        output = ""
        for i in range(len(self.getEdges())):
            output += str(self.getEdges()[i]) + "\n"

        return output + str(self.getPlane())

    def getPoints(self):
        points = []
        for i in range(len(self.edges)):
            points.append(self.edges[i].getStart())

        return points

    def signedLinksWithhelper1(self, that):
        alignedIntersections = 0
        antiAlignedIntersections = 0
        for edge in self.edges:
            if that.intersectedBy(edge):
                if that.alignmentWithNormal(edge) == 1:
                    alignedIntersections += 1
                elif that.alignmentWithNormal(edge) == -1:
                    antiAlignedIntersections += 1
        return alignedIntersections - antiAlignedIntersections

    def signedLinksWithhelper2(self, that):
        output = 0
        for triangle in that.getTriangles():
            output += self.signedLinksWith(triangle)
        return output

    def signedLinksWith(self, that):
        if isinstance(that, Polygon):
            return self.signedLinksWithhelper1(that)
        else:
            return self.signedLinksWithhelper2(that)

    @dispatch(Segment)
    def alignmentWithNormal(self, edge):
        product = dotProduct(edge.vector(), self.plane.getNorm())
        #print("alignment: "+str(product))
        if compare(product, 0) == 0:
            return 0
        elif compare(product, 0) < 0:
            return -1
        elif compare(product, 0) > 0:
            return 1

        return 2

    @dispatch(Point)
    def alignmentWithNormal(self, edge):
        product = dotProduct(edge, self.plane.getNorm())
        #print("alignment: "+str(product))
        if compare(product, 0) == 0:
            return 0
        elif compare(product, 0) < 0:
            return -1
        elif compare(product, 0) > 0:
            return 1

        return 2
    def get_intersection(self, that):
        for edge in self.edges:
            if that.intersectedBy(edge):
                return intersection(edge,that.getPlane())


@dispatch(Point, Point)
def sum(v, u):
    return Point(v.x + u.x, v.y + u.y, v.z + u.z)


@dispatch(Point, Point)
def dotProduct(v, u):
    return v.x * u.x + v.y * u.y + v.z * u.z


@dispatch(Segment, Segment)
def dotProduct(seg1, seg2):
    return dotProduct(seg1.vector(), seg2.vector())
