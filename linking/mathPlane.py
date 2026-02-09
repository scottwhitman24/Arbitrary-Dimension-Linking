from multipledispatch import *
from segment import *
from collections import *


@dispatch(Point, Point)
def dotProduct(v, u):
    return Point(v.x * u.x, v.y * u.y, v.z * u.z)


@dispatch(Segment, Segment)
def dotProduct(seg1, seg2):
    return dotProduct(seg1.vector(), seg2.vector())


class mathPlane:
    norm = Point(0, 0, 0)  # the normal vector to the plane
    k = 0  # the constant on the right hand side of the equation induced

    # by the normal and a point on the plane
    @dispatch(Point, float)
    def __init__(self, n, k):
        self.norm = n
        self.k = k

    @dispatch(Segment, Segment)
    def __init__(self, segment1, segment2):
        #print("segments")
        #print(segment1)
        #print(segment2)
        # Constructor that takes as input two segments lying inside the plane
        self.norm = crossProduct(segment1.vector(), segment2.vector())
        #print(self.norm)
        self.k = dotProduct(self.norm, segment1.getStart())
        #print(self)

    def getNorm(self):
        return self.norm

    def setNorm(self, norm):
        sum = self.norm.x + self.norm.y + self.norm.z
        self.norm = norm
        if compare(self.k, 0) == 0:
            self.k = dotProduct(self.norm, Point(0, 0, 0))
            return
        self.k = dotProduct(self.norm, Point(self.k / sum, self.k / sum, self.k / sum))

    def getK(self):
        return self.k

    def setK(self, k):
        self.k = k

    @dispatch(Segment)
    def contains(self, segment):
        # determines if the input segment is contained in the plane
        # print("Contains:\n"+this+"\n"+segment)
        orthogonality = compare(dotProduct(segment.vector(), self.getNorm()), 0)
        startContained = compare(dotProduct(segment.getStart(), self.getNorm()), self.getK())
        # print(orthogonality);
        # print(startContained);
        return orthogonality == 0 and startContained == 0

    @dispatch(Point)
    def contains(self, point):
        # determines if the input segment is contained in the plane
        # System.out.println("Contains:");
        # System.out.println(Main.dotProduct(point, this.getNorm()));
        # System.out.println(this.getK());
        return compare(dotProduct(point, self.getNorm()), self.getK()) == 0

    def __str__(self):
        return (str(self.getNorm()[0]) + "x+" + str(self.getNorm()[1])
                + "y+" + str(self.getNorm()[2]) + "z=" + str(self.getK()))
