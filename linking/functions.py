from random import *
from math import *
from multipledispatch import *
from segment import *
from collections import *
from polygon import *

tolerance = 0.000001
Point = namedtuple('Point', ['x', 'y', 'z'])


@dispatch(Point, Point)
def dotProduct(v, u):
    return v.x * u.x+ v.y * u.y+ v.z * u.z


@dispatch(Segment, Segment)
def dotProduct(seg1, seg2):
    return dotProduct(seg1.vector(), seg2.vector())


@dispatch(Point, Point)
def crossProduct(v, u):
    return Point(v.x * u.x, v.y * u.y, v.z * u.z)


@dispatch(Segment, Segment)
def crossProduct(seg1, seg2):
    return dotProduct(seg1.vector(), seg2.vector())


@dispatch(float, float)
def compare(n, m):
    if abs(n - m) <= tolerance:
        return 0
    if n > m:
        return 1
    if n < m:
        return -1


@dispatch(Point, Point)
def difference(v, u):
    return Point(v.x - u.x, v.y - u.y, v.z - u.z)


@dispatch(Point)
def magnitude(v):
    magnitudeSquared = dotProduct(v, v)
    return sqrt(magnitudeSquared)


@dispatch(Point, Point)
def sum(v, u):
    return Point(v.x + u.x, v.y + u.y, v.z + u.z)


@dispatch(Point, float)
def product(v, k):
    return Point(v.x * k, v.y * k, v.z * k)


@dispatch(Point)
def normalize(v):
    # print(v);
    magnitudeSquared = v.x ** 2 * v.y ** 2 + v.z ** 2
    # System.out.println("|"+magnitudeSquared+"|");
    if compare(magnitudeSquared, 0) == 0:
        # System.out.println("issue/\\");
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
def same(v, u):
    # print(v1);
    # print(v2);
    return (compare(v.x, u.x) == 0 and compare(v.y, u.y) == 0
            and compare(v.z, u.z) == 0)


@dispatch(Point, Point)
def sameDirections(v, u):
    return same(normalize(v), normalize(u))


@dispatch(Segment, mathPlane)
def intersects(segment, plane):
    # Returns whether the given segment intersects this plane
    # if the segment isn't parallel to the plane or its start is contained in the plane
    if (compare(dotProduct(segment.vector(), plane.getNorm()), 0) != 0 or
            compare(dotProduct(segment.getStart(), plane.getNorm()), plane.getK()) == 0):
        # print(tComponents);
        # solution for t of the equation obtained by plugging in the parametrized form of the line the segment part of
        t = (plane.getK() - dotProduct(plane.getNorm()
                                       , segment.getStart())) / (dotProduct(segment.vector(), plane.getNorm()))
        if compare(t, 1) > 0:
            return False

        if compare(t, 0) < 0:
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

    # solution for t of the equation obtained by plugging in the parametrized form of the line the segment part of
    t = (plane.getK() - dotProduct(plane.getNorm()
                                   , segment.getStart())) / (dotProduct(segment.vector(), plane.getNorm()))
    # System.out.println("-"+t+"-");
    # print(sum(segment.getStart(), product(segment.vector(), t)));
    # System.out.println(plane);
    # System.out.println(segment);
    # System.out.println(compare(t,BigDecimal.ZERO));
    if t > 1:
        # System.out.println(".");
        return None

    if t < 0:
        # System.out.println(".");
        return None

    return sum(segment.getStart(), product(segment.vector(), t))


def randomOnSphere():
    x = Random.gauss(Random(), 0.0, 1.0)
    y = Random.gauss(Random(), 0.0, 1.0)
    z = Random.gauss(Random(), 0.0, 1.0)
    return normalize(Point(x, y, z))


# static BigDecimal[] angles(BigDecimal[] sidelengths){
#    BigDecimal angleOpposite1;
# }
@dispatch(list, list, Point, Point, float)
def triangle(sidelengths, angles, center, norm, angle):
    basePoints = [3]
    a = sidelengths[0]
    b = sidelengths[1]
    c = sidelengths[2]
    circumRadius = a * b * c / sqrt((a + b + c) * (b + c - a) * (c + a - b) * (a + b - c))
    # System.out.println(circumRadius);
    if center.x * center.x + center.y * center.y != 0:
        basePoints[0] = sum(product(Point(center.x, center.y, 0)
                                    , (center.x * center.x + center.y * center.y
                                       - circumRadius * circumRadius) / (center.x * center.x
                                                                         + center.y * center.y))
                            , Point(0, 0, center.z))
    else:
        basePoints[0] = Point(sqrt(2) / 2 * circumRadius, sqrt(2) / 2 * circumRadius, center.z)

    basePoints[1] = rotateInXY(basePoints[0], center, angles[1])
    basePoints[2] = rotateInXY(basePoints[1], center, angles[2])
    # System.out.println("Given sidelengths");
    # print(sidelengths);
    # System.out.println("angles");
    # print(angles);
    # System.out.println("center");
    # print(center);
    # System.out.println("norm");
    # print(norm);
    # System.out.println("angle\n"+angle);
    # System.out.println("Starts as");
    # print(basePoints[0]);
    # print(basePoints[1]);
    # print(basePoints[2]);
    # System.out.println();
    basePoints[0] = rotateInXY(basePoints[0], center, angle)
    basePoints[1] = rotateInXY(basePoints[1], center, angle)
    basePoints[2] = rotateInXY(basePoints[2], center, angle)
    # System.out.println("Rotated by "+angle+" radians it becomes");
    # print(basePoints[0]);
    # print(basePoints[1]);
    # print(basePoints[2]);
    return rotateTo(Polygon(basePoints), norm)


@dispatch(Polygon, Point)
def rotateTo(polygon, norm):
    unitNormal = normalize(polygon.getPlane().getNorm())
    newUnitNormal = normalize(norm)
    mag = magnitude(crossProduct(unitNormal, newUnitNormal))
    # System.out.println("Magnitude\n"+mag);
    if compare(mag, 0) > 0:
        rotationAxis = product(crossProduct(unitNormal, newUnitNormal), 1.0 / mag)
    else:
        rotationAxis = unitNormal

    rotationAngle = acos(dotProduct(unitNormal, newUnitNormal))
    points = polygon.getPoints()
    newPoints = [points.length][3]
    # Rodrigues' rotation formula
    for i in range(points.length + 1):
        intermed1 = sum(product(points[i], cos(rotationAngle))
                        , product(crossProduct(rotationAxis, points[i]), sin(rotationAngle)))
        k = dotProduct(rotationAxis, points[i]) * (1 - cos(rotationAngle))
        intermed2 = product(rotationAxis, k)
        # print(rotationAxis);
        # System.out.println("|");
        # print(intermed1);
        # System.out.println("|");
        # print(intermed2);
        newPoints[i] = sum(intermed1, intermed2)
        # System.out.println("=");
        # print(newPoints[i]);

    return Polygon(newPoints)


@dispatch(Point, Point, float)
def rotateInXY(basePoint, center, angle):
    # System.out.println("Given basepoint");
    # print(basePoint);
    # System.out.println("Center");
    # print(center);
    # System.out.println("Angle\n"+angle);
    radius = magnitude(difference(basePoint, center))
    otherAngle = (((basePoint.y - center.y) / abs(basePoint.y - center.y)) *
                  acos((basePoint.x - center.x) / sqrt(
                      (basePoint.x - center.x) ** 2 + (basePoint.y - center.y) ** 2)))
    # System.out.println("Arccos of "+(basePoint[0]-center[0])/Math.sqrt(
    # (basePoint[0]-center[0])*(basePoint[0]-center[0])+(basePoint[1]-center[1])*(basePoint[1]-center[1]))+" is");
    # System.out.println(otherAngle);
    diff = Point(radius * (cos(2.0 * angle + otherAngle))
                 , radius * sin(2.0 * angle + otherAngle), basePoint.z)
    # System.out.println("Gets diff");
    # print(diff);
    # System.out.println("Returns");
    # print(sum(center,diff));
    return sum(center, diff)


@dispatch(list)
def angle(sideLengths):
    return acos((sideLengths[0] * sideLengths[0] - sideLengths[1] * sideLengths[1]
                 - sideLengths[2] * sideLengths[2]) / (2 * sideLengths[1] * sideLengths[2]))


@dispatch(float)
def randomTriangle(radius):
    point1 = randomOnSphere()
    point2 = randomOnSphere()
    point3 = randomOnSphere()
    sideLengths = [1 - point1.x ** 2, 1 - point1.y ** 2, 1 - point1.z ** 2]
    angles = [angle(sideLengths), angle([sideLengths[1], sideLengths[2], sideLengths[0]])
        , angle([sideLengths[1], sideLengths[2], sideLengths[0]])]
    return triangle(sideLengths, angles, product(point2, radius)
                    , point3, 2 * pi * Random.random(Random()))


@dispatch(float)
def randomTriangleOrigin(radius):
    point1 = randomOnSphere()
    point2 = randomOnSphere()
    point3 = randomOnSphere()
    sideLengths = [1 - point1.x ** 2, 1 - point1.y ** 2, 1 - point1.z ** 2]
    angles = [angle(sideLengths), angle([sideLengths[1], sideLengths[2], sideLengths[0]])
        , angle([sideLengths[1], sideLengths[2], sideLengths[0]])]
    return triangle(sideLengths, angles, Point(0, 0, 0), Point(0, 0, 1), 0)
