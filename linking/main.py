import random

import numpy as np
import numpy.linalg
import ursina.color
from ursina import *

from orientedBoundedHyperPlane import orientedBoundedHyperPlane
from polygon import *
from pointcurve import *
from multipledispatch import *
from collections import *
from random import *

tolerance = 0.000001

popup_text = None
intersect_text = None
visible = False

this = None
that = None


def randomOnSphere():
    x = Random.gauss(Random(), 0.0, 1.0)
    y = Random.gauss(Random(), 0.0, 1.0)
    z = Random.gauss(Random(), 0.0, 1.0)
    return normalize(Point(x, y, z))


def randomOnHyperSphere(dimension):
    while True:
        coords = [0.0 for i in range(dimension)]
        for i in range(dimension):
            coords[i] = Random.gauss(Random(), 0.0, 1.0)
        try:
            out = normalize(coords)
        except(ZeroDivisionError):
            continue
        return out


# static BigDecimal[] angles(BigDecimal[] sidelengths){
#    BigDecimal angleOpposite1;
# }
@dispatch(list, list, Point, Point, float)
def triangle(sidelengths, angles, center, norm, angle):
    basePoints = [0, 0, 0]
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
    newPoints = [[0, 0, 0] for i in range(len(points))]
    # Rodrigues' rotation formula
    for i in range(len(points)):
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
    if abs((sideLengths[0] * sideLengths[0] - sideLengths[1] * sideLengths[1]
            - sideLengths[2] * sideLengths[2]) / (2 * sideLengths[1] * sideLengths[2])) > 1:
        print("issue")
        print(sideLengths)
        print((sideLengths[0] * sideLengths[0] - sideLengths[1] * sideLengths[1]
               - sideLengths[2] * sideLengths[2]) / (2 * sideLengths[1] * sideLengths[2]))
        raise ZeroDivisionError
    try:
        return acos((sideLengths[0] * sideLengths[0] - sideLengths[1] * sideLengths[1]
                     - sideLengths[2] * sideLengths[2]) / (2 * sideLengths[1] * sideLengths[2]))
    except ValueError:
        print("uncaught issue")
        print(sideLengths)
        print((sideLengths[0] * sideLengths[0] - sideLengths[1] * sideLengths[1]
               - sideLengths[2] * sideLengths[2]) / (2 * sideLengths[1] * sideLengths[2]))
        raise ZeroDivisionError


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


class ursina_polygon(Button):
    def __init__(self, mesh, color, polygon):
        super().__init__(
            parent=scene,
            model=mesh,
            collider='box',
            rotation=(0, 0, 0),
            color=color,
            scale=(1, 1, 1))
        self.polygon = polygon
        self.infotext = Text(text="None", position=(0.6, 0.5, 0))

    def input(self, key):
        global this
        if self.hovered:
            if key == 'left mouse down':
                if self.alpha_getter() == 1:
                    self.alpha_setter(0.1)
                else:
                    self.alpha_setter(1)
            elif key == 'right mouse down':
                global visible
                if visible:
                    # destroy(popup_text)
                    # popup_text = Text(text=str(self.polygon), position=(0.6, 0.5, 0))
                    self.infotext.text = str(self.polygon)
                    visible = True
                else:
                    # popup_text = Text(text=str(self.polygon), position=(0.6, 0.5, 0))
                    visible = True
                    # destroy(popup_text, delay=3)
                    self.infotext.text = str(self.polygon)
                    visible = False
            elif key == 's':
                this = self.polygon
            elif key == 'e':
                global that
                that = self.polygon
                intersect = this.get_intersection(that)
                # intersect_text = Text(text=str(intersect), position=(0.5, 0.4, 0))
                self.infotext.text = str(intersect)
                if intersect != None:
                    ball = Entity(model='sphere', scale=0.2, position=(intersect.x, intersect.y, intersect.z))
                    destroy(ball, delay=3)
                # destroy(intersect_text, delay=3)


def makeMesh(polygon):
    vertices = [[p.x, p.y, p.z] for p in polygon.getPoints()]
    face = []
    facer = []
    for i in range(len(vertices)):
        face.append(i)
        facer.append(len(vertices) - i - 1)
    triangles = [face, facer]
    return [vertices, triangles]


def makeMeshs(pointcurve):
    meshs = []
    for triangle in pointcurve.getTriangles():
        meshs.append(makeMesh(triangle))
    return meshs





def adjoin(boundary, interior, vertex, points=None):
    c = choice(range(len(boundary)))
    # face = choice(boundary)
    face = boundary[c]
    # boundary.remove(face)
    boundary = boundary[:c] + boundary[c + 1:]
    interior.append([face[0][:face[1]] + [vertex] + face[0][face[1]:], -face[2]])
    boundary.extend([[face[0][:j] + [vertex] + face[0][j + 1:], j, interior[-1][1]] for j in range(len(face[0]))])
    if not points is None:
        points.append(vertex)
    return (boundary, interior)


def linking(boundary, interior, dimension):
    link = 0
    for simplex in interior:
        for face in boundary:
            link += orientedBoundedHyperPlane(simplex[0], simplex[1], dimension).signedLinksWith(
                orientedBoundedHyperPlane(face[0], face[2] ** face[1], dimension))
    return link


def sphereApprox(center, dimension, ambientdimension, flag, scale=1.0):
    if flag == 2:
        hedron = (np.diag([1 * scale for i in range(ambientdimension)]).tolist()[:dimension])
        hedron.extend(np.diag([-1 * scale for i in range(ambientdimension)]).tolist()[:dimension])
    else:
        hedron = (np.diag([1 * scale for i in range(ambientdimension)]).tolist()[ambientdimension - dimension:])
        hedron.extend(np.diag([-1 * scale for i in range(ambientdimension)]).tolist()[ambientdimension - dimension:])
    directions = [[0, 0] for j in range(dimension-1)]
    for j in range(dimension-1):
        if flag == 1:
            zero = [0 for k in range(j + ambientdimension - dimension)]
            zero.extend(randomOnHyperSphere(dimension - j))
            directions[j] = [np.array(hedron[j]), np.array(zero)]
        else:
            zero = randomOnHyperSphere(dimension - j)
            zero.extend([0 for k in range(j + ambientdimension - dimension)])
            directions[j] = [np.array(hedron[j]), np.array(zero)]
    rotation = np.identity(ambientdimension)
    for j in range(dimension-1):
        temp = getRotationMatrix(directions[j][0].tolist(), directions[j][1].tolist(), ambientdimension)
        rotation = np.linalg.matmul(temp, rotation)
        for i in range(dimension-1):
            directions[i][0] = np.linalg.matmul(temp, directions[i][0])
            directions[i][1] = np.linalg.matmul(temp, directions[i][1])
    if abs(np.linalg.det(rotation))>2.0:
        print("rotation")
        print(np.linalg.det(rotation))
        print(rotation)
    for j in range(2 * dimension):
        hedron[j] = np.linalg.matmul(rotation, hedron[j])
    hedronlist = [hedron[i] for i in range(len(hedron))]
    hedronlist = [[hedronlist[i][j] + center[j] for j in range(len(hedronlist[i]))] for i in range(len(hedronlist))]
    hedronInterior = [[hedronlist[:dimension + 1], 1]]
    hedronBoundary = [[hedronInterior[0][0][:j] + hedronInterior[0][0][j + 1:], j, hedronInterior[-1][1]] for j in
                      range(dimension + 1)]
    for i in range(dimension + 1):
        (hedronBoundary, hedronInterior) = adjoin(hedronBoundary, hedronInterior, hedronlist[dimension + i - 1])
    return (hedronBoundary, hedronInterior)


def getRotationMatrix(initial, final, dimension):
    rotation = np.identity(dimension)
    x = np.array(initial)
    x.shape = (1, dimension)
    y = np.array(final)
    y.shape = (1, dimension)
    if np.matmul(x, x.transpose())[0][0] != 0:
        u = x * (1 / sqrt(np.matmul(x, x.transpose())[0][0]))
    else:
        u = x
    u.shape = (1, dimension)
    v = y - (np.matmul(np.matmul(u, y.transpose()), u))
    # v=y
    if np.matmul(v, v.transpose()) != 0:
        v = v * (1 / sqrt(np.matmul(v, v.transpose())[0][0]))
    v.shape = (1, dimension)
    uv = np.array([u.tolist()[0], v.tolist()[0]])
    if np.matmul(u, u.transpose()) * np.matmul(v, v.transpose()) != 0:
        cosine = (np.matmul(x, y.transpose())[0][0] / sqrt(
            np.matmul(x, x.transpose())[0][0] * np.matmul(y, y.transpose())[0][0]))
    else:
        cosine = 1
    sine = sqrt(1 - min(abs(cosine), 1.0) ** 2)
    twod = [[cosine, -sine], [sine, cosine]]
    rotation = rotation - np.matmul(u.transpose(), u) - np.matmul(v.transpose(), v) + np.matmul(uv.transpose(),
                                                                                                np.matmul(
                                                                                                    np.array(twod), uv))
    return rotation


def spheresim(dim1, dim2, dimension, reps):
    linkingsum = 0
    for i in range(reps):
        center1 = [0 for i in range(dimension - dim1)]
        center1.extend(randomOnHyperSphere(dim1))
        center1 = [center1[i] * 2 for i in range(len(center1))]
        center2 = randomOnHyperSphere(dim2)
        center2.extend([0 for i in range(dimension - dim2)])
        center2 = [center2[i] * 2 for i in range(len(center2))]
        # print(center1)

        inout = 0
        for i in range(400):
            hedron1 = sphereApprox(center1, dim1, dimension, 1)
            hedron2 = sphereApprox(center2, dim2, dimension, 2, sqrt(dim2))
            if linking(hedron1[0], hedron2[1], dimension) == 1:
                inout = 1
                # print("link")
                break
        outin = 0
        for i in range(400):
            hedron1 = sphereApprox(center1, dim1, dimension, 1, sqrt(dim1))
            hedron2 = sphereApprox(center2, dim2, dimension, 2)
            if linking(hedron1[0], hedron2[1], dimension) == 1:
                outin = 1
                # print("link")
                break
        linkingNum = 0
        if inout == 1 and outin == 1:
            for i in range(250):
                hedron1 = sphereApprox(center1, dim1, dimension, 1)
                hedron2 = sphereApprox(center2, dim2, dimension, 2)
                if linking(hedron1[0], hedron2[1], dimension) == 1:
                    linkingNum = 1
                    # print("link")
                    break
        linkingsum += linkingNum
    return linkingsum

def hedronsim(dim1, dim2, dimension, reps):
    linkingsum = 0
    center1 = [0 for i in range(dimension)]
    center2 = [0 for i in range(dimension)]
    center2[0]=1
    # print(center1)

    linkingsum=0
    for i in range(reps):
        linkingNum=0
        hedron1 = sphereApprox(center1, dim1, dimension, 1, 0.9)
        print(hedron1)
        hedron2 = sphereApprox(center2, dim2, dimension, 2, 0.9)
        if linking(hedron1[0], hedron2[1], dimension) == 1:
            linkingNum=1
        linkingsum += linkingNum
    return linkingsum




# points = [[1 if j == i else 0 for j in range(4)] for i in range(4)]
'''points=[[-2,-2,-2],[1,1,1],[2,-1,0]]
boundary = [[points[:j] + points[j + 1:],j, 1] for j in range(len(points))]
points=[[1,0,0],[0,1,0],[0,0,1]]
interior = [[points[:], 1]]
#adjoin(points, boundary, interior, [1, 1, 1, 1])
print("Linking of")
print(boundary)
print("and")
print(interior)
print("is")
print(linking(boundary, interior, 3))'''

'''def sinint(n):
    if n%2==1:
        sum=0
        for i in range(n+1):
            if i%2==1:
                sum+=-math.comb(n,i)/(2*i-n)
            else:
                sum+=math.com(n,i)/(2*i-n)
    else:
        return math.pi/(2**(n+1))*math.comb(n,i)'''

'''vals=[0 for i in range(100)]
for i in range(100):
    vals[i]=spheresim(3,4,6,10)
    print(vals[i])
print("sum")
sum=0
for i in range(100):
    sum+=vals[i]
print(sum)
# print(spheresimupper(2,2,3,100))'''
print(hedronsim(2,2,3,1000))



'''center1 = [0 for i in range(dimension - dim1)]
center1.extend(randomOnHyperSphere(dim1))
center1 = [center1[i] * 2 for i in range(len(center1))]
center2 = randomOnHyperSphere(dim2)
center2.extend([0 for i in range(dimension - dim2)])
center2 = [center2[i] * 2 for i in range(len(center2))]
center1 = [0, 0, 0]
center2 = [0, 2.1, 0]
inout = 0
for k in range(10):
    for i in range(300):
        hedron1 = sphereApprox(center1, dim1, dimension, 1, 1.0)
        hedron2 = sphereApprox(center2, dim2, dimension, 2, sqrt(dim2))
        if linking(hedron1[0], hedron2[1], dimension) == 1:
            inout = 1
            # print("link")
            break
    outin = 0
    for i in range(300):
        hedron1 = sphereApprox(center1, dim1, dimension, 1, sqrt(dim1))
        hedron2 = sphereApprox(center2, dim2, dimension, 2, 1.0)
        if linking(hedron1[0], hedron2[1], dimension) == 1:
            outin = 1
            # print("link")
            break
    linkingNum = 0
    if inout == 1 and outin == 1:
        print("both")
        for i in range(100):
            hedron1 = sphereApprox(center1, dim1, dimension, 1, 1.0)
            hedron2 = sphereApprox(center2, dim2, dimension, 2, 1.0)
            if linking(hedron1[0], hedron2[1], dimension) == 1:
                linkingNum = 1
                print("link")
                break
    print(linkingNum)
'''
'''
app = Ursina()
# mesh = Mesh(vertices=[[0, 0, 0], [0, 0, 1], [2, -1, 0]], triangles=[[0, 1, 2], [2, 1, 0]])
# mesh = makeMesh(Polygon([Point(0, 0, 0), Point(0, 0, 1), Point(2, -1, 0)]))
# entity = Entity(model=Mesh(vertices=[[0,0,0],[0,0,1],[2,-1,0]],triangles=[[0,1,2]]))
# entity = Entity(model=Mesh(vertices=mesh[0], triangles=mesh[1]))
p1 = PointCurve([Point(0.5, 0, 0), Point(0, 0, 1), Point(2, -1, 0), Point(-0.3, 2, 0)])
t1 = Polygon([Point(0.5, 0, 0), Point(2, -1, 0), Point(-0.3, 2, 0)])
meshs = makeMeshs(p1)
# entities = []
polygons = []
colors = [ursina.color.white, ursina.color.red, ursina.color.blue, ursina.color.cyan]
i = 0
for mesh in meshs:
    # entities.append(Entity(model=Mesh(mesh[0], mesh[1]), color=colors[i % (len(colors))]))
    ursina_polygon(Mesh(mesh[0], mesh[1]), colors[i % (len(colors))], p1.getTriangles()[meshs.index(mesh)])
    i += 1
p2 = PointCurve([Point(0, 0, -1), Point(0, 0, -3), Point(0, 2, 1), Point(-3, 0, -1)])
t2 = Polygon([Point(0, 0, -1), Point(0, 0, -3), Point(0, 2, 1)])
meshs = makeMeshs(p2)
for mesh in meshs:
    # entities.append(Entity(model=Mesh(mesh[0], mesh[1]), color=colors[i % (len(colors))]))
    ursina_polygon(Mesh(mesh[0], mesh[1]), colors[i % (len(colors))], p2.getTriangles()[meshs.index(mesh)])
    i += 1
points = [Vec3(-5, 0, 0), Vec3(5, 0, 0)]
entity = Entity(model=Mesh(vertices=points, mode='line'), color=colors[(len(colors) - 3)])
# entities.append(entity)
points = [Vec3(0, -5, 0), Vec3(0, 5, 0)]
entity2 = Entity(model=Mesh(vertices=points, mode='line'), color=colors[(len(colors) - 2)])
# entities.append(entity)
points = [Vec3(0, 0, -5), Vec3(0, 0, 5)]
entity3 = Entity(model=Mesh(vertices=points, mode='line'), color=colors[(len(colors) - 1)])
# entities.append(entity)
print(p1.signedLinkWith(p2))
print(p2.signedLinkWith(p1))
EditorCamera()
app.run()
'''
'''app = Ursina()
flag=True
for r in range(1, 50):
    signed = 0
    absolute = 0
    for i in range(100000):
        try:
            triangle1 = randomTriangle(1+r / 10)
            triangle2 = randomTriangle(1+r / 10)
            linking=triangle1.signedLinksWith(triangle2)
            signed += linking
            absolute += abs(linking)
            if linking==0:
                print("polys")
                ursina_polygon(makeMesh(triangle1),ursina.color.red, triangle1)
                ursina_polygon(makeMesh(triangle2),ursina.color.blue, triangle2)
                flag=False
                break
        except ZeroDivisionError:
            continue
    if flag==False:
        break
    print(str(5+r/10)+" signed: " + str(signed) + " "+str(5+r/10)+" absolute: " + str(absolute)+" "+str(5+r/10))
points = [Vec3(-5, 0, 0), Vec3(5, 0, 0)]
entity = Entity(model=Mesh(vertices=points, mode='line'), color=colors[(len(colors) - 3)])
# entities.append(entity)
points = [Vec3(0, -5, 0), Vec3(0, 5, 0)]
entity2 = Entity(model=Mesh(vertices=points, mode='line'), color=colors[(len(colors) - 2)])
# entities.append(entity)
points = [Vec3(0, 0, -5), Vec3(0, 0, 5)]
entity3 = Entity(model=Mesh(vertices=points, mode='line'), color=colors[(len(colors) - 1)])
EditorCamera()
app.run()'''
