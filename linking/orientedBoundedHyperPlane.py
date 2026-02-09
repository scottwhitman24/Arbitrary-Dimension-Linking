import numpy.linalg

from gensegment import *
from hyperPlane import *

class orientedBoundedHyperPlane:

    '''def __init__(self, plane, points, orientation):
        self.plane = plane
        self.points = points
        self.dimension = self.plane.dimension
        self.n = len(self.points)
        self.orientation=orientation'''

    def __init__(self, points, orientation, dimension):
        #self.plane = HyperPlane(points,dimension)
        self.points = points
        self.dimension = dimension
        self.n = len(self.points)-1
        self.orientation = orientation
    def getFaces(self):
        faces = []
        for i in range(len(self.points)):
            face = []
            for j in range(len(self.points)):
                if j != i:
                    face.append(self.points[j])
            faces.append(face)
        return faces
    def getBasis(self):
        return [difference(self.points[i+1],self.points[i]) for i in range(len(self.points)-1)]

    # def contains(self, point):
    #     faces = self.getFaces()
    #     allpos = True
    #     allneg = True
    #     for face in faces:
    #         normal = crossProduct(self.getBasis() + [difference(face[i+1], face[i]) for i in range(len(face)-1)],
    #                               self.dimension)
    #         dist = dotProduct(self.plane.getNorm(), point) / dotProduct(self.plane.getNorm(), normal)
    #         if compare(dist, 0) > 0:
    #             allneg = False
    #         elif compare(dist, 0) < 0:
    #             allpos = False
    #     return allneg | allpos

    def intersection(self, that):
        if self.n + that.n != self.dimension:
            return None
        vectors = [difference(self.points[i], self.points[0]) for i in range(self.n+1)] + [
            -difference(that.points[i], that.points[0]) for i in range(that.n+1)]
        equations = np.array([[vectors[i][j] for i in range(self.dimension)] for j in range(self.dimension)])
        selfequations = np.array([[vectors[i][j] for i in range(self.n+1)] for j in range(self.dimension)])
        try:
            sol = np.linalg.solve(equations, np.array([0 for i in range(self.dimension)]))
        except ValueError:
            return None
        return np.matmul(selfequations, sol)

    def intersectionContained(self, that):
        if self.n + that.n != self.dimension:
            return False
        vectors = ([difference(self.points[i], self.points[0]) for i in range(1,self.n+1)]
                   + [difference(that.points[i], that.points[0]) for i in range(1,that.n+1)])
        equations = np.array([[vectors[i][j] for i in range(self.dimension)] for j in range(self.dimension)])
        #selfequations = np.array([[vectors[i][j] for i in range(self.n)] for j in range(self.dimension)])
        try:
            #u,s,vh = numpy.linalg.svd(equations)
            #print(vh)
            #sol=numpy.compress(s<=0.0001,vh,axis=0)
            #print(s)
            sol = np.linalg.solve(equations, np.array(difference(that.points[0],self.points[0])))
            #print(sol)
        except ValueError:
            #print("Error")
            return False
        if sol.size==0:
            return 0
        solution=sol.tolist()
        if compare(min(solution[:self.n]),0)>=0:
            if compare(sum(solution[:self.n]),1)<=0:
                return True
        return False

    def intersectsWith(self, that):
        return self.intersectionContained(that) and that.intersectionContained(self)

    def signedLinksWith(self, that):
        if self.intersectsWith(that):
            return compare((np.linalg.det(np.array(self.getBasis()+that.getBasis())))*self.orientation*that.orientation,0)
        else:
            return 0


    # @dispatch(genSegment)
    # def alignmentWithNormal(self, edge):
    #     product = dotProduct(edge.vector(), self.plane.getNorm())
    #     # print("alignment: "+str(product))
    #     if compare(product, 0) == 0:
    #         return 0
    #     elif compare(product, 0) < 0:
    #         return -1
    #     elif compare(product, 0) > 0:
    #         return 1
    #
    #     return 2
    #
    # @dispatch(list)
    # def alignmentWithNormal(self, edge):
    #     product = dotProduct(edge, self.plane.getNorm())
    #     if compare(product, 0) == 0:
    #         return 0
    #     elif compare(product, 0) < 0:
    #         return -1
    #     elif compare(product, 0) > 0:
    #         return 1
    #
    #     return 2
    #
    # def get_intersection(self, that):
    #     for edge in that.edges:
    #         if self.intersectedBy(edge):
    #             return intersection(edge, self.plane)


# def compare(norm, face, point, dimension):
#     if dimension == 1:
#         return alignment(difference(point, face[0]), norm)
#     normal = crossProduct(norm + [difference(face[i], face[0]) for i in range(1, len(face))], dimension)
#     t = dotProduct(norm, point) / dotProduct(norm, normal)
#     val = compare(dist, 0)
#     if val > 0:
#         return 1
#     elif val < 0:
#         return -1
#     else:
#         return 0


def intersects(segment, bounded):
    if (compare(dotProduct(segment.vector(), bounded.plane.getNorm()), 0) != 0 or
            compare(dotProduct(segment.getStart(), bounded.plane.getNorm()), bounded.plane.getK()) == 0):
        if compare(dotProduct(segment.vector(), bounded.plane.getNorm()), 0) == 0:
            return False
        t = (bounded.plane.getK() - dotProduct(bounded.plane.getNorm(), segment.getStart())) / (
            dotProduct(segment.vector(), bounded.plane.getNorm()))
        if compare(t, 1) > 0:
            return False
        if compare(t, 0) < 0:
            return False

        return True

    return False


def intersection(bounded, segment):
    if compare(dotProduct(segment.vector(), bounded.plane.getNorm()), 0) == 0:
        if compare(dotProduct(segment.getStart(), bounded.plane.getNorm()), bounded.plane.getK()) == 0:
            return segment.getStart()

        else:
            print(",")
            return None

    t = (bounded.plane.getK() - dotProduct(bounded.plane.getNorm(), segment.getStart())) / (
        dotProduct(segment.vector(), bounded.plane.getNorm()))
    if compare(t, 1) > 0:
        return None

    if compare(t, 0) < 0:
        return None
    return sum(segment.getStart(), product(segment.vector(), t))


def alignment(v, u):
    return compare(dotProduct(v, u), 0)


def getFaces(points):
    faces = []
    for i in range(len(points)):
        face = []
        for j in range(len(points)):
            if j != i:
                face.append(points[j])
        faces.append(face)
    return faces
