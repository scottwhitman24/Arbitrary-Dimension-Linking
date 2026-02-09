import numpy as np
class simplicialObject:
    def __init__(self,interior,boundary,simplexdimension,ambientdimension):
        self.interior=interior
        self.boundary=boundary
        self.simplexdimension = simplexdimension
        self.ambientdimension = ambientdimension

    def __init__(self,simplex,simplexdimension,ambientdimension):
        self.simplexdimension = simplexdimension
        self.ambientdimension = ambientdimension
        self.interior=[[simplex,1]]
        self.boundary=[[np.delete(self.interior[0][0],j,axis=0), j, 1] for j in range(self.simplexdimension + 1)]

    def refine(self):
        face=self.boundary[0]
        newvertex=face[0].sum(axis=0)*1/(self.simplexdimension)
        newvertex=newvertex/np.linalg.norm(newvertex)
        self.interior.append([np.insert(face[0],-1,newvertex,axis=0),-face[2]])
        self.boundary.pop(0)
        for i in range(self.simplexdimension):
            newface=np.insert(face[0],i+1,newvertex,axis=0)
            self.boundary.append([newface,face[1],face[2],self.interior[-1][1]])
    def __str__(self):
        return "Simplicial Object of dimension "+str(self.simplexdimension)+" in ambient dimension "+str(self.ambientdimension)+"\nInterior: "+str(self.interior)+"\nBoundary: "+str(self.boundary)


object=simplicialObject(np.eye(4,6),3,6)
print(object)
for i in range(30):
    object.refine()
print(object)
print(object.boundary[50])