import numpy as np
from orientedBoundedHyperPlane import orientedBoundedHyperPlane
from polygon import *
from pointcurve import *
from multipledispatch import *
from collections import *
from random import *
dimension=4
initial=[1,0,0,0]
final=[0,0,1,0]
rotation=np.identity(dimension)
x=np.array(initial)
x.shape=(1,4)
y=np.array(final)
y.shape=(1,4)
if np.matmul(x,x.transpose())[0][0]!=0:
    u=x*(1/sqrt(np.matmul(x,x.transpose())[0][0]))
else:
    u=x
u.shape=(1,4)
v=y-(np.matmul(np.matmul(u,y.transpose()),u))
if np.matmul(v,v.transpose())!=0:
    v=v*(1/sqrt(np.matmul(v,v.transpose())))
v.shape=(1,4)
uv=np.array([u.tolist()[0],v.tolist()[0]])
if np.matmul(u,u.transpose())*np.matmul(v,v.transpose())!=0:
    cosine=(np.matmul(u,v.transpose())/sqrt(np.matmul(u,u.transpose())*np.matmul(v,v.transpose())))[0][0]
else:
    cosine=1
sine=sqrt(1-cosine**2)
twod=[[cosine,-sine],[sine,cosine]]
rotation=rotation-np.matmul(u.transpose(),u)-np.matmul(v.transpose(),v)+np.matmul(uv.transpose(),np.matmul(np.array(twod),uv))
print(rotation)