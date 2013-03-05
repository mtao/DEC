#/usr/bin/env python3

import sys
from math import sqrt
from numpy import *



class Sphere:
    def __init__(self, depth):
        self.max_depth = depth
        self.indices = []
        self.edges = {}
        self.vertices = []

    def add_edge(self, edge):
        if edge[0] > edge[1]:
            edge = (edge[1],edge[0])
        if edge in self.edges:
            return self.edges[edge]
        else:
            self.edges[edge] = len(self.vertices)
            v = (self.vertices[edge[0]] + self.vertices[edge[1]])/2.0
            v = v / linalg.norm(v)
            self.vertices.append(v)
            return len(self.vertices)-1


    def triforce(self, indices, depth):
        if depth == 0:
            self.indices.append(indices)
            return
        else:
            e01 = self.add_edge( (indices[0],indices[1]) )
            e12 = self.add_edge( (indices[1],indices[2]) )
            e02 = self.add_edge( (indices[0],indices[2]) )
            self.triforce([indices[0],e01,e02],depth-1)
            self.triforce([indices[1],e12,e01],depth-1)
            self.triforce([indices[2],e02,e12],depth-1)
            self.triforce([e01,e12,e02],depth-1)
            

    def write(self, filename):
        outfile = open(filename, "w")

        for v in self.vertices:
            outfile.write("v " + str(v[0]) + " " + str(v[1]) + " " + str(v[2]) + "\n")

        for i in self.indices:
            outfile.write("f " + str(i[0]+1) + " " + str(i[1]+1) + " " + str(i[2]+1) + "\n")


class OctahedralSphere(Sphere):
    def __init__(self,depth):
        super(OctahedralSphere,self).__init__(depth)
        self.vertices = [\
                array([1,0,0]),\
                array([0,1,0]),\
                array([-1,0,0]),\
                array([0,-1,0]),\
                array([0,0,1]),\
                array([0,0,-1]),\
                ]
        self.triforce([0,1,4],depth)
        self.triforce([1,0,5],depth)
        self.triforce([1,2,4],depth)
        self.triforce([2,1,5],depth)
        self.triforce([2,3,4],depth)
        self.triforce([3,2,5],depth)
        self.triforce([3,0,4],depth)
        self.triforce([0,3,5],depth)

class TetrahedralSphere(Sphere):
    def __init__(self,depth):
        super(TetrahedralSphere,self).__init__(depth)
        self.vertices = [\
                array([1,0,-.5**.5]),\
                array([-1,0,-.5**.5]),\
                array([0,1,.5**.5]),\
                array([0,-1,.5**.5]),\
                ]
        for i in range(len(self.vertices)):
            self.vertices[i] = self.vertices[i] / linalg.norm(self.vertices[i])

        self.triforce([0,1,2],depth)
        self.triforce([1,0,3],depth)
        self.triforce([0,2,3],depth)
        self.triforce([2,1,3],depth)

class IcosahedralSphere(Sphere):
    def __init__(self,depth):
        super(IcosahedralSphere,self).__init__(depth)
        self.vertices = [\
                array([  0          ,    -0.525731 , 0.850651]),\
                array([  0.850651   ,    0         ,  0.525731]),\
                array([  0.850651   ,    0         ,  -0.525731]),\
                array([  -0.850651  ,    0         ,  -0.525731]),\
                array([  -0.850651  ,    0         ,  0.525731]),\
                array([  -0.525731  ,    0.850651  ,  0]),\
                array([  0.525731   ,    0.850651  ,  0]),\
                array([  0.525731   ,    -0.850651 ,  0]),\
                array([  -0.525731  ,    -0.850651 ,  0]),\
                array([  0          ,    -0.525731 ,  -0.850651]),\
                array([  0          ,    0.525731  ,  -0.850651]),\
                array([  0          ,    0.525731  ,  0.850651])\
                ]
        for i in range(len(self.vertices)):
            self.vertices[i] = self.vertices[i] / linalg.norm(self.vertices[i])

        self.triforce([ 1 , 2 , 6 ],depth)
        self.triforce([ 1 , 7 , 2 ],depth)
        self.triforce([ 3 , 4 , 5 ],depth)
        self.triforce([ 4 , 3 , 8 ],depth)
        self.triforce([ 6 , 5 , 11 ],depth)
        self.triforce([ 5 , 6 , 10 ],depth)
        self.triforce([ 9 , 10 , 2 ],depth)
        self.triforce([ 10 , 9 , 3 ],depth)
        self.triforce([ 7 , 8 , 9 ],depth)
        self.triforce([ 8 , 7 , 0 ],depth)
        self.triforce([ 11 , 0 , 1 ],depth)
        self.triforce([ 0 , 11 , 4 ],depth)
        self.triforce([ 6 , 2 , 10 ],depth)
        self.triforce([ 1 , 6 , 11 ],depth)
        self.triforce([ 3 , 5 , 10 ],depth)
        self.triforce([ 5 , 4 , 11 ],depth)
        self.triforce([ 2 , 7 , 9 ],depth)
        self.triforce([ 7 , 1 , 0 ],depth)
        self.triforce([ 3 , 9 , 8 ],depth)
        self.triforce([ 4 , 8 , 0 ],depth)


for m in range(7):
    OctahedralSphere(m).write("oct"+str(m)+".obj")
    TetrahedralSphere(m).write("tet"+str(m)+".obj")
    IcosahedralSphere(m).write("icosa"+str(m)+".obj")


