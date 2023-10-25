#Code by GVV Sharma
#October 2, 2023
#released under GNU GPL
#Angle Bisectors of a triangle
#Incentre and Incircle

import sys  #for path to external scripts
sys.path.insert(0, '/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/Angle_bisector/codes/CoordGeo')#path to my scripts
import numpy as np
import scipy.linalg as SA
import mpmath as mp
import numpy.linalg as LA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import *


#if using termux
import subprocess
import shlex
#end if

A=np.array([-5,-4])
B=np.array([3,-3])
C=np.array([4,0])

#values of m,n,p
a=np.linalg.norm(C-B)
b=np.linalg.norm(C-A)
c=np.linalg.norm(A-B)

#Orthogonal matrix
omat = np.array([[0,1],[-1,0]]) 

#finding Incentre 
t1 = norm_vec(B,C) 
n1 = t1/np.linalg.norm(t1) #unit normal vector
t2 = norm_vec(C,A)
n2 = t2/np.linalg.norm(t2)
t3 = norm_vec(A,B)
n3 = t3/np.linalg.norm(t3)

I=line_intersect(n1-n3,B,n1-n2,C)  #Incentre
#Incentre
print("I = ",I)

#finding k for E_3 and F_3
k1=((I-A)@(A-B))/((A-B)@(A-B))
k2=((I-A)@(A-C))/((A-C)@(A-C))
k3=((I-B)@(C-B))/((C-B)@(C-B))

#finding E_3 and F_3
F3=A+(k1*(A-B))
E3=A+(k2*(A-C))
D3=B+(k3*(C-B))
print("k1 = ",k1)
print("k2 = ",k2)
print("k3 = ",k3)
print("E3 = ",E3)
print("F3 = ",F3)
print("D3 = ",D3)

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
I = I.reshape(-1,1)
E3 = E3.reshape(-1,1)
F3 = F3.reshape(-1,1)
D3 = D3.reshape(-1,1)

#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_ID_3 = line_gen(I,D3)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_ID_3[0,:],x_ID_3[1,:],label='$ID_3$')


tri_coords = np.block([[A,B,C,I,E3,F3,D3]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','I','E3','F3','D3']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$X-axis$')
plt.ylabel('$Y-axis$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig('/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/Angle_bisector/figs/Incentre_I.png')
