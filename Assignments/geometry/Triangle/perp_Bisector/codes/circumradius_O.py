import numpy as np
from sympy import init_printing
init_printing(use_unicode=True)
from sympy.matrices import Matrix

import sys                            
sys.path.insert(0, '/sdcard/Download/Internship/FWC-2/Assignments/geometry/Triangle/perp_Bisector/codes/CoordGeo')       
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#local imports
from line.funcs import *

from triangle.funcs import *
from conics.funcs import circ_gen
A=np.array([-5,-4])
B=np.array([3,-3])
C=np.array([4,0])

O=np.array([-1.4,0.15])
X = A - O
print("value of X:",X)
radius = np.linalg.norm(X)
#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_OA = line_gen(O,A)
#Generating the circumcirclecircle
[O,r] = ccircle(A,B,C)
x_ccirc= circ_gen(O,radius)
#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_OA[0,:],x_OA[1,:],label='$OA$')

#Plotting the circumcircle
plt.plot(x_ccirc[0,:],x_ccirc[1,:],label='$circumcircle$')

A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
O = O.reshape(-1,1)
#Labeling the coordinates
tri_coords = np.block([[A,B,C,O]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','O']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$X-axis$')
plt.ylabel('$Y-axis$')
plt.legend(loc='upper left')
plt.grid()
plt.axis('equal')
plt.savefig('/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/perp_Bisector/figs/circumradius_o.png')
