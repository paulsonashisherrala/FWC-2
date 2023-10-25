import sys
sys.path.insert(0, '/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/perp_Bisector/codes/CoordGeo')  #path to my scripts
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen

#if using termux
import subprocess
import shlex 
#end if
 
# enter vectors A,B & C
A=np.array([-5,-4])
B=np.array([3,-3])
C=np.array([4,0])

# direction vector along line joining A & B
AB = dir_vec(A,B)
# direction vector along line joining A & C
AC = dir_vec(A,C)

# midpoint of A & B is F
F = (A+B)/2
# midpoint of A & C is E
E = (A+C)/2

# O is the point of intersection of perpendicular bisectors of AB and AC
O = line_intersect(AB,F,AC,E)
print("O",O)

G=(C+B)/2
#Generating all lines 
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
#x_OE = line_gen(O,E)
#x_OF = line_gen(O,F)
x_OG=line_gen(O,G)

#plotting all lines 
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_OG[0,:],x_OG[1,:],label='$OG$')

A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
O = O.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
G= G.reshape(-1,1)
tri_coords = np.block([[A,B,C,O,E,F,G]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','O','E','F','G']
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
plt.savefig('/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/perp_Bisector/figs/O_satisfy.png')
