import sys                                       
sys.path.insert(0, '/sdcard/Download/Internship/FWC-2/Assignments/geometry/Triangle/Altitude/codes/CoordGeo')     
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen

import subprocess
import shlex
X = np.array([19/10,-63/10])    #X is point of intersection of line AH and BC

# Given input arrays
A = np.array([-5, -4])
B = np.array([3, -3])
C = np.array([4, 0])

# H is the point of intersection of altitudes on side AB and AC from point C and B respectively...

H = np.array([113/23, -168/23])  # Reference from Question 1.3.4
 
result = int(((A - H).T) @ (B - C))    # Checking orthogonality condition...
 # printing output
if result == 0:
    print("(A - H)^T (B - C) = 0\nHence Verified...")
else:
    print("(A - H)^T (B - C)) != 0\nHence the given statement is wrong...")

#Generating all lines
x_AH = line_gen(A,H)
x_BC = line_gen(B,C)
x_BX = line_gen(B,X)

#Plotting all lines
plt.plot(x_AH[0,:],x_AH[1,:],label='$AH$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_BX[0,:],x_BX[1,:],linestyle='dotted',label='$BX$')

A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
H = H.reshape(-1,1)
X = X.reshape(-1,1)
tri_coords = np.block([[A,B,C,H,X]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','H','X']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$X-axis$')
plt.ylabel('$Y-axis$')
plt.legend(loc='upper left')
plt.grid() # minor
plt.axis('equal')
plt.show()
#if using termux
#plt.savefig('tri_sss.pdf')
plt.savefig('/sdcard/Download/Internship/FWC-2/Assignments/geometry/Triangle/Altitude/figs/AHBC_verify.png')
