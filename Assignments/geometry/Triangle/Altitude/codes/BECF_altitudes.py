import sys
sys.path.insert(0, '/sdcard/Download/Internship/FWC-2/Assignments/geometry/Triangle/Altitude/codes/CoordGeo') #path to my script    s
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

A = np.array([-5, -4])
B = np.array([3, -3])
C = np.array([4, 0])
omat = np.array([[0, 1], [-1, 0]])

D =  alt_foot(A,B,C)
E =  alt_foot(B,C,A)
F =  alt_foot(C,A,B)

# Print altitude foot points
print("Altitude foot point for A:", D)
print("Altitude foot point for B:", E)
print("Altitude foot point for C:", F)

BE_norm = norm_vec(E, B)
CF_norm = norm_vec(F, C)

print(BE_norm)
print(CF_norm)

print(f"{BE_norm}[x-B]=0" )
print(f"{CF_norm}[x-C]=0" )

#Normal vectors of AD and BE
n1 = B-C
n2 = C-A

#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_AD = line_gen(D,A)
x_BE = line_gen(B,E)
x_CF = line_gen(C,F)
x_BD = line_gen(B,D)
x_BF = line_gen(B,F)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')
plt.plot(x_BD[0,:],x_BD[1,:],label='$BD$',linestyle='dotted')
plt.plot(x_BF[0,:],x_BF[1,:],label='$BF$',linestyle='dotted')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
tri_coords = np.block([[A,B,C,D,E,F]])
#tri_coords = np.vstack((A,B,C,alt_foot(A,B,C),alt_foot(B,A,C),alt_foot(C,A,B),H)).T
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','E','F']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,9), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or c
plt.xlabel('$X-axis$')
plt.ylabel('$Y-axis$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig('/sdcard/Download/Internship/FWC-2/Assignments/geometry/Triangle/Altitude/figs/BECF_altitudes.png')
#subprocess.run(shlex.split("termux-open ./figs/tri_sss.pdf"))
#else
# image = mpimg.imread('tri_sss.png')
# plt.imshow(image)
plt.show()
