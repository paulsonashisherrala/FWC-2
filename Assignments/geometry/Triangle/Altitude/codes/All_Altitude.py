import sys
sys.path.insert(0, '/sdcard/Download/Internship/FWC-2/Assignments/geometry/Triangle/Altitude/codes/CoordGeo') #path to my scripts
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

#vertices of triangle
A = np.array([-5,-4])
B = np.array([3,-3])
C = np.array([4,0])

#orthogonal matrice
omat = np.array([[0, 1], [-1, 0]])

#altitude foot
D = alt_foot(A,B,C)
E = alt_foot(B,A,C)
F = alt_foot(C,A,B)
print(D,E,F)

#norm 
BE_norm= norm_vec(E,B)
CF_norm= norm_vec(F,C)
AD_norm= norm_vec(D,A)

#direction vecrors
n1=B-C
n2=C-A
n3=B-A

#Finding orthocentre
H = line_intersect(norm_vec(B,E),E,norm_vec(C,F),F)
print(H)

#Generating all lines
x_AB = line_gen(A,B)  
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_AD = line_gen(A,alt_foot(A,B,C))
x_AE = line_gen(A,alt_foot(B,A,C))
x_BE = line_gen(B,alt_foot(B,A,C))
x_CF = line_gen(C,alt_foot(C,A,B))
x_AF = line_gen(A,alt_foot(C,A,B))
x_BD = line_gen(B,D)
x_DH = line_gen(D,H)
x_BH = line_gen(B,H)
x_FH = line_gen(F,H)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE_1$')
#plt.plot(x_AE[0,:],x_AE[1,:],linestyle = 'dashed',label='$AE_1$')
plt.plot(x_CF[0,:],x_CF[1,:],label='$CF_1$')
plt.plot(x_AF[0,:],x_AF[1,:],linestyle = 'dashed',label='$AF_1$')
plt.plot(x_DH[0,:],x_DH[1,:],linestyle='dashed',label='$DH$')
plt.plot(x_BH[0,:],x_BH[1,:],linestyle='dashed',label='$BH$')
plt.plot(x_FH[0,:],x_FH[1,:],linestyle ='dashed',label='$FH$')
plt.plot(x_BD[0,:],x_BD[1,:],linestyle='dashed',label='$BD$')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
H = H.reshape(-1,1)
tri_coords = np.block([[A,B,C,D,E,F,H]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','$D_1$','$E_1$','$F_1$','H']
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
plt.savefig('/sdcard/Download/Internship/FWC-2/Assignments/geometry/Triangle/Altitude/figs/All_Altitude.pdf')
plt.show()
