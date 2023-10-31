#Code by GVV Sharma
#Revised October 1, 2023
#Revised October 2, 2023
#Revised October 20, 2023
#released under GNU GPL
#Matrix Algebra

import sys                                          #for path to external scripts
sys.path.insert(0, '/sdcard/Download/Internship/FWC-2/Assignments/geometry/Triangle/Matrices/codes/CoordGeo')        #path to my scripts
import numpy as np
import math
import numpy.linalg as LA
import scipy.linalg as SA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen
from plotting.funcs import label_pts


#if using termux
import subprocess
import shlex
#end if

#-----------------Vectors-------------------------------

#Input parameters from excel file
df= pd.read_excel('/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/Matrices/tables/vertices.xlsx')
print(df)

#Triangle Vertices
G_v= df.to_numpy()[:,:]
print(G_v)

#Direction vector circulant matrix
C_m= SA.circulant([1,0,-1]).T


#Direction vector Matrix
G_dir = G_v@C_m
print(G_dir)

#Normal vector matrix
G_n = R_o@G_v@C_m
print(G_n)

#Find the line constants
cmat = np.diag(G_n.T@G_v).reshape(-1,1)

#line matrix
linmat = np.block([G_n.T,cmat])
#print(linmat)

#sides vector
C_dis= SA.circulant([0,1,0]).T
dis = C_dis@np.linalg.norm(G_dir, axis=0).reshape(-1,1)
#print(a,b,c)
'''
#Finding the angles of the triangle
dmat = np.diag(1/d)
G_dnorm = G_dir@dmat
G_dgram = G_dnorm.T@G_dnorm
#print(np.degrees(np.arccos(G_dgram)))
'''


#-----------------Vectors Ends-------------------------------

#-----------------Medians-------------------------------
#Median circulant matrix
C_mid = SA.circulant([0,1,1]).T
#print(C_mid)

#Mid point matrix
G_mid = 0.5*G_v@C_mid
#print(G_mid)

#Median direction circulant matrix
C_mid_dir = SA.circulant([1,-0.5,-0.5])

#Median direction matrix
G_med_dir = G_v @ C_mid_dir
#print(G_med_dir)

#Normal vector matrix
G_n_med = R_o@G_med_dir

#Find the line constants
cmat_med = np.diag(G_n_med.T@G_v).reshape(-1,1)

#median  matrix
linmat_med = np.block([G_n_med.T,cmat_med])
#print(linmat_med)

#Find the centroid
G  = LA.lstsq(G_n_med.T,cmat_med,rcond=None)

#-----------------Median Ends-------------------------------

#-----------------Altitude-------------------------------

#Circulant matrix
C_alt= SA.circulant([0,-1,1]).T

#Normal Matrix
G_dir_alt = G_v@C_alt

#Find the line constants
cmat_alt = np.diag(G_dir_alt.T@G_v).reshape(-1,1)

#altitude matrix
linmat_alt= np.block([G_dir_alt.T,cmat_alt])

#Find the orthocentre
H  = LA.lstsq(G_dir_alt.T,cmat_alt,rcond=None)

#-----------------Altitude Ends-------------------------------

#-----------------Perpendicular Bisector-------------------------------
#Find the line constants
cmat_perp_bis= np.diag(G_dir_alt.T@G_mid).reshape(-1,1)
#print( np.block([G_dir_alt.T,cmat_perp_bis]))

#Find the Circumcentre
O = LA.lstsq(G_dir_alt.T,cmat_perp_bis,rcond=None)
#-----------------Perpendicular Bisector Ends -------------------------------

#-----------------Angle Bisector-------------------------------
#Incircle circulant matrix
C_in = SA.circulant([1,1,0]).T
#m,n,p
secvec = LA.inv(C_in)@dis
m,n,p  = secvec 
a,b,c= dis
#bisector matrix
cont_mat =np.array([1/a*np.block([0,n,m]), 1/b*np.block([n,0,p]),1/c*np.block([m,p,0])]).T
G_ang_bis = G_v@cont_mat 

#Normal Matrix
G_dir_alt_ang = G_ang_bis@C_alt
#Mid point matrix
G_mid_ang = 0.5*G_ang_bis@C_mid
#Find the line constants
cmat_ang_bis= np.diag(G_dir_alt_ang.T@G_mid_ang).reshape(-1,1)
print( np.block([G_dir_alt.T,cmat_perp_bis]))

#Find the incentre
I = LA.lstsq(G_dir_alt_ang.T,cmat_ang_bis,rcond=None)
print("I",I[0])
print("G_ang_bis",G_ang_bis)
#-----------------End Angle Bisector-------------------------------
#print(G)
#print(H)
#print(O)
#-----------------Plotting-------------------------------
# Generating all sides
x = np.zeros((3, 2, 10))
for i in range(3):
    x[i, :, :] = line_gen(G_v[:, i], G_v[:, (i + 1) % 3])

# Plotting the sides and labeling them
lines = ['AB', 'BC', 'CA']
for i in range(3):
    plt.plot(x[i, 0, :], x[i, 1, :], label=lines[i])
    
# Labeling the coordinates
plt.scatter(G_v[0, :], G_v[1, :])
plt.scatter(G_ang_bis[0, :], G_ang_bis[1, :])
vert_labels = ['A', 'B', 'C']
label_pts(G_v, vert_labels)


#labeling the angle bisectors
angle_bisector_labels = ['$D_3$', '$E_3$', '$F_3$']
label_pts(G_ang_bis, angle_bisector_labels)

## Generate the equation of the incircle
#ID3 = 1.091423664# Radius of the incircle
#I_center = I[0]  # Center of the incircle
#incircle_eq = circ_gen(I_center, ID3)
#
## Create a set of points on the incircle for plotting
#theta = np.linspace(0, 2 * np.pi, 100)
#x_incircle = I_center[0] + ID3 * np.cos(theta)
#y_incircle = I_center[1] + ID3 * np.sin(theta)
#
## Plot the incenter (I) and the incircle
#plt.scatter(I[0][0], I[0][1], color='red', marker='o'); plt.text(I[0][0], I[0][1], 'I', ha='left', va='bottom')
#plt.plot(x_incircle, y_incircle, label='Incircle')

plt.xlabel('$X-axis$')
plt.ylabel('$Y-axis$') 
plt.legend(loc='upper left')
plt.grid()
plt.axis('equal')
plt.show()
#if using termux
plt.savefig('/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/Matrices/figs/D3E3F3.png')
#subprocess.run(shlex.split("termux-open ./figs/tri_sss.pdf"))
#else
plt.show()
