import sys                               #for path to external scripts
sys.path.insert(0, './CoordGeo')        #path to my scripts
import numpy as np
import numpy.linalg as LA
import scipy.linalg as SA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen

#if using termux
import subprocess
import shlex
#end if

B = np.array([0,0]).reshape(-1,1)#B(0,0)Origin
C = np.array([6,0]).reshape(-1,1)#C(6,0) on same axis 
l = LA.norm(B-C) #length of BC
D = np.array([0,l]).reshape(-1,1) #(X=0,Y=l) taking 'l' as Y-coordinate
A = np.array([l,l]).reshape(-1,1) #(X=l,Y=l) taking 'l'as x,y coordinate
M = midpoint(A,B) #midpoint of (A,B) 
#print(" A = {} \n B ={} \n C= {} \n D= {} \n M= {} \n l={} ".format(A,B,C,D,M,l))

print("*************************************** MATH COMPUTING ****************************************** \n")
#print("1.To prove △ AMC ≅ △ BMD")
DM = round(LA.norm(D-M),5)    #length of AD
CM = round(LA.norm(C-M),5)   #length of AC  
BM = round(LA.norm(B-M),5)    #length of BM
AM = round(LA.norm(A-M),5)    #length of AM
#print("DM = {} \n CM = {} \n BM = {} \n AM = {}".format(DM,CM,BM,AM))

#finding angle CMA
dot_C = ((C - M).T) @ (A - M)
norm_C = np.linalg.norm(C - M) * np.linalg.norm(A - M)
cos_theta_C = dot_C / norm_C
angle_CMA =np.round(np.degrees(np.arccos(cos_theta_C)),2)
#print("angle CMA = ",angle_CMA)

#finding angle BMD
dot_B = ((B - M).T) @ (D - M)
norm_B = np.linalg.norm(B - M) * np.linalg.norm(D - M)
cos_theta_B = dot_B / norm_B
angle_BMD =np.round(np.degrees(np.arccos(cos_theta_B)),2)
#print("angle BMD = ",angle_BMD)


if((DM==CM) and (BM==AM) and (angle_CMA==angle_BMD)):
    print("△ AMC ≅ △ BMD (congruent By SAS Congruency)")
else:
    print("△ AMC ≇ △ BMD (is NOt congruent) ")


#print("\n 2. for proving ∠ B is right angle ")
#finding the angle DBC 
dot_D = ((D - B).T) @ (C - B)
norm_D = np.linalg.norm(D - B) * np.linalg.norm(C - B)
cos_theta_D = dot_D / norm_D
angle_DBC =np.round(np.degrees(np.arccos(cos_theta_D)),2)
#print("angle_DBC = ",angle_DBC)

if(angle_DBC == 90):
    print("△ DBC IS right angled At B ")
else:
    print("△ DBC IS NOT right angled At B")

#print("\n 3. ∆ DBC ≅ ∆ ACB ")
#finding angle BCA
dot_A = ((B - C).T) @ (A - C)
norm_A = np.linalg.norm(B - C) * np.linalg.norm(A - C)
cos_theta_A = dot_A / norm_A
angle_BCA =np.round(np.degrees(np.arccos(cos_theta_A)),2)
#print("angle_BCA = ",angle_BCA)

CB = LA.norm(B-C)   #length of AD  
BC = LA.norm(C-B)   #length of AC  
DB = LA.norm(B-D)   #length of BM
AC = LA.norm(C-A)   #length of AM
#print("CB = {} \n BC = {} \n DB = {} \n AC ={}".format(CB,BC,DB,AC))


if((CB==BC) and (DB==AC) and (angle_DBC==angle_BCA)):
    print("△ DBC ≅ △ ACB (congruent By SAS Congruency)")
else:
    print("△ DBC ≇ △ ACB (is NOt congruent) ")

#print("4.CM = 1/2  AB")
AB = round(LA.norm(B-A),5) #length of AB
if((CM == AB/2) or (2*CM == AB)):
    print("CM = AB/2")
else:
    print("CM ≠ AB/2")

print("*************************************** MATH COMPUTING ****************************************** \n")

#plotting the graph
x_AB = line_gen(A,B)
x_BD = line_gen(B,D)
x_BC = line_gen(B,C)
x_AC = line_gen(A,C)
x_CD = line_gen(C,D)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BD[0,:],x_BD[1,:],label='$BD$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_AC[0,:],x_AC[1,:],label='$AC$')
plt.plot(x_CD[0,:],x_CD[1,:],label='$CD$') 

#Labeling the coordinates
tri_coords = np.block([[A,B,C,D,M]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','M']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y) 
                 ha='center') # horizontal alignmnt can be left, right or center
plt.xlabel('$X-axis$')
plt.ylabel('$Y-axis$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig('/home/ashishroy007/Desktop/FWC-2/Assignments/Math_Computing/figs/python_plot.png')
plt.show()

