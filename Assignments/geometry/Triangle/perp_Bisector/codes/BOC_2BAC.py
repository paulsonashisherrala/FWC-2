import numpy as np
A = np.array([-5,-4])
B = np.array([3, -3])
C = np.array([4, 0])
# O is a point of intersection of perpendicular bisectors of AB and AC
O = np.array([-1.4565 ,0.1521])

#To find angle AOC
dot_pt_O = (B - O) @ ((C - O).T)
norm_pt_O = np.linalg.norm(B - O) * np.linalg.norm(C - O)
cos_theta_O = dot_pt_O / norm_pt_O
angle_BOC = round(np.degrees(np.arccos(cos_theta_O)),2)  #Round is used to round of number till 5 decimal places
print("angle BOC = ",angle_BOC)

#To find angle BAC
dot_pt_B = (B - A) @ ((C - A).T)
norm_pt_B = np.linalg.norm(B - A) * np.linalg.norm(C - A)
cos_theta_B = dot_pt_B / norm_pt_B
angle_BAC = round(np.degrees(np.arccos(cos_theta_B)),2)  #Round is used to round of number till 5 decimal places
print("angle BAC = ",angle_BAC)
#To check whether the answer is correct
if angle_BOC == 2 * angle_BAC:
  print("\nangle AOC = 2 times angle ABC\nHence the give statement is correct")
else:
  print("\nangle AOC ≠ 2 times angle ABC\nHence the given statement is wrong")
import sys                                        
sys.path.insert(0, '/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/perp_Bisector/codes/CoordGeo')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen
import subprocess
import shlex
#end if

O = np.array([-67/46, 7/46])
A = np.array([-5, -4])
B = np.array([3, -3])
C = np.array([4, 0])

#Generating all lines
x_AB = line_gen(A,B)
x_AC = line_gen(A,C)
x_OB = line_gen(O,B)
x_OC = line_gen(O,C)
x_OA = line_gen(O,A)
x_BC = line_gen(B,C)

#Generating the circumcircle
[O,R] = ccircle(A,B,C)
x_circ= circ_gen(O,R)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_OC[0,:],x_OC[1,:],label='$OC$')
plt.plot(x_OA[0,:],x_OA[1,:],label='$OA$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_AC[0,:],x_AC[1,:],linestyle='dashed',label='$AC$')
plt.plot(x_OB[0,:],x_OB[1,:],linestyle='dashed',label='$OB$')



#Plotting the circumcircle
plt.plot(x_circ[0,:],x_circ[1,:],label='$circumcircle$')


#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
O = O.reshape(-1,1)
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
plt.savefig('/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/perp_Bisector/figs/BOC_2BAC.png')

