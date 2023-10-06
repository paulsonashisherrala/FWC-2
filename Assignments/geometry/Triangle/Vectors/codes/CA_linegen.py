import sys
sys.path.insert(0,'/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/Vectors/codes/CoordGeo') #for path to external scripts
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import math

from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen

#given points
A=np.array([-5, -4])
B=np.array([3, -3])
C=np.array([4, 0])

#getting the equation of line
omat = np.array([[0,1],[-1,0]])
m = A-C  #direction vector
n = omat@m    #normal vector
c = n@C
eqn = f"{n}x = {c}"
print("The equation of line CA is",eqn)

#plotting the line AB
x_CA = line_gen(C, A)
C = C.reshape(-1,1)
A = A.reshape(-1,1)
tri_coords = np.block([C,A])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
vert_labels=['C','A']     #for labelling points A and B
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,7), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$X-axis$')
plt.ylabel('$Y-axis$')
plt.legend(loc='best')
plt.grid()
plt.savefig('/storage/self/primary/Download/Internship/FWC-2/Assignments/geometry/Triangle/Vectors/figs/CA_line.png')
