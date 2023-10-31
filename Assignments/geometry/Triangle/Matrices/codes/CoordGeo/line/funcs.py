#Code by GVV Sharma
#December 7, 2019
#Revised July 15, 2020
#released under GNU GPL
#Functions related to line
import numpy as np
import mpmath as mp
#from line.params import *
from params import *


def dir_vec(A,B):
  return B-A

def norm_vec(A,B):
  return omat@dir_vec(A,B)
  #return np.matmul(omat, dir_vec(A,B))

def ang_vec(m1,m2):
    return mp.acos(float((m1.T@m2)/(np.linalg.norm(m1)*np.linalg.norm(m2))))

#Generate line points
def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB

def line_dir_pt(m,A,k1,k2):
  len = 10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(k1,k2,len)
  for i in range(len):
    temp1 = A + lam_1[i]*m
    x_AB[:,i]= temp1.T
  return x_AB


#Intersection of two lines
def line_intersect(n1,A1,n2,A2):
  N=np.block([n1,n2]).T
  p = np.zeros((2,1))
  p[0] = n1.T@A1
  p[1] = n2.T@A2
  #Intersection
  P=np.linalg.solve(N,p)
  return P


#Foot of the perpendicular
def perp_foot(n,cn,P):
  m = omat@n
  N=np.block([[n],[m]])
  p = np.zeros(2)
  p[0] = cn
  p[1] = m@P
  #Intersection
  x_0=np.linalg.solve(N,p)
  return x_0

