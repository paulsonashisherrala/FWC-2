#Code by GVV Sharma
#December 7, 2019
#Revised November 12, 2023
#Revised July 15, 2020
#released under GNU GPL
#Functions related to conics

import numpy as np
import numpy.linalg as LA
from params import *

#Generating points on a circle
def circ_gen(O,r):
	len = 50
	theta = np.linspace(0,2*np.pi,len)
	x_circ = np.zeros((2,len))
	x_circ[0,:] = r*np.cos(theta)
	x_circ[1,:] = r*np.sin(theta)
	x_circ = (x_circ + O)
	return x_circ

#Generating points on an ellipse
def ellipse_gen(a,b):
	len = 50
	theta = np.linspace(0,2*np.pi,len)
	x_ellipse = np.zeros((2,len))
	x_ellipse[0,:] = a*np.cos(theta)
	x_ellipse[1,:] = b*np.sin(theta)
	return x_ellipse

#Generating points on a parabola
def parab_gen(y,a):
	x = y**2/a
	return x

#Generating points on a standard hyperbola 
def hyper_gen(y):
	x = np.sqrt(1+y**2)
	return x

#Points of contact for a pair of tangents
def contact(V,u,f,h):
    #intermediate
    gh = h.T@V@h+2*u.T@h+f 
    
    #matrix of tangents
    sigmat = (V@h+u)@(V@h+u).T-gh*V
    
    
    #Spectral decomposition
    D, P = LA.eig(sigmat)
    
    u1 = np.array(([np.sqrt(np.abs(D[1])),np.sqrt(np.abs(D[0]))]))
    u2 = np.array(([np.sqrt(np.abs(D[1])),-np.sqrt(np.abs(D[0]))]))
    
    u1 = u1.reshape(-1,1)
    u2 = u2.reshape(-1,1)
    
    #direction vectors
    m1 = P@u1
    m2 = P@u2
    #print(m1,m2)
    # Converting 1D array to a 2D numpy array of incompatible shape will cause error
    m1= np.reshape(m1, (2, 1))
    m2= np.reshape(m2, (2, 1))
    mu1 = -(m1.T@(V@h+u))/(m1.T@V@m1)
    mu2 = -(m2.T@(V@h+u))/(m2.T@V@m2)
    #print(mu1,mu2)
    x1 = h + mu1*m1
    x2 = h + mu2*m2
    return(x1,x2)
    
