//Code by G V V Sharma
//October 27, 2023
//matrix operations using arrays
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libs/matfun.h"


int  main()
{
FILE *fp; //file pointer
double **A,**B,**C;//declare vertices
double **D,**E,**F;//mid_points
double **m1,**m2,**m3,**dir_vec,**m4,**m5,**m6,**mdir_vec;//direction vectors
double **n1,**n2,**n3,**norm_vec,**n4,**n5,**n6,**mnorm_vec,**anorm_vec;//normal vectors
double **temp;//temporary array
double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;//constant vectors
double **l_vec,**t_vec,**c_vec,**m_vec,**ml_vec,**mt_vec,**mc_vec,**mm_vec,**ac_vec,
       **at_vec,**alc_vec,**am_vec,**pc_vec,**plc_vec,**pm_vec;
double a,b,c,pv,mv,nv;//triangle sides
double **cx1,**cy1;
double **D3,**E3,**F3;
double **s_vec,**ss_vec;
double **midpt_vec,**centroid,**o_center,**circum_c;
int m =2, n=3, i, j;
double **mat =createMat(m,n);//vertices matrix 
mat = loadMat("vertices.dat",m, n);//loading matrix from file
//temp= loadMat("circ.dat",3, 1);//loading matrix from file
temp= loadMat("matex.dat",2, 2);//loading matrix from file

//Extracting vertices
A = Matcol(mat,m,0);
B = Matcol(mat,m,1);
C = Matcol(mat,m,2);
printf("Vertices Matrix = \n");
printMat(mat,2,3);
printf("\n");

//***************************************** VECTORS *****************************************
printf("**************************** Vectors ********************************** \n");
//Direction vectors
m1 = Matsub(A,B,m,1);
m2 = Matsub(B,C,m,1);
m3 = Matsub(C,A,m,1);
dir_vec=matrix_merge(m1,m2,m3,2,3);
printf("Direction Matrix = \n");
printMat(dir_vec,2,3);
printf("\n");

//Line Parameters
//Normal vectors
n1 = normVec(m1);
n2 = normVec(m2);
n3 = normVec(m3);
norm_vec=matrix_merge(n1,n2,n3,2,3);
printf("Normal Matrix = \n");
printMat(norm_vec,2,3);
printf("\n");

//Line constants
c1 = Matdot(n1,A,2);
c2 = Matdot(n2,B,2);
c3 = Matdot(n3,C,2);
l_vec = matrix(c1,c2,c3);
printf("Constant Matrix = \n");
//printf("%lf %lf %lf \n",c1,c2,c3);
printMat(l_vec,1,3);
printf("\n");

//Sides lengths
a = Matnorm(m2,m);
b = Matnorm(m3,m);
c = Matnorm(m1,m);
s_vec = matrix(a,b,c);
printf("Distance matrix = \n");
printMat(s_vec,1,3);
//printf("\nDistance:\n %lf , %lf ,%lf \n",a,b,c);
printf("\n");

//Line Matrix
t_vec =transposeMat(norm_vec,2,3);
//printMat(t_vec,3,2);
c_vec =transposeMat(l_vec,1,3);
//printMat(c_vec,3,1);
m_vec =matrix_2Merge(t_vec,3, 2,c_vec, 3, 1);
printf("Line matrix = \n");
printMat(m_vec, 3, 3);
printf("\n");

//Rotation matrix
//printMat(rotMat(M_PI/2),2,2);

//Circulant matrix
//circulantMat(temp, 3);
//printMat(temp,3,3);

//Matrix inversion
//printMat(Matinv(temp, 2),2,2);

//***************************************END OF VECTORS*************************************

//******************************************MEDIAN******************************************
printf("**************************** Medians ********************************** \n");

//Midpoint matrix
D= Matsec(C,B,2,1);
//printMat(D,2,1);

E= Matsec(A,C,2,1);
//printMat(E,2,1);

F= Matsec(A,B,2,1);
//printMat(F,2,1);

midpt_vec=matrix_merge(D,E,F,2,3);
printf("Midpoint Matrix = \n");
printMat(midpt_vec,2,3);
printf("\n");

//Direction vectors
m4 = Matsub(A,D,m,1);
m5 = Matsub(B,E,m,1);
m6 = Matsub(C,F,m,1);
mdir_vec=matrix_merge(m4,m5,m6,2,3);
printf("Median Direction Matrix = \n");
printMat(mdir_vec,2,3);
printf("\n");


//Normal vectors
n4 = normVec(m4);
n5 = normVec(m5);
n6 = normVec(m6);
mnorm_vec=matrix_merge(n4,n5,n6,2,3);
printf("Median Normal Matrix = \n");
printMat(mnorm_vec,2,3);
printf("\n");

//Line constants
c4 = Matdot(n4,D,2);
c5 = Matdot(n5,E,2);
c6 = Matdot(n6,F,2);
ml_vec = matrix(c4,c5,c6);
printf("Median Constant Matrix = \n");
//printf("%lf %lf %lf \n",c4,c5,c6);
printMat(ml_vec,1,3);
printf("\n");

//centroid
cx1=Matcol(transposeMat(mnorm_vec,2,3),m,0);
cy1=Matcol(transposeMat(mnorm_vec,2,3),m,1);

//Line Matrix
mt_vec =transposeMat(mnorm_vec,2,3);
//printMat(t_vec,3,2);
mc_vec =transposeMat(ml_vec,1,3);
//printMat(c_vec,3,1);
mm_vec =matrix_2Merge(mt_vec,3, 2,mc_vec, 3, 1);
printf("Median Line matrix = \n");
printMat(mm_vec, 3, 3);
printf("\n");

centroid=line_intersect(mm_vec,3,3);
printf("Centroid = \n");
printMat(centroid,1,2);
printf("\n");
//***************************************END OF MEDIANS*************************************


//***************************************ALTITUDE******************************************
printf("**************************** Altitude ********************************** \n");

//Normal Matrix
printf("Altitude Normal Matrix = \n");
anorm_vec= matrix_merge(m2,m3,m1,2,3);
printMat(anorm_vec,2,3);
printf("\n");

//constant Matrix
c7 = Matdot(m2,A,2);
c8 = Matdot(m3,B,2);
c9 = Matdot(m1,C,2);
ac_vec = matrix(c7,c8,c9);
printf("Altitude Constant Matrix = \n");
//printf("%lf %lf %lf \n",c7,c8,c9);
printMat(ac_vec,1,3);
printf("\n");

at_vec =transposeMat(anorm_vec,2,3);
//printMat(at_vec,3,2);
alc_vec =transposeMat(ac_vec,1,3);
//printMat(alc_vec,3,1);
//merging two into one
am_vec =matrix_2Merge(at_vec,3, 2,alc_vec, 3, 1);
printf("Altitude Line matrix = \n");
printMat(am_vec, 3, 3);
printf("\n");

o_center=line_intersect(am_vec,3,3);
printf("Orthocenter = \n");
printMat(o_center,1,2);
printf("\n");

//***************************************END OF ALTITUDE************************************

//***************************************PERPENDICULAR BISECTOR*****************************
printf("**************************** Perpendicular Bisector ********************************** \n");
c10 = Matdot(m2,D,2);
c11 = Matdot(m3,E,2);
c12 = Matdot(m1,F,2);

pc_vec = matrix(c10,c11,c12);
printf("Perpendicular Constant Matrix = \n");
//printf("%lf %lf %lf \n",c10,c11,c12);
printMat(pc_vec,1,3);
printf("\n");

//merging two matrices into oneprintf("\nIncentre = \n");
plc_vec = transposeMat(pc_vec,1,3);
pm_vec = matrix_2Merge(at_vec,3,2,plc_vec,3,1);
printf("Perpendicular Line matrix = \n");
printMat(pm_vec, 3, 3);
printf("\n");

circum_c=line_intersect(pm_vec,3,3);
printf("Circumcentre = \n");
printMat(circum_c,1,2);
printf("\n");

//********************************END OF PERPENDICULAR BISECTOR*****************************

//*************************************** ANGULAR BISECTOR *****************************
printf("**************************** Angular Bisector ********************************** \n");
//Direction vectors
//m1 = Matsub(A,B,m,1);
//m2 = Matsub(B,C,m,1);
//m3 = Matsub(C,A,m,1);

//Sides lengths
//a = Matnorm(m2,m);
//b = Matnorm(m3,m);
//c = Matnorm(m1,m);

//calculations for p,m,n
pv = ((c+b)-a)/2;
mv = ((a+c)-b)/2;
nv = ((a+b)-c)/2;
ss_vec = matrix(pv,mv,nv); 
printf("p,m,n Values = \n");
printMat(ss_vec,1,3);
printf("\n");

//double **C_in= createMat(3,3);
//C_in=loadMat("C_in.dat",3,3);
//m,n,p values
//printf("m,n,p values = \n");
//double **mnp=Matscale(Matmul(s_vec,C_in, 1,3,3),1,3,0.5);
//printMat(mnp,1,3);

//angular direction matrix
printf("Angular Bisector Direction Matrix = \n");
double **dir_Mat= createMat(3,3);
//i_con[0][0]=  -1;      
//i_con[0][1]= mnp[0][1]/s_vec[0][1];    
//i_con[0][2]= mnp[0][0]/s_vec[0][2] ;
//i_con[1][0]= mnp[0][1]/s_vec[0][0] ;
//i_con[1][1]= -1;
//i_con[1][2]= mnp[0][2]/s_vec[0][2];
//i_con[2][0]= mnp[0][0]/s_vec[0][0];
//i_con[2][1]= mnp[0][2]/s_vec[0][1];
//i_con[2][2]= -1;

dir_Mat[0][0]=  -1 ;      
dir_Mat[0][1]= nv/b;    
dir_Mat[0][2]= mv/c;
dir_Mat[1][0]= nv/a;
dir_Mat[1][1]= -1 ;
dir_Mat[1][2]= pv/c;
dir_Mat[2][0]= mv/a;
dir_Mat[2][1]= pv/b;
dir_Mat[2][2]= -1 ;

double **a_dir=Matmul(mat,dir_Mat,2,3,3);
printMat(a_dir,2,3);
printf("\n");

//Angular Normal matrix
printf("Angular Bisector Normal matrix = \n");
double **a_n1,**a_n2,**a_n3,**a_nor,**at_nor;
a_n1=normVec(Matcol(a_dir,2,0));
a_n2=normVec(Matcol(a_dir,2,1));
a_n3=normVec(Matcol(a_dir,2,2));
a_nor=matrix_merge(a_n1,a_n2,a_n3,2,3);
printMat(a_nor,2,3);
printf("\n");
//printf("transpose \n");
//at_nor =transposeMat(a_nor,2,3);
//printMat(at_nor,3,2);

//constant matrix
printf("Angular Bisector constant matrix= \n");
double x1,x2,x3,**a_c_m,**at_c_m,**aam_vec;
x1 = Matdot(a_n1,A,2);
x2 = Matdot(a_n2,B,2);
x3 = Matdot(a_n3,C,2);
a_c_m=matrix(x1,x2,x3);
printMat(a_c_m,1,3);
printf("\n");
//printf("transpose = \n");
//at_c_m = transposeMat(a_c_m,1,3);
//printMat(at_c_m,3,1);

//printing line matrix
at_nor =transposeMat(a_nor,2,3);
at_c_m = transposeMat(a_c_m,1,3);
aam_vec = matrix_2Merge(at_nor,3,2,at_c_m,3,1);
printf("Angular Bisector Line matrix = \n");
printMat(aam_vec, 3, 3);
printf("\n");

//Angular intersection points
double **in_cen=line_intersect(aam_vec,3,3);
printf("Angular Bisector Intersection points = \n");
printMat(in_cen,1,2);
printf("\n");

//Angular Contact points
double **adir_Mat= createMat(3,3);
adir_Mat[0][0]= 0;
adir_Mat[0][1]= nv/b;    
adir_Mat[0][2]= mv/c;
adir_Mat[1][0]= nv/a;
adir_Mat[1][1]= 0;
adir_Mat[1][2]= pv/c;
adir_Mat[2][0]= mv/a;
adir_Mat[2][1]= pv/b;
adir_Mat[2][2]= 0;
double **aa_dir=Matmul(mat,adir_Mat,2,3,3);
printf("Angular Bisector Contact Points = \n");
printMat(aa_dir,2,3);
printf("\n");
printf("**************************** THE END ********************************** \n");
//********************************END OF ANGULAR BISECTOR*****************************

}
