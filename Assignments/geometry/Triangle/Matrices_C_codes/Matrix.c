#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"coeffs-mat.h"
int main() {
    int m = 2;  // You can set the dimensions of your matrix
    int n = 3;

// Create a matrix
	double **G_v = createMat(m, n);
	double **C_m= createMat(3,3);
	double **C_mid= createMat(3,3);
	double **R_o=createMat(2,2);
	double **C_mid_dir= createMat(3,3);
	double **C_alt= createMat(3,3);
	double **C_in= createMat(3,3);
	double **cont_mat= createMat(3,3);
	double **i_con= createMat(3,3);

// Reading matrix data from .dat file
G_v=loadtxt("vertices.dat",m,n);
C_m=loadtxt("circulant_matrix.dat",3,3);
C_mid=loadtxt("midpoint_matrix.dat",3,3);
C_mid_dir=loadtxt("median_direction_matrix.dat",3,3);
R_o=loadtxt("rotation_matrix.dat",2,2);
C_alt=loadtxt("altitude_circulant_matrix.dat",3,3);
C_in=loadtxt("angle_bisector.dat",3,3);



	// **********************VECTORS*****************************
	double **G_dir=matmul(G_v,C_m,2,3,3);
	double **G_n=matmul(R_o,G_dir,2,2,3);
	double **G_con=diag(matmul(transpose(G_n,2,3),G_v,3,2,3),3,3);
	double **G_dis=sqrt_diag(matmul(transpose(G_dir,2,3),G_dir,3,2,3),3,3);
	double **G_line=h_concat(transpose(G_n,2,3),transpose(G_con,1,3),3,2,3,1);
	//***********************MEDIANS*******************************
	double **G_mid=nmatmul(0.5,matmul(G_v,C_mid,2,3,3),2,3);
	double **G_med_dir=matmul(G_v,C_mid_dir,2,3,3);
	double **G_n_med = matmul(R_o,G_med_dir,2,2,3);
	double **cmat_med=diag(matmul(transpose(G_n_med,2,3),G_v,3,2,3),3,3);
	double **linemat_med=h_concat(transpose(G_n_med,2,3),transpose(cmat_med,1,3),3,2,3,1);	
	double **G_G=line_intersect(linemat_med,3,3);
	//***********************ALTITUDE*******************************
	double **G_n_alt=matmul(G_v,C_alt,2,3,3);
	double **cmat_alt=diag(matmul(transpose(G_n_alt,2,3),G_v,3,2,3),3,3);
	double **linemat_alt=h_concat(transpose(G_n_alt,2,3),transpose(cmat_alt,1,3),3,2,3,1);	
	double **G_H=line_intersect(linemat_alt,3,3);
	//******************PERPENDICULAR BISECTOR**********************
	double **cmat_perp_bis=diag(matmul(transpose(G_n_alt,2,3),G_mid,3,2,3),3,3);
	double **linemat_perp_bis=h_concat(transpose(G_n_alt,2,3),transpose(cmat_perp_bis,1,3),3,2,3,1);	 
	double **G_O=line_intersect(linemat_perp_bis,3,3);
	//********************ANGULAR  BISECTOR************************
	double **secvec=nmatmul(0.5, matmul(G_dis,C_in, 1,3,3),1,3);
	i_con[0][0]=  -1 ;      i_con[0][1]=secvec[0][2]/G_dis[0][2];     i_con[0][2]= secvec[0][1]/G_dis[0][0] ;
	i_con[1][0]=secvec[0][2]/G_dis[0][1] ;       i_con[1][1]= -1 ;       i_con[1][2]=secvec[0][0]/G_dis[0][0];
	i_con[2][0]= secvec[0][1]/G_dis[0][1];      i_con[2][1]= secvec[0][0]/G_dis[0][2] ;    i_con[2][2]= -1 ;
	double **a_dir=matmul(G_v,i_con,2,3,3);
	double **a_nor=matmul(R_o,a_dir,2,2,3);
	double **a_cof=diag(matmul(transpose(a_nor,2,3),G_v,3,2,3),3,3);
	double **a_line=h_concat(transpose(a_nor,2,3),transpose(a_cof,1,3),3,2,3,1);
	double **a_i = line_intersect(a_line,3,3);


        cont_mat[0][0]=  0 ;    cont_mat[0][1]=secvec[0][2]/G_dis[0][2];     cont_mat[0][2]= secvec[0][1]/G_dis[0][0] ;
        cont_mat[1][0]=secvec[0][2]/G_dis[0][1] ;       cont_mat[1][1]= 0 ;      cont_mat[1][2]=secvec[0][0]/G_dis[0][0];
        cont_mat[2][0]= secvec[0][1]/G_dis[0][1];      cont_mat[2][1]= secvec[0][0]/G_dis[0][2] ;    cont_mat[2][2]= 0 ;
        double **G_i=matmul(G_v,cont_mat,2,3,3);
        double **G_imid=nmatmul(0.5,matmul(G_i,C_mid,2,3,3),2,3);
        double **G_idir_alt=matmul(G_i,C_alt,2,3,3);
        double **cmat_iperp_bis=diag(matmul(transpose(G_idir_alt,2,3),G_imid,3,2,3),3,3);
        double **linemat_imed=h_concat(transpose(G_idir_alt,2,3),transpose(cmat_iperp_bis,1,3),3,2,3,1);
        double **G_I=line_intersect(linemat_imed,3,3);


// Printing matrices . 	
printf("**************************** Vectors ********************************** \n");
printf("direction matrix= \n");
print(G_dir,2,3);
printf("normal matrix= \n");
print(G_n,2,3);
printf("constant matrix= \n");
print(G_con,1,3);
printf("distance matrix= \n");
print(G_dis,1,3);
printf("line  matrix= \n");
print(G_line,3,3);

printf("**************************** Medians ********************************** \n");
printf("midpoint matrix= \n");
print(G_mid,2,3);
printf("median direction  matrix= \n");
print(G_med_dir,2,3);
printf("median normal matrix= \n");
print(G_n_med,2,3);
printf("median constant matrix= \n");
print(cmat_med,1,3);
printf("median line matrix= \n");
print(linemat_med,3,3);
printf("Centroid = \n");
print(G_G,1,2);

printf("**************************** Altitude ********************************** \n");
printf("altitude normal matrix= \n");
print(G_n_alt,2,3);
printf("altitude constant matrix= \n");
print(cmat_alt,1,3);
printf("Altitude line matrix= \n");
print(linemat_alt,3,3);
printf("Orthocentre = \n");
print(G_H,1,2);

printf("**************************** Perpendicular Bisector ********************************** \n");
printf("perp_bisect constant  matrix= \n");
print(cmat_perp_bis,1,3);
printf("perp_bis line matrix= \n");
print(linemat_perp_bis,3,3);
printf("Circumcentre = \n");
print(G_O,1,2);

printf("**************************** Angular Bisector ********************************** \n");
printf("m,n,p values = \n");
print(secvec,1,3);
printf("ang_bis direction matrix= \n");
print(a_dir,2,3);
printf("ang_bis normal matrix= \n");
print(a_nor,2,3);
printf("ang_bis constant matrix= \n");
print(a_cof,1,3);
printf("ang_bis line  matrix= \n");
print(a_line,3,3);
printf("ang_bis intersection points= \n");
print(a_i,1,2);
printf("Incentre = \n");
print(G_I,1,2);

printf("contact points = \n");
print(G_i,2,3);
printf("************************************************************** \n");
}



