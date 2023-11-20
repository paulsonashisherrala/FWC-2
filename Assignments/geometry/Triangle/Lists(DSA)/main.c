#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libs/listgen.h"
#include "libs/listfun.h"

int  main()
{
FILE *fp; //file pointer
double val;//for reading file data
int m =2, n=3, i, j;
//load matrix from file
avyuh *G_v= loadList("vertices.dat", m, n);
avyuh *C_mid=loadList("C_mid.dat",3,3);
avyuh *C_mid_dir=loadList("C_mid_dir.dat",3,3);
avyuh *C_alt=loadList("C_alt.dat",3,3);
avyuh *C_in=loadList("C_in.dat",3,3);

//Triangle vertices
avyuh *A = Listcol(G_v,0);
avyuh *B = Listcol(G_v,1);
avyuh *C = Listcol(G_v,2);
//Direction vectors
avyuh *m1 = Listsub(A, B);
avyuh *m2 = Listsub(B, C);
avyuh *m3 = Listsub(C, A);
//Triangle sides;
double c = Listnorm(m1);
double a = Listnorm(m2);
double b = Listnorm(m3);
//Rotation matrix
avyuh *R_o = rotList(M_PI/2);

//************************************VECTORS***************************************
		avyuh *G_dir=VertToList(m1,m2,m3);
		avyuh *G_n=Listmul(R_o,G_dir);
		avyuh *G_con= Listdiag(Listmul(transposeList(G_n),G_v));
		avyuh *G_dis= Listsqrtdiag(Listmul(transposeList(G_dir),G_dir));
		avyuh *G_line=h_stkList(G_n,G_con);
//************************************MEDIANS***************************************
		avyuh *G_mid=Listmul(G_v,C_mid);
		avyuh *G_med_dir=Listmul(G_v,C_mid_dir);
		avyuh *G_n_med = Listmul(R_o,G_med_dir);
		avyuh *cmat_med= Listdiag(Listmul(transposeList(G_n_med),G_v));
		avyuh *linemat_med=h_stkList(G_n_med,cmat_med);
		avyuh *G_G=line_intersect(linemat_med);
//************************************ALTITUDE**************************************
		avyuh *G_n_alt = Listmul(G_v,C_alt);
                avyuh *cmat_alt= Listdiag(Listmul(transposeList(G_n_alt),G_v));
                avyuh *linemat_alt=h_stkList(G_n_alt,cmat_alt);
                avyuh *G_H=line_intersect(linemat_alt);
//******************************PERPENDICULAR BISECTOR******************************
		avyuh *cmat_perp_bis= Listdiag(Listmul(transposeList(G_n_alt),G_mid));	
	 	avyuh *linemat_perp_bis=h_stkList(G_n_alt,cmat_perp_bis);
	 	avyuh *G_O=line_intersect(linemat_perp_bis);
//*****************************ANGULAR BISECTOR*************************************
		avyuh *secvec=Listmul(G_dis,C_in);
		avyuh *i_con = createList(3,3);
        sadish *r1,*r2,*r3;
        r1=i_con->vector;
        r2=i_con->next->vector;
        r3=i_con->next->next->vector;
r1->data= -1;  
r1->next->data=secvec->vector->next->next->data/G_dis->vector->next->next->data;
r1->next->next->data=secvec->vector->next->data/G_dis->vector->data;
r2->data=secvec->vector->next->next->data/G_dis->vector->next->data;
r2->next->data=-1;
r2->next->next->data=secvec->vector->data/G_dis->vector->data;
r3->data=secvec->vector->next->data/G_dis->vector->next->data;
r3->next->data=secvec->vector->data/G_dis->vector->next->next->data;
r3->next->next->data=-1;
		avyuh *a_dir=Listmul(G_v,i_con);
		avyuh *a_nor=Listmul(R_o,a_dir);
                avyuh *a_cof= Listdiag(Listmul(transposeList(a_nor),G_v));
                avyuh *a_line=h_stkList(a_nor,a_cof);
                avyuh *a_i=line_intersect(a_line);
		avyuh *cont_mat=i_con;
cont_mat->vector->data=0;
cont_mat->next->vector->next->data=0;
cont_mat->next->next->vector->next->next->data=0;
		avyuh *G_i=Listmul(G_v,cont_mat);
		avyuh *G_imid=Listmul(G_i,C_mid);
		avyuh *G_idir_alt = Listmul(G_i,C_alt);
		avyuh *cmat_iperp_bis= Listdiag(Listmul(transposeList(G_idir_alt),G_imid));
		avyuh *linemat_imed=h_stkList(G_idir_alt,cmat_iperp_bis);
		avyuh *G_I=line_intersect(linemat_imed);

//writeidx(G_v,1,2,7)
// Printing lists . 	
printf("**************************** Vectors ********************************** \n");
printf("vertices matrix= \n");
printList(G_v);
printf("direction matrix= \n");
printList(G_dir);
printf("normal matrix= \n");
printList(G_n);
printf("constant matrix= \n");
printList(G_con);
printf("distance matrix= \n");
printList(G_dis);
printf("line  matrix= \n");
printList(G_line);



printf("**************************** Medians ********************************** \n");
printf("midpoint matrix= \n");
printList(G_mid);
printf("median direction  matrix= \n");
printList(G_med_dir);
printf("median normal matrix= \n");
printList(G_n_med);
printf("median constant matrix= \n");
printList(cmat_med);
printf("median line matrix= \n");
printList(linemat_med);
printf("Centroid = \n");
printList(G_G);


printf("**************************** Altitude ********************************** \n");
printf("altitude normal matrix= \n");
printList(G_n_alt);
printf("altitude constant matrix= \n");
printList(cmat_alt);
printf("Altitude line matrix= \n");
printList(linemat_alt);
printf("Orthocentre = \n");
printList(G_H);


printf("**************************** Perpendicular Bisector ********************************** \n");
printf("perp_bisect constant  matrix= \n");
printList(cmat_perp_bis);
printf("perp_bis line matrix= \n");
printList(linemat_perp_bis);
printf("Circumcentre = \n");
printList(G_O);
//printf("%lf,%lf",readidx(G_v,1,2),readidx(G_v,0,2));


printf("**************************** Angular Bisector ********************************** \n");
printf("m,n,p values = \n");
printList(secvec);
printf("ang_bis direction matrix= \n");
printList(a_dir);
printf("ang_bis normal matrix= \n");
printList(a_nor);
printf("ang_bis constant matrix= \n");
printList(a_cof);
printf("ang_bis line  matrix= \n");
printList(a_line);
printf("ang_bis intersection points= \n");
printList(a_i);
printf("Incentre = \n");
printList(G_I);
printf("contact points = \n");
printList(G_i);
return 0;
}
