#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libs/matfun.h"

int  main()
{

FILE *fp; //file pointer
double **A,**B,**C,**D,**M; //declaring matrices names
double DM,CM,BM,AM,AC,BD,BC,CD,AB,CB; //side lengths
int m =2, n=1; // (mxn) matrix
double l = 6; //length of a side 
double CMA,BMD,DBC,BCA;

//creating Matrix B and loading MATRIX from file 
B = createMat(m,n);
B = loadMat("B.dat",m, n);

//creating Matrix C and loading Matrix from file
C = createMat(m,n);
C = loadMat("C.dat",m,n);

//creating Matrix A and loading Matrix from file
A =createMat(m,n);
A =loadMat("A.dat",m, n);

//creating Matrix D and loading Matrix from file
D=createMat(m,n);
D =loadMat("D.dat",m,n);


//mid-point(A,B)
M = Matsec(A,B,m,n);


//printing the matrices
/*
printMat(B,2,1);
printMat(C,2,1);
printMat(A,2,1);
printMat(D,2,1);
printMat(M,2,1);
printf("l = %lf \n",l); 
*/


//lengths of All sides 
DM = Matnorm(Matsub(D,M,m,n),m);
CM = Matnorm(Matsub(C,M,m,n),m);
BM = Matnorm(Matsub(B,M,m,n),m);
AM = Matnorm(Matsub(A,M,m,n),m);
AC = Matnorm(Matsub(A,C,m,n),m);
BD = Matnorm(Matsub(B,D,m,n),m);
BC = Matnorm(Matsub(B,C,m,n),m);
CD = Matnorm(Matsub(C,D,m,n),m);
AB = Matnorm(Matsub(A,B,m,n),m);
CB = Matnorm(Matsub(C,B,m,n),m);


//printing side lengths
/*
printf("DM = %lf \n",DM);
printf("CM = %lf \n",CM);
printf("BM = %lf \n",BM);
printf("AM = %lf \n",AM);
printf("AC = %lf \n",AC);
printf("BD = %lf \n",BD);
printf("BC = %lf \n",BC);
printf("CD = %lf \n",CD);
printf("AB = %lf \n",AB);
printf("CB = %lf \n",CB);
*/

//finding the angles
CMA =angle(CM,AM,AC);
BMD =angle(BM,DM,BD);
DBC =angle(BC,BD,CD);
BCA =angle(BC,AC,AB);

/*
printf("∠ CMA = %lf \n",CMA);
printf("∠ BMD = %lf \n",BMD);
printf("∠ DBC = %lf \n",DBC);
printf("∠ BCA = %lf \n",BMD);
*/

//printf("1.For proving △ AMC ≅ △ BMD \n");
 if ((DM == CM) && (BM == AM) && (CMA == BMD)) {
        printf("△ AMC ≅ △ BMD (congruent By SAS Congruency)\n");
    } else {
        printf("△ AMC ≇ △ BMD (is NOT congruent)\n");
    };

//printf("2.△ DBC IS right angled At ∠ B \n");
      if (round(DBC) == 90) {
        printf("△ DBC IS right angled At ∠ B\n");
    } else {
        printf("△ DBC IS NOT right angled At ∠ B\n");
    }

//printf("3.△ DBC ≅ △ ACB \n");
    if (CB == BC && BD == AC && DBC == BCA) {
    	printf("△ DBC ≅ △ ACB (congruent By SAS Congruency)\n");
    } else {
        printf("△ DBC ≇ △ ACB (is NOT congruent)\n");
    }

//printf("4.CM = AB/2 \n");
    if (CM == AB / 2 || 2 * CM == AB) {
        printf("CM = AB/2\n");
    } else {
        printf("CM ≠ AB/2\n");
    }
return 0;
}
