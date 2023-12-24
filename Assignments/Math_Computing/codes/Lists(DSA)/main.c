#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libs/listgen.h"
#include "libs/listfun.h"

int  main()
{
double DM,CM,BM,AM,AC,BD,BC,CD,AB,CB; //side and their lengths
int m =2, n=1; //(mxn) matrices
double l=6; //lengths of a side 
double CMA,BMD,DBC,BCA; //angles 

//load matrix from file
avyuh *A=loadList("A.dat",m,n);
avyuh *B=loadList("B.dat",m,n);
avyuh *C=loadList("C.dat",m,n);
avyuh *D=loadList("D.dat",m,n);

//print length of side 
avyuh *M=Listsec(A,B,1);

//printing the lists
/*
printList(A);
printList(B);
printList(C);
printList(D);
printList(M);
*/

//lengths of All sides 
DM = Listnorm(Listsub(D,M));
CM = Listnorm(Listsub(C,M));
BM = Listnorm(Listsub(B,M));
AM = Listnorm(Listsub(A,M));
AC = Listnorm(Listsub(A,C));
BD = Listnorm(Listsub(B,D));
BC = Listnorm(Listsub(B,C));
CD = Listnorm(Listsub(C,D));
AB = Listnorm(Listsub(A,B));
CB = Listnorm(Listsub(C,B));

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
