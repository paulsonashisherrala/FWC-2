#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libs/listgen.h"
#include "libs/listfun.h"

int  main()
{
avyuh *A,*B,*C,*D,*M;
double l=6; //lengths of a side 
int m =2,n=1,k=4; //(mxn) matrices
double CMA,BMD,DBC,BCA; //angles 

//load matrix from file
avyuh *vert=loadList("vertices.dat",m,k);

//Triangle vertices
A = Listcol(vert,0);
B = Listcol(vert,1);
C = Listcol(vert,2);
D = Listcol(vert,3);

//Midpoint
M = Listsec(A,B,1);

//printing the lists
//printList(vert);
printList(A);
printList(B);
printList(C);
printList(D);
printList(M);

// calculating side lengths using loop
avyuh *points[] = {M,A,B,C,D};
char *sideNames[] = {"AM", "BM", "CM", "DM", "AB", "AC", "AD", "BC", "BD", "CD"};
    int z=0;
    double arr[100];
    for (int i = 0; i < 4;i++){
            for (int j=0;j<=4;j++){
                    if(i<j){
        arr[z] = Listnorm(Listsub(points[i],points[j]));
        z++;
        }
         }
    }
   for(int l=z-1;l>=0;l--){
          printf("%s: %.5f\n",sideNames[l],arr[l]);
  }

//finding the angles 
CMA = angle(arr[2],arr[0],arr[5]);
BMD = angle(arr[1],arr[3],arr[8]);
DBC = angle(arr[7],arr[6],arr[9]);
BCA = angle(arr[7],arr[5],arr[4]);

//printing the angles
printf("\n∠ CMA = %lf \n",CMA);
printf("∠ BMD = %lf \n",BMD);
printf("∠ DBC = %lf \n",DBC);
printf("∠ BCA = %lf \n",BCA);

//printf("1.For proving △ AMC ≅ △ BMD \n");
 if ((arr[3] == arr[2]) && (arr[1] == arr[0]) && (CMA == BMD)) {
        printf("\n△ AMC ≅ △ BMD (congruent By SAS Congruency)\n");
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
    if (arr[7]== arr[7] && arr[8] == arr[5] && DBC == BCA) {
        printf("△ DBC ≅ △ ACB (congruent By SAS Congruency)\n");
    } else {
        printf("△ DBC ≇ △ ACB (is NOT congruent)\n");
    }

//printf("4.CM = AB/2 \n");
    if (arr[2] == arr[4] / 2 || 2 * arr[2] == arr[4]) {
        printf("CM = AB/2\n");
    } else {
        printf("CM ≠ AB/2\n");
    }
return 0;
}
