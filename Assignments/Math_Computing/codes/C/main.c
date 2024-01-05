#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libs/matfun.h"

int main() {
	FILE *fp; //file pointer
	double **vert;
	double **A,**B,**C,**D,**M; //declaring matrices names
	//double DM,CM,BM,AM,AC,BD,BC,CD,AB,CB; //side lengths
	int m =2,k=4, n=1; // (mxn),(mxk) matrix
	double l = 6; //length of a side 
	double CMA,BMD,DBC,BCA;//angles of the triangle

// creating Matrix MAT and loading MATRIX from file
vert = createMat(m, k);
vert = loadMat("vertices.dat", m, k);

// Extracting vertices
A = Matcol(vert, m, 0);
B = Matcol(vert, m, 1);
C = Matcol(vert, m, 2);
D = Matcol(vert, m, 3);

// mid-point(A,B)
M = Matsec(A, B, m, n);

//Define side names and corresponding matrices
char *sideNames[] = { "AM", "BM", "CM", "DM", "AB","AC", "AD", "BC", "BD", "CD"};//arr[0]=AM,arr[1]=BM,arr[2]=CM,arr[3]=DM,arr[4]=AB,arr[5]=AC,arr[6]=AD,arr[7]=BC,arr[8]=BD,arr[9]=CD
double ***matrices[] = {&M ,&A, &B,&C,&D};
    
    int z=0;
    double arr[100];
    // Print side lengths using a loop
    for (int i = 0; i <= 4; i++){ 
	    for (int j=0;j<=4;j++){
		    if(i<j){
         arr[z] = Matnorm(Matsub(*matrices[i], *matrices[j], m, n), m);
        z++;
    }}}
   for (int l=z-1;l>=0;l--){
	printf("%s: %.5f\n", sideNames[l],arr[l]);
   }

//finding the angles
CMA =angle(arr[2],arr[0],arr[5]);
BMD =angle(arr[1],arr[3],arr[8]);
DBC =angle(arr[7],arr[6],arr[9]);
BCA =angle(arr[7],arr[5],arr[4]);


//printing the angles
printf("∠ CMA = %lf \n",CMA);
printf("∠ BMD = %lf \n",BMD);
printf("∠ DBC = %lf \n",DBC);
printf("∠ BCA = %lf \n",BCA);


//printf("1.For proving △ AMC ≅ △ BMD \n");
 if ((arr[3] == arr[2]) && (arr[1] == arr[0]) && (CMA == BMD)) {
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

