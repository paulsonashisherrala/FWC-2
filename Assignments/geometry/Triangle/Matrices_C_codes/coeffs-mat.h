//Function declaration
double **createMat(int m,int n);
void readMat(double **p, int m,int n);
void print(double **p,int m,int n);
double **loadtxt(char *str,int m,int n);
double linalg_norm(double **a, int m);
double **linalg_sub(double **a, double **b, int m, int n);
double **linalg_inv(double **mat, int m);
double **matmul(double **a, double **b, int m, int n, int p);
double **transpose(double **a,  int m, int n);
void uniform(char *str, int len);
void gaussian(char *str, int len);
double mean(char *str);
double **diag(double **a,int m,int n);
double **sqrt_diag(double **a,int m,int n);
double **nmatmul(double num ,double **a,int m,int n);
double **h_concat(double **a, double **b, int m, int n, int p, int q);
double **line_intersect(double **a,int m,int n);
//End function declaration


//Defining the function for matrix creation
double **createMat(int m,int n)
{
 int i;
 double **a;
 
 //Allocate memory to the pointer
a = (double **)malloc(m * sizeof( *a));
    for (i=0; i<m; i++)
         a[i] = (double *)malloc(n * sizeof( *a[i]));

 return a;
}
//End function for matrix creation


//Defining the function for reading matrix 
void readMat(double **p, int m,int n)
{
 int i,j;
 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
   scanf("%lf",&p[i][j]);
  }
 }
}
//End function for reading matrix

//Read  matrix from file
double **loadtxt(char *str,int m,int n)
{
FILE *fp;
double **a;
int i,j;


a = createMat(m,n);
fp = fopen(str, "r");

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
   fscanf(fp,"%lf",&a[i][j]);
  }
 }
//End function for reading matrix from file

fclose(fp);
 return a;

}


//Defining the function for printing
void print(double **p, int m,int n)
{
 int i,j;

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  printf("%lf ",p[i][j]);
 printf("\n");
 }
 printf("\n");
}
//End function for printing

//Defining the function for norm

double linalg_norm(double **a, int m)
{
int i;
double norm=0.0;

 for(i=0;i<m;i++)
 {
norm = norm + a[i][0]*a[i][0];
}
return sqrt(norm);

}
//End function for norm

//Defining the function for difference of matrices

double **linalg_sub(double **a, double **b, int m, int n)
{
int i, j;
double **c;
c = createMat(m,n);

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
c[i][j]= a[i][j]-b[i][j];
  }
 }
return c;

}
//End function for difference of matrices

//Defining the function for inverse of 2x2 matrix


double **linalg_inv(double **mat, int m)
{
double **c, det;
c = createMat(m,m);

det = mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];

c[0][0] = mat[1][1]/det;
c[0][1] = -mat[1][0]/det;
c[1][0] = -mat[0][1]/det;
c[1][1] = mat[0][0]/det;

return c;

}
// End  function for inverse of 2x2 matrix


//Defining the function for difference of matrices

double **matmul(double **a, double **b, int m, int n, int p)
{
int i, j, k;
double **c, temp =0;
c = createMat(m,p);

 for(i=0;i<m;i++)
 {
  for(k=0;k<p;k++)
  {
    for(j=0;j<n;j++)
    {
	temp= temp+a[i][j]*b[j][k];
    }
	c[i][k]=temp;
	temp = 0;
  }
 }
return c;

}
//End function for difference of matrices

// function for horizontal concatination of a mxn & mx1 matrix 
double **h_concat(double **a, double **b, int m, int n, int p,int q)
{
	if(m==p)
	{
	int i,j;
	double **c=createMat(m,n+q);
	for(i=0;i<m;i++)
		{
		for(j=0;j<n;j++)
			{
			c[i][j]=a[i][j];
			}
		}

	for(i=0;i<m;i++)
		{
		c[i][n]=b[i][0];
	
		}

	return c;
	}
	else
		{
		printf("Dimensional error");
		}
}
// end of horizontal concatination function

//Defining the function for transpose of matrix

double **transpose(double **a,  int m, int n)
{
int i, j;
double **c;
//printf("I am here");
c = createMat(n,m);

 for(i=0;i<n;i++)
 {
  for(j=0;j<m;j++)
  {
c[i][j]= a[j][i];
//  printf("%lf ",c[i][j]);
  }
 }
return c;

}
//End function for transpose of matrix

//defining function for generating diagonal matrix
double **diag(double **a ,int m ,int n)
{
if(m==n)
{
int i;
double **d;
d=createMat(1,m);
	for(i=0;i<m;i++)
	{
	d[0][i]=a[i][i];
	}
	return d;
}
else
{
	printf("Dimensional error");
}
}
//end function for diagonal matrix generation

//defining function for generating square root diagonal matrix
double **sqrt_diag(double **a ,int m ,int n)
{
if(m==n)
{
int i;
double **d;
d=createMat(1,m);
        for(i=0;i<m;i++)
        {
        d[0][i]=sqrt(a[i][i]);
        }
        return d;
} 
else
{ printf("Dimensional error");
}
}
//end of  function for diagonal matrix generation



// defining a function for finding the intersection point between two lines
double **line_intersect(double **a,int m,int n)
{

double **c=createMat(1,2);
a[0][2]=-a[0][2];
a[1][2]=-a[1][2];
c[0][0]= (a[0][1]*a[1][2]-a[1][1]*a[0][2])/(a[0][0]*a[1][1]-a[1][0]*a[0][1]);
c[0][1]= (a[0][2]*a[1][0]-a[1][2]*a[0][0])/(a[0][0]*a[1][1]-a[1][0]*a[0][1]);

return c;
}
//End of function for finding point of  intersection of lines 

// defining a function to multiply a number n with matrix mxn
double **nmatmul(double num,double **a ,int m ,int n)
{

int i,j;
double **d;
d=createMat(m,n);
        for(i=0;i<m;i++)
        {
         	for(j=0;j<n;j++)
		{
		d[i][j]=num*a[i][j];
		}
        }
        return d;

}
//end of function for number matrix multiplication



//Defining the function for generating uniform random numbers
void uniform(char *str, int len)
{
int i;
FILE *fp;
fp = fopen(str,"w");
//Generate numbers
for (i = 0; i < len; i++)
{
fprintf(fp,"%lf\n",(double)rand()/RAND_MAX);
}
fclose(fp);

}
//End function for generating uniform random numbers

//Defining the function for calculating the mean of random numbers
double mean(char *str)
{
int i=0,c;
FILE *fp;
double x, temp=0.0;

fp = fopen(str,"r");
//get numbers from file
while(fscanf(fp,"%lf",&x)!=EOF)
{
//Count numbers in file
i=i+1;
//Add all numbers in file
temp = temp+x;
}
fclose(fp);
temp = temp/(i-1);
return temp;

}
//End function for calculating the mean of random numbers

//Defining the function for generating Gaussian random numbers
void gaussian(char *str, int len)
{
int i,j;
double temp;
FILE *fp;

fp = fopen(str,"w");
//Generate numbers
for (i = 0; i < len; i++)
{
temp = 0;
for (j = 0; j < 12; j++)
{
temp += (double)rand()/RAND_MAX;
}
temp-=6;
fprintf(fp,"%lf\n",temp);
}
fclose(fp);

}
//End function for generating Gaussian random numbers
