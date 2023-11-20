
//Function declaration
double Listnorm(avyuh *a);//norm of a vector                        //
double Listdot(avyuh *a, avyuh * b);//inner product                 //
double ListVecdot(sadish *a, sadish * b);//inner product            //
avyuh *Listsub(avyuh *a, avyuh *b);//subtract two matrices          //
avyuh *Listadd(avyuh *a, avyuh *b);//add two matrices               //
avyuh *Listscale(avyuh *a, int m, int n, double k);//scale matrix
avyuh *Listinv(avyuh *mat, int m);//invert an m  x m matrix, m <=3
avyuh *Listmul(avyuh *a, avyuh *b);//multiply matrices a and b      // 
avyuh *rotList(double theta); //rotation matrix                     //
avyuh *normVec(avyuh *a); //normal vector                           //
void circulantList(avyuh *c, int m);
avyuh *Listsec(avyuh *a, avyuh * b, int m, double k);//section formula


avyuh *VertToList(avyuh *a,avyuh *b,avyuh *c);// points to list conversion
avyuh *Listdiag(avyuh *a);//diagonal matrix
avyuh *Listsqrtdiag(avyuh *a);//square root of elements in a diagonal matrix
sadish *xtrtdiag(avyuh *a);//to extract diagonal vector from a matrix
sadish *xtrtsqrtdiag(avyuh *a);//to extract sqrt diagonal vector from a matrix
avyuh *h_stkList(avyuh *a,avyuh *b);//concates list a,b horizontally and transposes the result
avyuh *line_intersect(avyuh *a);//to find intersection point of lines (solving pair of linear equations)
//double readidx(avyuh *a,int m,int n);//to read the value of an element at index (m,n)
//End function declaration


//inner product
double ListVecdot(sadish *a, sadish *b){
	double val = 0;
	sadish  *tempb=b;
	 for(sadish *tempa=a;tempa!=NULL;tempa=tempa->next){
		val += tempa->data*tempb->data;
		//tempa = tempa->next;
		tempb = tempb->next;
	}
	return val;
}
double Listdot(avyuh *a, avyuh *b){
	return ListVecdot(a->vector, b->vector);
}
sadish *ListVecsub(sadish *a, sadish *b){
	sadish *head = (sadish *)malloc(sizeof(sadish)), *c,  *tempb=b;
	c = head; 
	head->next = NULL;
	for(sadish *tempa=a;tempa!=NULL;tempa=tempa->next){
		c->data = tempa->data-tempb->data;
		//tempa = tempa->next;
		tempb = tempb->next;
	if(tempa!=NULL){
		c->next = (sadish *)malloc(sizeof(sadish));
		c->next->next=NULL;
		c= c->next;
	}
	}
	return head;
}

//subtract two matrices
avyuh *Listsub(avyuh *a, avyuh *b){
	avyuh *c= (avyuh *)malloc(sizeof(avyuh)),  *tempb = b, *head; 
	c->next = NULL;
	head = c;
	 for(avyuh *tempa=a;tempa!=NULL;tempa=tempa->next){
		c->vector = ListVecsub(tempa->vector,tempb->vector);
		//tempa = tempa->next;
		tempb = tempb->next;
	if(tempa!=NULL){
		c->next = (avyuh *)malloc(sizeof(avyuh));
		c->next->next=NULL;
		c= c->next;
	}
	}
	return head;
}


sadish *ListVecadd(sadish *a, sadish *b){
        sadish *head = (sadish *)malloc(sizeof(sadish)), *c, *tempb=b;
        c = head; 
        head->next = NULL;
         for(sadish *tempa=a;tempa!=NULL;tempa=tempa->next){
                c->data = tempa->data+tempb->data;
                //tempa = tempa->next;
                tempb = tempb->next;
        if(tempa!=NULL){
                c->next = (sadish *)malloc(sizeof(sadish));
                c->next->next=NULL;
                c= c->next;
        }
        }
        return head;
}

//add two matrices
avyuh *Listadd(avyuh *a, avyuh *b){
        avyuh *c= (avyuh *)malloc(sizeof(avyuh)), *tempb = b, *head; 
        c->next = NULL;
        head = c;
         for(avyuh *tempa=a;tempa!=NULL;tempa=tempa->next){
                c->vector = ListVecadd(tempa->vector,tempb->vector);
                //tempa = tempa->next;
                tempb = tempb->next;
        if(tempa!=NULL){
                c->next = (avyuh *)malloc(sizeof(avyuh));
                c->next->next=NULL;
                c= c->next;
        }
        }
        return head;
}





//norm of a vector
double Listnorm(avyuh *a){
	return sqrt(Listdot(a,a));
}

//rotation matrix
avyuh *rotList(double theta){ 
avyuh *head = createList(2,2), *temp;//create empty 2 x 2 matrix 
sadish *row1, *row2;
double c = cos(theta), s = sin(theta);

row1 = head->vector;
row2 = head->next->vector;
row1->data = c;
row1->next->data = -s;
row2->data = s;

row2->next->data = c;
return head;

}

avyuh *VertToList(avyuh *a,avyuh *b,avyuh *c)
{
	int i;
avyuh *head = createList(2,3);
sadish *r1,*r2;
r1      =head->vector;
r2      =head->next->vector;

r1->data             =Vecind(a->vector,0)->data;
r1->next->data       =Vecind(b->vector,0)->data;
r1->next->next->data =Vecind(c->vector,0)->data;

r2->data             =Vecind(a->vector,1)->data;
r2->next->data       =Vecind(b->vector,1)->data;
r2->next->next->data =Vecind(c->vector,1)->data;

return head;
}


//Matrix multiplication
avyuh *Listmul(avyuh *a, avyuh *b){
	avyuh *c=(avyuh *)malloc(sizeof(avyuh)),*atemp, *btemp;
	avyuh *head = c;
	sadish *cvec=(sadish *)malloc(sizeof(sadish));
	head->next = NULL;
	for(atemp=a; atemp !=NULL; atemp=atemp->next){
		cvec=(sadish *)malloc(sizeof(sadish));
		cvec->next = NULL;
		c->vector = cvec;  
		for(btemp=transposeList(b); btemp !=NULL; btemp=btemp->next){
			cvec->data=ListVecdot(atemp->vector, btemp->vector);//inner product
		if(btemp->next!=NULL){
			cvec->next = (sadish *)malloc(sizeof(sadish));
			cvec->next->next=NULL;
			cvec= cvec->next;
		}
		}

		if(atemp->next!=NULL){
			c->next = (avyuh *)malloc(sizeof(avyuh));
			c->next->next=NULL;
			c= c->next;
		}
	}
return head;
}


//function to find the diagonal matrix
avyuh *Listdiag(avyuh *a){
	avyuh *head= (avyuh *)malloc(sizeof(avyuh)), *alist;
	head->vector= xtrtdiag(a);
	head->next= NULL;
return head;
}
//function to extract diagonal vector from the matrix
sadish *xtrtdiag(avyuh *a){
	int i = 0,j = 0;//dummy integers
	int m = 2;
	//avyuh *head= (avyuh *)malloc(sizeof(avyuh)), *alist;
	sadish *head=(sadish *)malloc(sizeof(sadish)),*temp, *btemp;
	btemp = head;
	head->next= NULL;
	for (avyuh *alist = a; alist != NULL; alist= alist->next){
		temp = Vecind(alist->vector,i);//getting address of the nth column
		btemp->data = temp->data;
		i++;
	if(alist->next !=NULL){
		btemp->next = (sadish *)malloc(sizeof(sadish));
		btemp->next->next=NULL;
		btemp= btemp->next;
	}
	}
return head;
}



// function to find sqrt diagonal matrix
avyuh *Listsqrtdiag(avyuh *a){
        avyuh *head= (avyuh *)malloc(sizeof(avyuh)), *alist;
        head->vector= xtrtsqrtdiag(a);
        head->next= NULL;
return head;   
} 
//function to extract sqrt diagonal vector from a matrix
sadish *xtrtsqrtdiag(avyuh *a){
        int i = 0,j = 0;//dummy integers
        int m = 2;
        //avyuh *head= (avyuh *)malloc(sizeof(avyuh)), *alist;
        sadish *head=(sadish *)malloc(sizeof(sadish)),*temp, *btemp;
        btemp = head;
        head->next= NULL;
        for (avyuh *alist = a; alist != NULL; alist= alist->next){
                temp = Vecind(alist->vector,i);//getting address of the nth column
                btemp->data = sqrt(temp->data);
                i++;
        if(alist->next !=NULL){
                btemp->next = (sadish *)malloc(sizeof(sadish));
                btemp->next->next=NULL;
                btemp= btemp->next;
        }
        }
return head;
}



//function to concate lists a and b horizontally and transpose the result
avyuh *h_stkList(avyuh *a, avyuh *b){
	avyuh *c= (avyuh *)malloc(sizeof(avyuh)), *tempb = b, *head; 
	c->next = NULL;
	head = c;
	 for(avyuh *tempa=a;tempa!=NULL;tempa=tempa->next){
		c->vector = tempa->vector;
		//tempa = tempa->next;
	
		c->next = (avyuh *)malloc(sizeof(avyuh));
		c->next->next=NULL;
		c= c->next;

	
	}
	c->vector=b->vector;
	return transposeList(head);
}



// function to find the point of intersection of lines
avyuh  *line_intersect(avyuh *a)
{
avyuh *temp = transposeList(a);
double a1,a2,b1,b2,c1,c2;
a1=temp->vector->data;
a2=temp->vector->next->data;
b1=temp->next->vector->data;
b2=temp->next->vector->next->data;
c1=-temp->next->next->vector->data;
c2=-temp->next->next->vector->next->data;
//printf("%lf,%lf,%lf,%lf,%lf,%lf,",a1,a2,b1,b2,c1,c2);
avyuh *head = createList(1,2);
head->vector->data=(b1*c2-b2*c1)/(a1*b2-a2*b1);
head->vector->next->data=(c1*a2-c2*a1)/(a1*b2-a2*b1);
return head;

}


/*
//function to read the value of element at index (m,n)
double readidx(avyuh *a,int m,int n)
{
avyuh *head=a;

for(int i=0;i<m;i++)
{
head=head->next;
}
//for (int j=0;j<n;j++)
//{
//head->vector=head->vector->next;
//}
double val = Vecind(head->vector,n)->data;
return val;
}
*/


