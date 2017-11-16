#include <stdio.h>
#include "function.h"
#include "global.h"


gauss(n,a,p,q)
int n;
double a[MAXG][MAXG];
int p[MAXG],q[MAXG];
{
FILE *fp;
int finished,pivot_found,it;
int i,j,k,itemp,last_row,last_col;
double s[MAXG];
double z,dtemp,vtemp[MAXG];
coordinate epslon;

epslon  = 0.000001;


last_row = 0;
last_col = 0;

  for (i=0;i<n;i++)
   {
    p[i]=i;
    q[i]=i;
    s[i]=0;
   }

  for (i=0;i<n;i++)
  {
    do
    {
     for (j=0;j<n;j++)
     if (fabs(a[p[i]][j]) > s[p[i]])
	s[p[i]] = fabs(a[p[i]][j]);
     if (s[p[i]] < epslon) 
	{
	p[n-1-last_row]=p[i];
	p[i] = n-1-last_row;
	s[p[n-1-last_row]]=1; /* hack to avoid divided by 0 */
	last_row += 1;
	}
    }
    while ((s[p[i]] < epslon)  && (i< n-1-last_row));
   }
  

  k=0;
  finished = 0;
  while (!finished)  /* for k = 1 to n- 1 in the book */
  {
  
     pivot_found = 0;

     while ( !pivot_found && (k <= n-1-last_col))
     {   
	 j =k;
         dtemp=0;
	 for (i=k;i<n;i++)
	    if ((fabs(a[p[i]][k])/s[p[i]]) > dtemp)
	       {
		   dtemp = fabs(a[p[i]][k]);
		   j=i;    /* j is the row number of the pivot */
	       }
	 if (fabs(a[p[j]][k]) > epslon)
	    {
		itemp = p[k];
		p[k] = p[j];
		p[j] = itemp;
		pivot_found = 1;
	    }
	 else
	    {
	 /*    if (k != n-1-last_col) */
	      {
		for (it= 0; it< n; it++)
		   {
		       vtemp[it] = a[p[it]][k];
		       a[p[it]][k] = a[p[it]][n-1-last_col];
		       a[p[it]][n-1-last_col] =  vtemp[it];
		   }
		q[n-1-last_col] = q[k];
		q[k] = n-1-last_col; /* no indirection here because last col is */
	      }				/* never moved twice */
		last_col += 1;
	    }
     }
     if (pivot_found)
     for (i=k+1; i<n;i++)
	   {
	       z = a[p[i]][k]/a[p[k]][k];
	       a[p[i]][k] =0;
	       for (j=k+1;j<n;j++)
		  a[p[i]][j] -= z*a[p[k]][j];
	   }
     if ((k >  n-1-last_col ) || (k> n-1-last_row))
	finished = 1;
/*       for (i=0;i<n;i++)
	  {

	      for (j=0; j< n;j++)
		 printf("%f   ", a[p[i]][j]);
	      printf("\n");
	  }
       
       for (i=0;i< n;i++)
	  printf("p[i]  %d     q[i]  %d\n",p[i],q[i]);
*/
     k++;
 }

 fp = fopen ("gauss_result","w");
for (i=0;i<n;i++)
  {
      fprintf(fp,"{");
    for (j=0; j< n;j++)
     fprintf(fp,"%f,   ", a[p[i]][j]);
 fprintf(fp,"}\n");
   }

  for (i=0;i< n;i++)
   fprintf(fp,"p[i]  %d     q[i]  %d\n",p[i],q[i]);
 fclose(fp);
}

 /* main()
{
       int i,j,p[MAXG],q[MAXG];
       double a[MAXG][MAXG];
       float temp;

       int num;

       scanf("%d",&num);

       for (i=0;i<num;i++)
          for (j=0; j< num;j++)
             {
                 scanf("%f", &temp);
                 a[i][j] = temp;
             }

       for (i=0;i<num;i++)
          {
            printf("{");
              for (j=0; j< num;j++)
                 printf("%f,   ", a[i][j]);
            printf("}\n");
          }
printf("====================================\n");

       for (i=0;i<num;i++)
          {
              p[i] = i;
              q[i] = i;
          }
       gauss(num,a,p,q);

       for (i=0;i<num;i++)
          {
            printf("{");
              for (j=0; j< num;j++)
                 printf("%f,   ", a[p[i]][j]);
            printf("}\n");
          }

       for (i=0;i< num;i++)
          printf("p[i]  %d     q[i]  %d\n",p[i],q[i]);
}
*/


