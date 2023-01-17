#include "cstd.h"
#include <complex.h>

int main()
{
int n,i,j,iter;
double alpha,beta,rsold,rsnew,pAp;
float **A,*x,*b,*r,*p,*Ap;

n= 5;

A=alloc2float(n,n);
x=alloc1float(n);
b=alloc1float(n);
r=alloc1float(n);
p=alloc1float(n);
Ap=alloc1float(n);
 
for(j=0;j<n;j++) x[j]=j;
for(j=0;j<n;j++){
   b[j]=0;
   for(i=0;i<n;i++) {
   A[j][i]=fabs(i-j);
   b[j]+=A[j][i]*x[i];
   }
}

memset(x,0,n*sizeof(float));
for (rsold=0,i=0;i<n;i++){
   r[i]=b[i];
   p[i]=r[i];
   rsold+=r[i]*r[i];
   }

for(iter=0;iter<100;iter++){
   for(j=0;j<n;j++){
     Ap[j]=0;
     for(i=0;i<n;i++){
      Ap[j]+=A[j][i]*p[i];
     }
   }
  pAp=0;
  for(i=0;i<n;i++){
    pAp+=p[i]*Ap[i];
  }
    alpha=rsold/pAp;
    rsnew=0;
  for(i=0;i<n;i++){
    x[i]+=alpha*p[i];
    r[i]-=alpha*Ap[i];
    rsnew+=r[i]*r[i];
  }
  if(rsnew<1e-7)break;
  beta=rsnew/rsold;
  for(i=0;i<n;i++) p[i]=r[i]+beta*p[i];
  rsold=rsnew;
}
for(i=0;i<n;i++) printf("x=%e\n",x[i]);



free2float(A);
free1float(x);
free1float(b);
free1float(r);
free1float(p);
free1float(Ap);
return 0;

}
