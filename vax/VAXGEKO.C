#include <stdio.h>
#include <stdlib.h>
#include <math.h>




main()
{
FILE *cov,*mod;
int nlam=0,i,j,iq[5];
float xx[5],yy,w,qq[5][5],aa[5],q[5][5],a[5],z0=0.0,q0=0.0;

printf("Programa GEKO: Calculo de covarianzas generalizadas\n\n");
for(i=0;i<5;i++){aa[i]=0;for(j=0;j<5;j++)qq[i][j]=0;};
if((cov=fopen("cov.dat","r"))==NULL)
       printf("\nERROR no se puede abrir el fichero COV.DAT");
while((fscanf(cov,"%f %f %f %f %f %f %f ",xx,xx+1,xx+2,xx+3,xx+4,&yy,&w))>2)
  {
  nlam++;
  z0+=yy;
  q0+=w*yy*yy;
  for(i=0;i<5;i++){
    aa[i]+=yy*xx[i]*w;
    for(j=i;j<5;j++){
      qq[i][j]+=xx[i]*xx[j]*w;
      if(i!=j)qq[j][i]=qq[i][j];
      };
    };
  printf("\rProcesando GIK n§:%4d",nlam);
  };
/*for(i=0;i<5;i++)printf("\n\n %e %e %e %e %e : %e ",qq[i][0],qq[i][1],qq[i][2],qq[i][3],qq[i][4],aa[i]);*/
fclose(cov);
if((cov=fopen("modelos.dat","r"))==NULL)
       printf("\nERROR no se puede abrir el fichero MODELOS.DAT");
mod=fopen("modcov.dat","w");
while((fscanf(cov,"%d %d %d %d %d",iq,iq+1,iq+2,iq+3,iq+4))>2){
  for(i=0;i<5;i++){
     a[i]=iq[i]*aa[i];
     for(j=i;j<5;j++){
       q[i][j]=qq[i][j]*iq[i]*iq[j];
       q[j][i]=q[i][j];
       };
     if(iq[i]==0)q[i][i]=1.;
     };
  if((gaussj(q[0],a,5,1))==1){
      if(  ((a[4]==0.)&&
	    (a[0]>=0.)&&
	    (a[1]<=0.)&&
	    (a[3]<=0.)&&
	    (a[2]>=(float)(-10.*sqrt(1.*a[1]*a[3])/3.)) )||
	   ((a[3]==0.)&&
	    (a[0]>=0.)&&
	    (a[1]<=0.)&&
	    (a[2]>=0.)&&
	    (a[4]>=(float)(-1.5*sqrt(-1.*a[2]*a[1]))) )  )
	  {fprintf(mod,"\n%e %e %e %e %e",a[0],a[1],a[2],a[3],a[4]);};
    };
  };



/*fcloseall();*/
}

gaussj(a,y,ne,ns)
float *a,*y;
int ne,ns;
 {
float *(*ap),val,*tp,*(*yp);
int i,j,k,piv;

if((ap=   calloc(ne,sizeof(float *)))==NULL)
{printf("\nERROR: Insuficiente memoria para gaussj");abort();};
for(i=0;i<ne;i++)ap[i]=a+i*ne;
if((yp=   calloc(ne,sizeof(float *)))==NULL)
{printf("\nERROR: Insuficiente memoria para gaussj");abort();};
for(i=0;i<ne;i++)yp[i]=y+i*ns;

for(i=0;i<ne;i++)
/*para cada fila*/
{
/*buscar pivote maximo*/
piv=i;val=fabs(ap[i][i]);
for(j=i;j<ne;j++){if((fabs(ap[j][i]))>val){val=fabs(ap[j][i]);piv=j;}; };
if(val==0.){free(ap);free(yp);return(0);};

/*permutamos las filas*/
tp=ap[i];ap[i]=ap[piv];ap[piv]=tp;
tp=yp[i];yp[i]=yp[piv];yp[piv]=tp;
val=ap[i][i];
/*normalizamos ecuacion i*/
for(j=i;j<ne;j++)ap[i][j]/=val;
for(j=0;j<ns;j++)yp[i][j]/=val;
/*reducimos las demas*/
for(k=0;k<ne;k++)if(k!=i){
     val=ap[k][i];
     for(j=i;j<ne;j++){ap[k][j]-=val*ap[i][j];};
     for(j=0;j<ns;j++){yp[k][j]-=val*yp[i][j];};

			 };
		};

for(i=0;i<ne;i++){
for(j=0;j<ns;j++){ap[i][j]=yp[i][j];};
};
for(i=0;i<ne;i++){
for(j=0;j<ns;j++){*(y+i*ns+j)=ap[i][j];};
};


free(ap);free(yp);
return(1);}

