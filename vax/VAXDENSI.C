#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ntot=1,nmon[]={1,3,6},nork,*neig,*oct;
float *x,*y;
FILE *faik;



main()
{
char file[50];
FILE *in;
float a,b,zestim,xx,yy,xmin,xmax,ymin,ymax,h1,h2;
int i,ix,iy,nvec,bins,hist;

printf("Programa DENSI: Calculo de densidad por kernels\n");
printf("\n Fichero de datos:->");
scanf("%s",file);

in=fopen(file,"r");
if(in==NULL){printf("\nERROR! al abrir el fichero\n");abort();};
fscanf(in,"%f %f ",&a,&b);
xmin=xmax=a;ymin=ymax=b;
while((fscanf(in,"%f %f ",&a,&b))>0){
     ntot+=1;
     xmin=(xmin>a)? a:xmin;
     xmax=(xmax<a)? a:xmax;
     ymin=(ymin>b)? b:ymin;
     ymax=(ymax<b)? b:ymax;
     };
printf("\n%s contiene %d datos. \n(%g<X<%g) (%g<Y<%g)",file,ntot,xmin,xmax,ymin,ymax);
fclose(in);in=fopen(file,"r");

h1=.1*(xmax-xmin)/sqrt(1.*ntot);
h2=.1*(ymax-ymin)/sqrt(1.*ntot);

h1*=2;h2*=2;

/* tomamos memoria */

if((x=calloc((unsigned) ntot,sizeof(float)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};
if((y=calloc((unsigned) ntot,sizeof(float)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};

for(i=0;i<ntot;i++)fscanf(in,"%f %f",x+i,y+i);

fclose(in);

printf("\n\nCoordenadas X-Y de la region a estimar");
printf("\nEntre coordenadas de la primera esquina :");
scanf("%f %f",&xmin,&ymin);
printf("\nEntre coordenadas de la segunda esquina :");
scanf("%f %f",&xmax,&ymax);
printf("\n\nNumero de BINS e HIST :");
scanf("%d %d",&bins,&hist);



faik=fopen("densi.grd","w");
/*faike=fopen("faike.grd","w");*/
fprintf(faik,"DSAA");
fprintf(faik,"\n%d %d",bins,hist);
fprintf(faik,"\n%g %g",xmin,xmax);
fprintf(faik,"\n%g %g",ymin,ymax);
fprintf(faik,"\n1. 10.\n");

for(iy=0;iy<hist;iy++)for(ix=0;ix<bins;ix++){
	xx=xmin+ix*(xmax-xmin)/(bins-1);
	yy=ymin+iy*(ymax-ymin)/(hist-1);
	zestim=0.;
	for(i=0;i<ntot;i++){
	     zestim+=(float)(exp(-.5*(pow((xx-x[i])/h1,2.)
			+pow((yy-y[i])/h2,2.))));
			   };
     if((ix%10)==9){fprintf(faik,"%g\n",zestim);}
      else{fprintf(faik,"%g ",zestim);};

	printf("\r%5d %5d",ix+1,iy+1);		};
fclose(faik);}

