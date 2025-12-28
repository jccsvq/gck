#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ntot=1,nmon[]={1,3,6},nork,*neig,*oct;
double *x,*y,*z,*dist,*mat,*vec,*vec1,zk[5];
FILE *faik,*faike;



main()
{
char file[50];
FILE *in;
double a,b,c,d,zestim,sigma,xx,yy,xmin,xmax,ymin,ymax;
int i,ix,iy,nvec,bins,hist;
printf("PROGRAMA FAIK2: Krigeage universal por covarianza generalizada\n\r");
printf("\n\r\r Fichero de datos:->");
scanf("%s",file);

in=fopen(file,"r");
if(in==NULL){printf("\n\rERROR! al abrir el fichero\n\r");abort();};
fscanf(in,"%lf %lf %lf %lf",&a,&b,&c,&d);
xmin=xmax=c;ymin=ymax=d;
while((fscanf(in,"%lf %lf %lf %lf",&a,&b,&c,&d))>0){
     ntot+=1;
     xmin=(xmin>c)? c:xmin;
     xmax=(xmax<c)? c:xmax;
     ymin=(ymin>d)? d:ymin;
     ymax=(ymax<d)? d:ymax;
     };
printf("\n\r%s contiene %d datos. \n\r(%g<X<%g) (%g<Y<%g)",file,ntot,xmin,xmax,ymin,ymax);
fclose(in);in=fopen(file,"r");

/* tomamos memoria */

if((x=calloc((unsigned) ntot,sizeof(double)))==NULL)
	 {printf("\n\rERROR! No hay memoria\n\r");abort();};
if((y=calloc((unsigned) ntot,sizeof(double)))==NULL)
	 {printf("\n\rERROR! No hay memoria\n\r");abort();};
if((z=calloc((unsigned) ntot,sizeof(double)))==NULL)
	 {printf("\n\rERROR! No hay memoria\n\r");abort();};
if((dist=calloc((unsigned) ntot,sizeof(double)))==NULL)
	 {printf("\n\rERROR! No hay memoria\n\r");abort();};
if((oct=calloc((unsigned) ntot,sizeof(int)))==NULL)
	 {printf("\n\rERROR! No hay memoria\n\r");abort();};
if((neig=calloc((unsigned) ntot,sizeof(int)))==NULL)
	 {printf("\n\rERROR! No hay memoria\n\r");abort();};

for(i=0;i<ntot;i++)fscanf(in,"%lf %lf %lf %lf",&a,z+i,x+i,y+i);

fclose(in);

do{printf("\n\r\n\rEntre el orden k (0,1,2):->");
scanf("%d",&nork);}while((nork>2)||(nork<0));
do{printf("\n\rEntre el Numero de Vecinos (<=32) :->");
scanf("%d",&nvec);}while(nvec>32);
printf("\n\r\n\rEntre los parametros del modelo\n\r");
scanf("%lf %lf %lf %lf %lf",zk,zk+1,zk+2,zk+3,zk+4);
printf("\n\r\n\rCoordenadas X-Y de la region a estimar");
printf("\n\rEntre coordenadas de la primera esquina :");
scanf("%lf %lf",&xmin,&ymin);
printf("\n\rEntre coordenadas de la segunda esquina :");
scanf("%lf %lf",&xmax,&ymax);
printf("\n\r\n\rNumero de BINS e HIST :");
scanf("%d %d",&bins,&hist);

if((mat=calloc((unsigned) (nvec+nmon[nork])*(nvec+nmon[nork]),sizeof(double)))==NULL)
	 {printf("\n\rERROR! No hay memoria\n\r");abort();};
if((vec=calloc((unsigned) (nvec+nmon[nork]),sizeof(double)))==NULL)
	 {printf("\n\rERROR! No hay memoria\n\r");abort();};
if((vec1=calloc((unsigned) (nvec+nmon[nork]),sizeof(double)))==NULL)
	 {printf("\n\rERROR! No hay memoria\n\r");abort();};

faik=fopen("faik.grd","w");faike=fopen("faike.grd","w");
fprintf(faik,"DSAA");
fprintf(faik,"\n%d %d",bins,hist);
fprintf(faik,"\n%g %g",xmin,xmax);
fprintf(faik,"\n%g %g",ymin,ymax);
fprintf(faik,"\n1. 10.\n\r");
fprintf(faike,"DSAA");
fprintf(faike,"\n%d %d",bins,hist);
fprintf(faike,"\n%g %g",xmin,xmax);
fprintf(faike,"\n%g %g",ymin,ymax);
fprintf(faike,"\n1. 10.\n");

for(iy=0;iy<hist;iy++)for(ix=0;ix<bins;ix++){
	xx=xmin+ix*(xmax-xmin)/(bins-1);
	yy=ymin+iy*(ymax-ymin)/(hist-1);
	for(i=0;i<ntot;i++){
		a=xx-x[i];b=yy-y[i];
		dist[i]=(double)sqrt(1.*(a*a+b*b));
		oct[i]=octan(a,b);
			   };
		if(busca(nvec+1)!=1){
			     krige(nvec,&zestim,&sigma,xx,yy);}
			     else{zestim=1.70141e38;sigma=zestim;};
			     if((ix%10)==9){fprintf(faik,"%g\n",zestim);fprintf(faike,"%g\n",sigma);}
			      else{fprintf(faik,"%g ",zestim);fprintf(faike,"%g ",sigma);};
/*			     printf("\n%e %e %e %e",xx,yy,zestim,sigma);*/
/*			     if(zestim>0.)escribe_gki(nvec,vec,ix);*/

	printf("\r%5d %5d",ix+1,iy+1);		};
fclose(faik);fclose(faike);}

octan(ax,ay)
double ax,ay;
{
int oc;
if(ax<0.){
	if(ay<0.){
		if((-ax)>=(-ay)){oc=4;}else{oc=5;};
		}else{
		if((-ax)>=ay){oc=3;}else{oc=2;};
		};
	}else{
	if(ay<0.){
		if(ax>=(-ay)){oc=7;}else{oc=6;};
		}else{
		if(ax>=ay){oc=0;}else{oc=1;};
		};

	};
return(oc);
}


busca(nvec)
int nvec;
{
int i,j,k,ifrac=0;
double w;
for(i=0;i<ntot;i++)neig[i]=i;
for(i=0;i<(ntot-1);i++)for(j=(i+1);j<ntot;j++){
	if(dist[i]>dist[j]){
		w=dist[i];
		dist[i]=dist[j];
		dist[j]=w;
		k=neig[i];
		neig[i]=neig[j];
		neig[j]=k;
		k=oct[i];
		oct[i]=oct[j];
		oct[j]=k;
		};
	};
/*	quita(nvec+8);*/
ifrac=octavino(nvec);

return(ifrac);
}


quita(nvec)
int nvec;
{
int i;
for(i=0;i<(ntot-1);i++){
neig[i]=neig[i+1];dist[i]=dist[i+1];oct[i]=oct[i+1];
};
}

octavino(nvec)
int nvec;
{
int n[8],na[7][8],i,j,nnn,nonul,ifrac;
double b[7][8];
for(i=0;i<8;i++)n[i]=0;
for(i=0;i<7;i++)for(j=0;j<8;j++){na[i][j]=(-1);b[i][j]=0.;};
for(i=0;i<ntot;i++){
	if(n[oct[i]]<6){
	  n[oct[i]]+=1;
	  na[n[oct[i]]-1][oct[i]]=neig[i];
	  b[n[oct[i]]-1][oct[i]]=dist[i];
	  };
	};
nonul=0;
for(i=0;i<8;i++)if(n[i]!=0)nonul++;
if(nonul<=5)return(1);
nnn=(-1);
for(i=0;i<7;i++)for(j=0;j<8;j++)if(na[i][j]!=(-1)){
    nnn++;
    neig[nnn]=na[i][j];
    dist[nnn]=b[i][j];
    };
if(nnn>=nvec){ifrac=0;}else{ifrac=1;
			     };
return(ifrac);
}

gaussj(a,yz,ne,ns)
double *a,*yz;
int ne,ns;
 {
double *(*ap),val,*tp,*(*yp);
int i,j,k,piv;

if((ap= (double * *)  calloc((unsigned) ne,sizeof(double *)))==NULL)
{printf("\nERROR: Insuficiente memoria para gaussj");abort();};
for(i=0;i<ne;i++)ap[i]=a+i*ne;
if((yp= (double * *)  calloc((unsigned) ne,sizeof(double *)))==NULL)
{printf("\nERROR: Insuficiente memoria para gaussj");abort();};
for(i=0;i<ne;i++)yp[i]=yz+i*ns;

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

krige(nvec,zestim,sigma,xx,yy)
int nvec;
double *zestim,*sigma,xx,yy;
{
double *(*a),*b,*bb,c0,estruc_gki(),monomio();
int i,ii,j,jj,dim;
dim=nvec+nmon[nork];
if((a=  (double * *) calloc((unsigned) dim,sizeof(double *)))==NULL)
{printf("\nERROR: Insuficiente memoria para ecuaciones");abort();};
for(i=0;i<dim;i++)a[i]=mat+i*dim;
b=vec;bb=vec1;

/*establecemos las ecuaciones de krigea_gikge */

/* vector */

for(i=0;i<nvec;i++){
  b[i]=estruc_gki(dist[i]);
  bb[i]=b[i];
  };
for(i=0;i<nmon[nork];i++){
  b[nvec+i]=(-monomio(xx,yy,i));
  bb[nvec+i]=b[nvec+i];
  };

/* diagonal principal */

c0=estruc_gki(0.);
for(i=0;i<nvec;i++)a[i][i]=c0;

/* caja nula */

for(i=0;i<nmon[nork];i++)for(j=0;j<nmon[nork];j++)a[nvec+i][nvec+j]=0.;

/* parte simetrica de la caja */

for(i=0;i<(nvec-1);i++){
  ii=neig[i];
  for(j=i+1;j<nvec;j++){
    jj=neig[j];
    a[i][j]=estruc_gki((double)(sqrt(1.*((x[ii]-x[jj])*(x[ii]-x[jj])+(y[ii]-y[jj])*(y[ii]-y[jj])))));
    a[j][i]=a[i][j];
    };
  };

/* cajas laterales */

for(i=0;i<nvec;i++){
  ii=neig[i];
  for(j=0;j<nmon[nork];j++){
    a[i][nvec+j]=(-monomio(x[ii],y[ii],j));
    a[nvec+j][i]=a[i][nvec+j];
    };
  };

/*   resolucion de las ecuaciones  */

if((gaussj(mat,vec,dim,1))==1){
   (*zestim)=(*sigma)=0.;
   for(i=0;i<nvec;i++) *zestim+=b[i]*z[neig[i]];
   for(i=0;i<dim;i++) *sigma+=b[i]*bb[i];
   if((*sigma-=c0)>0.){*sigma=(double)sqrt(1.*(*sigma));}else{*sigma=9e-29;*zestim=(-*sigma);};
  }else{*sigma=9e-29;(*zestim)=(-*sigma);
				    };
free(a);}

double estruc_gki(h)
double h;
{
double est;
/*if(h>0.){est=zk[0]+h*zk[1]+h*h*h*zk[2]+h*h*h*h*h*zk[3]+h*h*(double)log((double)h);}*/
if(h>0.){est=h*zk[1]+h*h*h*zk[2]+h*h*h*h*h*zk[3]+h*h*zk[4]*(double)log((double)h);}
else{est=zk[0];};
return( (-est));
}

double monomio(ax,ay,i)
double ax,ay;
int i;
{
double mono;
switch(i){
  case 0 :{mono=1.;break;};
  case 1 :{mono=ax;break;};
  case 2 :{mono=ay;break;};
  case 3 :{mono=ax*ax;break;};
  case 4 :{mono=ay*ay;break;};
  case 5 :{mono=ax*ay;break;};
  };
return(mono);
}
