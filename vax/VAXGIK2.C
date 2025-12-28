#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ntot=0,nmon[]={1,3,6},nork,*neig,*oct,ngik=0;
float *x,*y,*z,*dist,*mat,*vec,*vec1;
FILE *gki,*cov;



main()
{
char file[50];
FILE *in;
float a,b,c,d,zestim,sigma,xx,yy;
int i,ix,nvec;

printf("PROGRAMA GIK2: Calculo de incrementos generalizados de orden k\n");
printf("\n Fichero de datos:->");
scanf("%s",&file);

in=fopen(file,"r");
if(in==NULL){printf("\nERROR! al abrir el fichero\n");abort();};
while((fscanf(in,"%f %f %f %f",&a,&b,&c,&d))>0){ntot+=1;};
printf("\n%s contiene %d datos.",file,ntot);
fclose(in);in=fopen(file,"r");

/* tomamos memoria */

if((x=calloc((size_t) ntot,sizeof(float)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};
if((y=calloc((size_t) ntot,sizeof(float)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};
if((z=calloc((size_t) ntot,sizeof(float)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};
if((dist=calloc((size_t) ntot,sizeof(float)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};
if((oct=calloc((size_t) ntot,sizeof(int)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};
if((neig=calloc((size_t) ntot,sizeof(int)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};

for(i=0;i<ntot;i++)fscanf(in,"%f %f %f %f",&a,z+i,x+i,y+i);

fclose(in);

do{printf("\n\nEntre el orden k (0,1,2):->");
scanf("%d",&nork);}while((nork>2)||(nork<0));
do{printf("\nEntre el Numero de Vecinos (<=32) :->");
scanf("%d",&nvec);}while(nvec>32);

if((mat=calloc((size_t) (nvec+nmon[nork])*(nvec+nmon[nork]),sizeof(float)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};
if((vec=calloc((size_t) (nvec+nmon[nork]),sizeof(float)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};
if((vec1=calloc((size_t) (nvec+nmon[nork]),sizeof(float)))==NULL)
	 {printf("\nERROR! No hay memoria\n");abort();};

gki=fopen("gki.dat","w");
cov=fopen("cov.dat","w");

for(ix=0;ix<ntot;ix++){
	xx=x[ix];yy=y[ix];
	for(i=0;i<ntot;i++){
		a=xx-x[i];b=yy-y[i];
		dist[i]=(float)sqrt(1.*(a*a+b*b));
		oct[i]=octan(a,b);
			   };
		if(busca(nvec+1)!=1){
			     krigea_gik(nvec,&zestim,&sigma,xx,yy);
			     if(zestim>0.)escribe_gki(nvec,vec,ix);
				};
	printf("\rix= %d",ix);		};
printf("\n\nCreado el fichero GKI.DAT con %d incrementos generalizados",ngik);
fclose(gki);fclose(cov);}

octan(ax,ay)
float ax,ay;
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
float w;
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
	quita(nvec+8);
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
float b[7][8];
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

krigea_gik(nvec,zestim,sigma,xx,yy)
int nvec;
float *zestim,*sigma,xx,yy;
{
float *(*a),*b,*bb,c0,estruc_gki(),monomio();
int i,ii,j,jj,dim;
dim=nvec+nmon[nork];
if((a=   calloc(dim,sizeof(float *)))==NULL)
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
    a[i][j]=estruc_gki((float)(sqrt(1.*((x[ii]-x[jj])*(x[ii]-x[jj])+(y[ii]-y[jj])*(y[ii]-y[jj])))));
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
   if((*sigma-=c0)>0.){*sigma=(float)sqrt(1.*(*sigma));}else{*sigma=9e-29;(*zestim)=(-*sigma);};
  }else{*sigma=9e-29;*zestim=(-*sigma);
				    };
free(a);}

float estruc_gki(h)
float h;
{

return((float) h);
}

float monomio(ax,ay,i)
float ax,ay;
int i;
{
float mono;
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
float covmodel(h,i)
float h;
int i;
{
float mono;
switch(i){
  case 0 :{mono=h;break;};
  case 1 :{mono=h*h*h;break;};
  case 2 :{mono=h*h*h*h*h;break;};
  case 3 :{mono=(h<=0.)? 0.:(float)h*h*log((double)h);break;};
  };
return(mono);
}

escribe_gki(nvec,ponde,ix)
int nvec,ix;
float *ponde;
{
int i,j,k;
float sq,dx,dy,cx[4],cy,cw;
sq=1.;
for(i=nvec;i>=1;i--){
  ponde[i]=(-ponde[i-1]);
  sq+=ponde[i]*ponde[i];
  };
ponde[0]=1.;
sq=(float)sqrt(1.*sq);
for(i=0;i<=nvec;i++)ponde[i]/=sq;
fprintf(gki,"%d %d\n%d ",nork,nvec+1,ix);
for(i=0;i<nvec;i++)fprintf(gki,"%d ",neig[i]);
fprintf(gki,"\n");
for(i=0;i<=nvec;i++)fprintf(gki,"%f ",ponde[i]);
fprintf(gki,"\n");
printf("      N§ de GIK's:%d       ",++ngik);
cy=ponde[0]*z[ix];
for(i=1;i<=nvec;i++)cy+=ponde[i]*z[neig[i-1]];cy=cy*cy;
for(i=nvec;i>0;i--)neig[i]=neig[i-1];
neig[0]=ix;
for(k=0;k<4;k++) cx[k]=0;
     for(i=0;i<=nvec;i++)for(j=0;j<=nvec;j++){
		dx=x[neig[i]]-x[neig[j]];
		dy=y[neig[i]]-y[neig[j]];
		cw=(float)sqrt(1.*(dx*dx+dy*dy));
		for(k=0;k<4;k++)cx[k]+=ponde[i]*ponde[j]*covmodel(cw,k);
			};
cw=1./(cx[1]*cx[1]);
fprintf(cov,"1. %e %e %e %e %e %e\n",cx[0],cx[1],cx[2],cx[3],cy,cw);
}
