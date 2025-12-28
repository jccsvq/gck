#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

int octan(float ax, float ay);
int busca(int nvec, int ntot, int *oct, float *dist, int *neig);
void quita(int ntot, int *neig, float *dist, int *oct);
int octavino(int nvec, int ntot, int *oct, float *dist, int *neig);
int gaussj(float *a, float *y, int ne, int ns);
float monomio(float ax, float ay, int i);
float estruc_gki(float h, float zk[5]);
float covmodel(float h, int i);
void krige(float *mat, float *vec, float *vec1, int nvec, int nork, float *dist,
           int *neig, float *x, float *y, float *z, float *zestim, float *sigma,
           float xx, float yy, float zk[5]);
void krigea_gik_1(float *mat, float *vec, float *vec1, int nvec, int nork,
                  float *dist, int *neig, float *x, float *y, float *z,
                  float *zestim, float *sigma, float xx, float yy);
float estruc_gki_1(float h);
float covmodel_1(float h, int i);
void escribe_gki(FILE *gki, FILE *cov, int nvec, int nork, int *neig, float *x,
                 float *y, float *z, float *ponde, int ix);

int main() {

  char file[50];
  FILE *in,*mod,*out;
  float a, b, c, d, zestim, sigma, xx, yy, sumcuad,*x,*y,*z,*dist,*mat,*vec,*vec1, zk[5];
  int i, ix, nvec, npoint,ntot,*oct,*neig,nork, nmon[] = {1, 3, 6};

  
  ntot = 0;
  printf("PROGRAMA CROSVAL: validacion cruzada de modelos\n");
  printf("\n Data filename:->");
  scanf("%s", file);

  in = fopen(file, "r");
  if (in == NULL) {
    printf("\nERROR! opening the file\n");
    abort();
  };
  while ((fscanf(in, "%f %f %f %f", &a, &b, &c, &d)) > 0) {
    ntot += 1;
  };
  printf("\n%s contiene %d datos.", file, ntot);
  fclose(in);
  in = fopen(file, "r");

  /*allocating memory*/

  if ((x = calloc((size_t)ntot, sizeof(float))) == NULL) {
    printf("\nERROR! Not enough memory\n");
    abort();
  };
  if ((y = calloc((size_t)ntot, sizeof(float))) == NULL) {
    printf("\nERROR! Not enough memory\n");
    abort();
  };
  if ((z = calloc((size_t)ntot, sizeof(float))) == NULL) {
    printf("\nERROR! Not enough memory\n");
    abort();
  };
  if ((dist = calloc((size_t)ntot, sizeof(float))) == NULL) {
    printf("\nERROR! Not enough memory\n");
    abort();
  };
  if ((oct = calloc((size_t)ntot, sizeof(int))) == NULL) {
    printf("\nERROR! Not enough memory\n");
    abort();
  };
  if ((neig = calloc((size_t)ntot, sizeof(int))) == NULL) {
    printf("\nERROR! Not enough memory\n");
    abort();
  };

  for (i = 0; i < ntot; i++)
    fscanf(in, "%f %f %f %f", &a, z + i, x + i, y + i);

  fclose(in);

  do {
    printf("\n\nEnter the order k (0,1,2):->");
    scanf("%d", &nork);
  } while ((nork > 2) || (nork < 0));
  do {
    printf("\nEnter the Number of Neighbors (<=32) :->");
    scanf("%d", &nvec);
  } while (nvec > 32);

  if ((mat = calloc((size_t)(nvec + nmon[nork]) * (nvec + nmon[nork]),
                    sizeof(float))) == NULL) {
    printf("\nERROR! Not enough memory\n");
    abort();
  };
  if ((vec = calloc((size_t)(nvec + nmon[nork]), sizeof(float))) == NULL) {
    printf("\nERROR! Not enough memory\n");
    abort();
  };
  if ((vec1 = calloc((size_t)(nvec + nmon[nork]), sizeof(float))) == NULL) {
    printf("\nERROR! Not enough memory\n");
    abort();
  };

  printf("\n Fitted models file:->");
  scanf("%s", file);

  mod = fopen(file, "r");
  if (mod == NULL) {
    printf("\nERROR! opening the file\n");
    abort();
  };
  out = fopen("cv.dat", "w");

  while ((fscanf(mod, "%f %f %f %f %f", zk, zk + 1, zk + 2, zk + 3, zk + 4)) >
         2) {
    sumcuad = 0.;
    npoint = 0;
    printf("\n");
    for (ix = 0; ix < ntot; ix++) {
      printf("\r ix=%d", ix);
      xx = x[ix];
      yy = y[ix];
      for (i = 0; i < ntot; i++) {
        a = xx - x[i];
        b = yy - y[i];
        dist[i] = (float)sqrt(1. * (a * a + b * b));
        oct[i] = octan(a, b);
      };
      /*int nvec, int ntot, int *oct, float *dist, int *neig*/
      if (busca(nvec + 1,ntot,oct,dist,neig) != 1) {
        krige(mat, vec, vec1, nvec, nork, dist, neig, x, y, z, &zestim, &sigma,
              xx, yy,zk);
        if (zestim > 0.) {
          sumcuad += (z[ix] - zestim) * (z[ix] - zestim);
          npoint++;
        };
      };
    };
    fprintf(out, "%g %d\n", sumcuad / npoint, npoint);
    fprintf(out, "%g %g %g %g %g\n", zk[0], zk[1], zk[2], zk[3], zk[4]);
    printf("%g %d\n", sumcuad / npoint, npoint);
    printf("%g %g %g %g %g\n", zk[0], zk[1], zk[2], zk[3], zk[4]);
  };
  printf("\n");
  fclose(mod);
  fclose(out);
}
