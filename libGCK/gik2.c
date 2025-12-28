#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* --- PROTOTIPOS DE FUNCIONES --- */
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
/* ------------------------------- */

int main() {
  char file[50];
  FILE *in, *gki, *cov;
  float a, b, c, d, zestim, sigma, xx, yy, *x, *y, *z, *dist, *mat, *vec, *vec1;
  int i, ix, nvec, ntot = 0, nmon[] = {1, 3, 6}, nork, *neig, *oct, ngik = 0;

  printf("GIK2 PROGRAM: Calculation of generalized increments of order k\n");
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

  gki = fopen("gki.dat", "w");
  cov = fopen("cov.dat", "w");

  for (ix = 0; ix < ntot; ix++) {
    xx = x[ix];
    yy = y[ix];
    for (i = 0; i < ntot; i++) {
      a = xx - x[i];
      b = yy - y[i];
      dist[i] = (float)sqrt(1. * (a * a + b * b));
      oct[i] = octan(a, b);
    };

    if (busca(nvec + 1, ntot, oct, dist, neig) != 1) {
      krigea_gik_1(mat, vec, vec1, nvec, nork, dist, neig, x, y, z, &zestim,
                   &sigma, xx, yy);
      if (zestim > 0.)
      /*FILE *gki, FILE *cov, int nvec, int nork, int *neig, float *x,
                 float *y, float *z, float *ponde, int ix*/
      {
        escribe_gki(gki, cov, nvec, nork, neig, x, y, z, vec, ix);
        ngik++;
      }
    };
    printf("\rix= %d", ix);
  };
  printf("\n\nGKI.DAT file created with %d generalized increments\n",
         ngik);
  fclose(gki);
  fclose(cov);
}
