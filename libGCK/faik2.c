#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* --- PROTOTIPOS DE FUNCIONES --- */
int octan(float ax, float ay);
int busca(int nvec, int ntot, int *oct, float *dist, int *neig);
int busca2(int nvec, int ntot, int *oct, float *dist, int *neig);
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

/* Modificado 'main()' a 'int main()' para cumplir con los estándares */
int main() {
  char file[50];
  FILE *in, *faik, *faike;
  float a, b, c, d, zestim, sigma, xx, yy, xmin, xmax, ymin, ymax, *x, *y, *z,
      *dist, *mat, *vec, *vec1, zk[5];
  int i, ix, iy, nvec, bins, hist, ntot = 0, nmon[] = {1, 3, 6}, nork, *neig,
                                   *oct;

  printf("FAIK2 PROGRAM: Generalized Covariance Universal Kriging\n\r");
  printf("\n\r\r Data filename:->");
  scanf("%s", file);

  in = fopen(file, "r");
  if (in == NULL) {
    printf("\n\rERROR! opening the file\n\r");
    abort();
  };
  fscanf(in, "%f %f %f %f", &a, &b, &c, &d);
  xmin = xmax = c;
  ymin = ymax = d;
  while ((fscanf(in, "%f %f %f %f", &a, &b, &c, &d)) > 0) {
    ntot += 1;
    xmin = (xmin > c) ? c : xmin;
    xmax = (xmax < c) ? c : xmax;
    ymin = (ymin > d) ? d : ymin;
    ymax = (ymax < d) ? d : ymax;
  };
  printf("\n\r%s contiene %d datos. \n\r(%g<X<%g) (%g<Y<%g)", file, ntot, xmin,
         xmax, ymin, ymax);
  fclose(in);
  in = fopen(file, "r");

  /*allocating memory*/

  if ((x = calloc((size_t)ntot, sizeof(float))) == NULL) {
    printf("\n\rERROR! Not enough memory\n\r");
    abort();
  };
  if ((y = calloc((size_t)ntot, sizeof(float))) == NULL) {
    printf("\n\rERROR! Not enough memory\n\r");
    abort();
  };
  if ((z = calloc((size_t)ntot, sizeof(float))) == NULL) {
    printf("\n\rERROR! Not enough memory\n\r");
    abort();
  };
  if ((dist = calloc((size_t)ntot, sizeof(float))) == NULL) {
    printf("\n\rERROR! Not enough memory\n\r");
    abort();
  };
  if ((oct = calloc((size_t)ntot, sizeof(int))) == NULL) {
    printf("\n\rERROR! Not enough memory\n\r");
    abort();
  };
  if ((neig = calloc((size_t)ntot, sizeof(int))) == NULL) {
    printf("\n\rERROR! Not enough memory\n\r");
    abort();
  };

  for (i = 0; i < ntot; i++)
    fscanf(in, "%f %f %f %f", &a, z + i, x + i, y + i);

  fclose(in);

  do {
    printf("\n\r\n\rEnter the order k (0,1,2):->");
    scanf("%d", &nork);
  } while ((nork > 2) || (nork < 0));
  do {
    printf("\n\rEnter the Number of Neighbors (<=32) :->");
    scanf("%d", &nvec);
  } while (nvec > 32);
  printf("\n\r\n\rEnter model parameters\n\r");
  scanf("%f %f %f %f %f", zk, zk + 1, zk + 2, zk + 3, zk + 4);
  printf("\n\r\n\rX-Y coordinates of the region to be estimated");
  printf("\n\rEnter the coordinates of the bottom left corner :");
  scanf("%f %f", &xmin, &ymin);
  printf("\n\rEnter the coordinates of the upper right corner :");
  scanf("%f %f", &xmax, &ymax);
  printf("\n\r\n\rNumero de BINS e HIST :");
  scanf("%d %d", &bins, &hist);

  if ((mat = calloc((size_t)(nvec + nmon[nork]) * (nvec + nmon[nork]),
                    sizeof(float))) == NULL) {
    printf("\n\rERROR! Not enough memory\n\r");
    abort();
  };
  if ((vec = calloc((size_t)(nvec + nmon[nork]), sizeof(float))) == NULL) {
    printf("\n\rERROR! Not enough memory\n\r");
    abort();
  };
  if ((vec1 = calloc((size_t)(nvec + nmon[nork]), sizeof(float))) == NULL) {
    printf("\n\rERROR! Not enough memory\n\r");
    abort();
  };

  faik = fopen("faik.grd", "w");
  faike = fopen("faike.grd", "w");
  fprintf(faik, "DSAA");
  fprintf(faik, "\n%d %d", bins, hist);
  fprintf(faik, "\n%g %g", xmin, xmax);
  fprintf(faik, "\n%g %g", ymin, ymax);
  fprintf(faik, "\n1. 10.\n");
  fprintf(faike, "DSAA");
  fprintf(faike, "\n%d %d", bins, hist);
  fprintf(faike, "\n%g %g", xmin, xmax);
  fprintf(faike, "\n%g %g", ymin, ymax);
  fprintf(faike, "\n1. 10.\n");

  for (iy = 0; iy < hist; iy++)
    for (ix = 0; ix < bins; ix++) {
      xx = xmin + ix * (xmax - xmin) / (bins - 1);
      yy = ymin + iy * (ymax - ymin) / (hist - 1);
      for (i = 0; i < ntot; i++) {
        a = xx - x[i];
        b = yy - y[i];
        dist[i] = (float)sqrt(1. * (a * a + b * b));
        oct[i] = octan(a, b);
      };
      if (busca2(nvec + 1, ntot, oct, dist, neig) != 1) {
        krige(mat, vec, vec1, nvec, nork, dist, neig, x, y, z, &zestim, &sigma,
              xx, yy, zk);
      } else {
        zestim = 1.70141e38;
        sigma = zestim;
      };
      if ((ix % 10) == 9) {
        fprintf(faik, "%g\n", zestim);
        fprintf(faike, "%g\n", sigma);
      } else {
        fprintf(faik, "%g ", zestim);
        fprintf(faike, "%g ", sigma);
      };
      /*			     printf("\n%e %e %e
       * %e",xx,yy,zestim,sigma);*/
      /*			     if(zestim>0.)escribe_gki(nvec,vec,ix);*/

      printf("\r%5d %5d", ix + 1, iy + 1);
    };
  printf("\n");
  fclose(faik);
  fclose(faike);
}
