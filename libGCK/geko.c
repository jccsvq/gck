#include <math.h>
#include <stdio.h>
#include <stdlib.h> // Añadido para 'calloc' y para estándares

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

/* Modificado 'main()' a 'int main()' para cumplir con los estándares */
int main() {
  FILE *cov, *mod;
  int nlam = 0, i, j, iq[5];
  float xx[5], yy, w, qq[5][5], aa[5], q[5][5], a[5], z0 = 0.0, q0 = 0.0;
  
  // ...
  
  printf("GEKO PROGRAM: Calculation of Generalized Covariances\n\n");
  for (i = 0; i < 5; i++) {
    aa[i] = 0;
    for (j = 0; j < 5; j++)
      qq[i][j] = 0;
  };
  if ((cov = fopen("cov.dat", "r")) == NULL)
    printf("\nERROR: Cannot open file COV.DAT");
  while ((fscanf(cov, "%f %f %f %f %f %f %f ", xx, xx + 1, xx + 2, xx + 3,
                 xx + 4, &yy, &w)) > 2) {
    nlam++;
    z0 += yy;
    q0 += w * yy * yy;
    for (i = 0; i < 5; i++) {
      aa[i] += yy * xx[i] * w;
      for (j = i; j < 5; j++) {
        qq[i][j] += xx[i] * xx[j] * w;
        if (i != j)
          qq[j][i] = qq[i][j];
      };
    };
    printf("\rProcessing GIK #:%4d", nlam);
  };
  /*for(i=0;i<5;i++)printf("\n\n %e %e %e %e %e : %e
   * ",qq[i][0],qq[i][1],qq[i][2],qq[i][3],qq[i][4],aa[i]);*/
  fclose(cov);
  printf("\n");
  if ((cov = fopen("modelos.dat", "r")) == NULL)
    printf("\nERROR: The file MODELOS.DAT cannot be opened.");
  mod = fopen("modcov.dat", "w");
  while ((fscanf(cov, "%d %d %d %d %d", iq, iq + 1, iq + 2, iq + 3, iq + 4)) >
         2) {
    for (i = 0; i < 5; i++) {
      a[i] = iq[i] * aa[i];
      for (j = i; j < 5; j++) {
        q[i][j] = qq[i][j] * iq[i] * iq[j];
        q[j][i] = q[i][j];
      };
      if (iq[i] == 0)
        q[i][i] = 1.;
    };
    if ((gaussj(q[0], a, 5, 1)) == 1) {
      if (((a[4] == 0.) && (a[0] >= 0.) && (a[1] <= 0.) && (a[3] <= 0.) &&
           (a[2] >= (float)(-10. * sqrt(1. * a[1] * a[3]) / 3.))) ||
          ((a[3] == 0.) && (a[0] >= 0.) && (a[1] <= 0.) && (a[2] >= 0.) &&
           (a[4] >= (float)(-1.5 * sqrt(-1. * a[2] * a[1]))))) {
        fprintf(mod, "\n%e %e %e %e %e", a[0], a[1], a[2], a[3], a[4]);
      };
    };
  };
  fclose(cov);
  fclose(mod);
}

