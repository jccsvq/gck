#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* --- PROTOTIPOS DE FUNCIONES --- */
int octan(float ax, float ay);
int busca(int nvec);
void quita();
int octavino(int nvec);
int gaussj(float *a, float *y, int ne, int ns);
void krige(int nvec, float *zestim, float *sigma, float xx, float yy);
float estruc_gki(float h);
float monomio(float ax, float ay, int i);
/* ------------------------------- */

int ntot = 1, nmon[] = {1, 3, 6}, nork, *neig, *oct;
float *x, *y, *z, *dist, *mat, *vec, *vec1, zk[5];
FILE *faik, *faike;

/* Modificado 'main()' a 'int main()' para cumplir con los estándares */
int main() {
  char file[50];
  FILE *in;
  float a, b, c, d, zestim, sigma, xx, yy, xmin, xmax, ymin, ymax;
  int i, ix, iy, nvec, bins, hist;
  
  // ...
  
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
  printf("\n\r%s has %d data rows. \n\r(%g<X<%g) (%g<Y<%g)", file, ntot, xmin,
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
      if (busca(nvec + 1) != 1) {
        krige(nvec, &zestim, &sigma, xx, yy);
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

int octan(float ax, float ay)
{
  int oc;
  
  if (ax < 0.) {
    if (ay < 0.) {
      if ((-ax) >= (-ay)) {
        oc = 4;
      } else {
        oc = 5;
      };
    } else {
      if ((-ax) >= ay) {
        oc = 3;
      } else {
        oc = 2;
      };
    };
  } else {
    if (ay < 0.) {
      if (ax >= (-ay)) {
        oc = 7;
      } else {
        oc = 6;
      };
    } else {
      if (ax >= ay) {
        oc = 0;
      } else {
        oc = 1;
      };
    };
  };
  return (oc);
}

int busca(int nvec)
{
  int i, j, k, ifrac = 0;
  float w;
  
  for (i = 0; i < ntot; i++)
    neig[i] = i;
  for (i = 0; i < (ntot - 1); i++)
    for (j = (i + 1); j < ntot; j++) {
      if (dist[i] > dist[j]) {
        w = dist[i];
        dist[i] = dist[j];
        dist[j] = w;
        k = neig[i];
        neig[i] = neig[j];
        neig[j] = k;
        k = oct[i];
        oct[i] = oct[j];
        oct[j] = k;
      };
    };
  /*	quita(nvec+8);*/
  ifrac = octavino(nvec);

  return (ifrac);
}

void quita()
{
  int i;
  for (i = 0; i < (ntot - 1); i++) {
    neig[i] = neig[i + 1];
    dist[i] = dist[i + 1];
    oct[i] = oct[i + 1];
  };
}

int octavino(int nvec)
{
  int n[8], na[7][8], i, j, nnn, nonul, ifrac;
  float b[7][8];
  
  for (i = 0; i < 8; i++)
    n[i] = 0;
  for (i = 0; i < 7; i++)
    for (j = 0; j < 8; j++) {
      na[i][j] = -1;
      b[i][j] = 0.;
    };
  for (i = 0; i < ntot; i++) {
    if (n[oct[i]] < 6) {
      n[oct[i]] += 1;
      na[n[oct[i]] - 1][oct[i]] = neig[i];
      b[n[oct[i]] - 1][oct[i]] = dist[i];
    };
  };
  nonul = 0;
  for (i = 0; i < 8; i++)
    if (n[i] != 0)
      nonul++;
  if (nonul <= 5)
    return (1);
  nnn = -1;
  for (i = 0; i < 7; i++)
    for (j = 0; j < 8; j++)
      if (na[i][j] != -1) {
        nnn++;
        neig[nnn] = na[i][j];
        dist[nnn] = b[i][j];
      };
  if (nnn >= nvec) {
    ifrac = 0;
  } else {
    ifrac = 1;
  };
  return (ifrac);
}

int gaussj(float *a, float *y, int ne, int ns)
{
  float *(*ap), val, *tp, *(*yp);
  int i, j, k, piv;
  
  if ((ap = calloc(ne, sizeof(float *))) == NULL) {
    printf("\nERROR: Insufficient memory for gaussj!");
    abort();
  };
  for (i = 0; i < ne; i++)
    ap[i] = a + i * ne;
  if ((yp = calloc(ne, sizeof(float *))) == NULL) {
    printf("\nERROR: Insufficient memory for gaussj!");
    abort();
  };
  for (i = 0; i < ne; i++)
    yp[i] = y + i * ns;

  for (i = 0; i < ne; i++)
  /*para cada fila*/
  {
    /*buscar pivote maximo*/
    piv = i;
    val = fabs(ap[i][i]);
    for (j = i; j < ne; j++) {
      if ((fabs(ap[j][i])) > val) {
        val = fabs(ap[j][i]);
        piv = j;
      };
    };
    if (val == 0.) {
      free(ap);
      free(yp);
      return (0);
    };

    /*permutamos las filas*/
    tp = ap[i];
    ap[i] = ap[piv];
    ap[piv] = tp;
    tp = yp[i];
    yp[i] = yp[piv];
    yp[piv] = tp;
    val = ap[i][i];
    /*normalizamos ecuacion i*/
    for (j = i; j < ne; j++)
      ap[i][j] /= val;
    for (j = 0; j < ns; j++)
      yp[i][j] /= val;
    /*reducimos las demas*/
    for (k = 0; k < ne; k++)
      if (k != i) {
        val = ap[k][i];
        for (j = i; j < ne; j++) {
          ap[k][j] -= val * ap[i][j];
        };
        for (j = 0; j < ns; j++) {
          yp[k][j] -= val * yp[i][j];
        };
      };
  };

  for (i = 0; i < ne; i++) {
    for (j = 0; j < ns; j++) {
      ap[i][j] = yp[i][j];
    };
  };
  for (i = 0; i < ne; i++) {
    for (j = 0; j < ns; j++) {
      *(y + i * ns + j) = ap[i][j];
    };
  };

  free(ap);
  free(yp);
  return (1);
}

void krige(int nvec, float *zestim, float *sigma, float xx,
           float yy)
{
  float *(*a), *b, *bb, c0;
  int i, ii, j, jj, dim;
  
  dim = nvec + nmon[nork];
  if ((a = calloc(dim, sizeof(float *))) == NULL) {
    printf("\nERROR: Insufficient memory for equations!");
    abort();
  };
  for (i = 0; i < dim; i++)
    a[i] = mat + i * dim;
  b = vec;
  bb = vec1;

  /*establecemos las ecuaciones de krigea_gikge */

  /* vector */

  for (i = 0; i < nvec; i++) {
    b[i] = estruc_gki(dist[i]);
    bb[i] = b[i];
  };
  for (i = 0; i < nmon[nork]; i++) {
    b[nvec + i] = -monomio(xx, yy, i);
    bb[nvec + i] = b[nvec + i];
  };

  /* diagonal principal */

  c0 = estruc_gki(0.);
  for (i = 0; i < nvec; i++)
    a[i][i] = c0;

  /* caja nula */

  for (i = 0; i < nmon[nork]; i++)
    for (j = 0; j < nmon[nork]; j++)
      a[nvec + i][nvec + j] = 0.;

  /* parte simetrica de la caja */

  for (i = 0; i < (nvec - 1); i++) {
    ii = neig[i];
    for (j = i + 1; j < nvec; j++) {
      jj = neig[j];
      a[i][j] =
          estruc_gki((float)(sqrt(1. * ((x[ii] - x[jj]) * (x[ii] - x[jj]) +
                                        (y[ii] - y[jj]) * (y[ii] - y[jj])))));
      a[j][i] = a[i][j];
    };
  };

  /* cajas laterales */

  for (i = 0; i < nvec; i++) {
    ii = neig[i];
    for (j = 0; j < nmon[nork]; j++) {
      a[i][nvec + j] = -monomio(x[ii], y[ii], j);
      a[nvec + j][i] = a[i][nvec + j];
    };
  };

  /* resolucion de las ecuaciones  */

  if ((gaussj(mat, vec, dim, 1)) == 1) {
    *zestim = *sigma = 0.;
    for (i = 0; i < nvec; i++)
      *zestim += b[i] * z[neig[i]];
    for (i = 0; i < dim; i++)
      *sigma += b[i] * bb[i];
    if ((*sigma -= c0) > 0.) {
      *sigma = (float)sqrt(1. * (*sigma));
    } else {
      *sigma = 9e-29;
      *zestim = -*sigma;
    };
  } else {
    *sigma = 9e-29;
    *zestim = -*sigma;
  };

  free(a);
}

float estruc_gki(float h)
{
  float est;
  /*if(h>0.){est=zk[0]+h*zk[1]+h*h*h*zk[2]+h*h*h*h*h*zk[3]+h*h*(float)log((double)h);}*/
  if (h > 0.) {
    est = h * zk[1] + h * h * h * zk[2] + h * h * h * h * h * zk[3] +
          h * h * zk[4] * (float)log((double)h);
  } else {
    est = zk[0];
  };
  return ((-est));
}

float monomio(float ax, float ay, int i)
{
  float mono;
  mono = 0.0;
  switch (i) {
  case 0: {
    mono = 1.;
    break;
  };
  case 1: {
    mono = ax;
    break;
  };
  case 2: {
    mono = ay;
    break;
  };
  case 3: {
    mono = ax * ax;
    break;
  };
  case 4: {
    mono = ay * ay;
    break;
  };
  case 5: {
    mono = ax * ay;
    break;
  };
  };
  return (mono);
}