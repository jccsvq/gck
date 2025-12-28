#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Determines the octant (0 to 7) for a given 2D vector (ax, ay).
 *
 * This function is crucial for organizing nearest neighbors in an octant-based
 * search strategy (like octant search in Kriging) to ensure even spatial coverage
 * around the point being estimated.
 * The octants are typically defined relative to the axes and the lines y=x and y=-x.
 *
 * @param ax X-component of the vector (x_estimate - x_neighbor).
 * @param ay Y-component of the vector (y_estimate - y_neighbor).
 * @return int The octant number (0-7).
 */
int octan(float ax, float ay) {
  int oc;
  if (ax < 0.) {
    // Left half-plane (Octants 2, 3, 4, 5)
    if (ay < 0.) {
      // Third quadrant (Octants 4, 5)
      if ((-ax) >= (-ay)) {
        oc = 4; // Between -X axis and y=-x line
      } else {
        oc = 5; // Between y=-x line and -Y axis
      };
    } else {
      // Second quadrant (Octants 2, 3)
      if ((-ax) >= ay) {
        oc = 3; // Between -X axis and y=x line
      } else {
        oc = 2; // Between y=x line and +Y axis
      };
    };
  } else {
    // Right half-plane (Octants 0, 1, 6, 7)
    if (ay < 0.) {
      // Fourth quadrant (Octants 6, 7)
      if (ax >= (-ay)) {
        oc = 7; // Between +X axis and y=-x line
      } else {
        oc = 6; // Between y=-x line and -Y axis
      };
    } else {
      // First quadrant (Octants 0, 1)
      if (ax >= ay) {
        oc = 0; // Between +X axis and y=x line
      } else {
        oc = 1; // Between y=x line and +Y axis
      };
    };
  };
  return (oc);
}

/**
 * @brief Shifts array elements to the left, effectively removing the first element.
 *
 * This is used to remove the first nearest neighbor (often the point itself in cross-validation)
 * or to simply shift a list of neighbors, distances, and octants.
 *
 * @param ntot Total number of elements in the arrays before shifting.
 * @param neig Array of neighbor indices.
 * @param dist Array of distances.
 * @param oct Array of octant numbers.
 */
void quita(int ntot, int *neig, float *dist, int *oct) {
  int i;
  for (i = 0; i < (ntot - 1); i++) {
    neig[i] = neig[i + 1];
    dist[i] = dist[i + 1];
    oct[i] = oct[i + 1];
  };
}

/**
 * @brief Selects the best neighbors by octant to ensure spatial dispersion.
 *
 * This function groups neighbors by their octant and selects up to a maximum
 * (currently hardcoded as 6) per octant. It checks if the minimum number of
 * non-empty octants is met (e.g., nonul <= 5 means insufficient spread).
 * Finally, it rebuilds the neighbor arrays with the selected neighbors.
 *
 * @param nvec The desired number of neighbors for Kriging.
 * @param ntot The total number of neighbors available.
 * @param oct Array of octant numbers for each neighbor.
 * @param dist Array of distances for each neighbor.
 * @param neig Array of neighbor indices.
 * @return int 1 if the neighbor selection process failed (e.g., too few non-null octants), 0 otherwise.
 */
int octavino(int nvec, int ntot, int *oct, float *dist, int *neig) {
  int n[8], na[7][8], i, j, nnn, nonul, ifrac;
  float b[7][8];
  // 1. Initialize octant counters and temporary neighbor/distance storage
  for (i = 0; i < 8; i++)
    n[i] = 0; // n[i] stores the count of neighbors in octant i
  for (i = 0; i < 7; i++)
    for (j = 0; j < 8; j++) {
      na[i][j] = -1;
      b[i][j] = 0.;
    };
  // 2. Distribute neighbors into octant bins (up to 6 per octant)
  for (i = 0; i < ntot; i++) {
    if (n[oct[i]] < 6) { // Hardcoded limit of 6 neighbors per octant
      n[oct[i]] += 1;
      na[n[oct[i]] - 1][oct[i]] = neig[i];
      b[n[oct[i]] - 1][oct[i]] = dist[i];
    };
  };
  // 3. Check for minimum spatial coverage (at least 6 non-empty octants)
  nonul = 0;
  for (i = 0; i < 8; i++)
    if (n[i] != 0)
      nonul++;
  if (nonul <= 5) // If 5 or less octants have neighbors, return failure (1)
    return (1);
  // 4. Rebuild the neighbor arrays (neig and dist) from the temporary octant storage
  nnn = -1;
  for (i = 0; i < 7; i++)
    for (j = 0; j < 8; j++)
      if (na[i][j] != -1) {
        nnn++;
        neig[nnn] = na[i][j];
        dist[nnn] = b[i][j];
      };
  // 5. Check if the total number of selected neighbors is sufficient
  if (nnn >= nvec) {
    ifrac = 0; // Success: sufficient neighbors selected
  } else {
    ifrac = 1; // Failure: not enough neighbors selected (although 5 octants check is primary)
  };
  return (ifrac);
}

/**
 * @brief Sorts neighbors by distance and performs neighbor selection (including octant check).
 *
 * Sorts the initial list of neighbors by distance (nearest first) and then
 * calls `quita` to remove the first element (e.g., the point itself in cross-validation)
 * before performing the octant-based neighbor selection via `octavino`.
 *
 * @param nvec The desired number of neighbors for Kriging.
 * @param ntot The total number of points in the dataset.
 * @param oct Array of octant numbers (output/input).
 * @param dist Array of distances (output/input).
 * @param neig Array of neighbor indices (output/input).
 * @return int The result of `octavino` (0 for success, 1 for failure).
 */
int busca(int nvec, int ntot, int *oct, float *dist, int *neig) {
  int i, j, k, ifrac = 0;
  float w;
  // Initialize neighbor indices
  for (i = 0; i < ntot; i++)
    neig[i] = i;
  // 1. Sort neighbors by distance using Bubble Sort (slow for large N)
  for (i = 0; i < (ntot - 1); i++)
    for (j = (i + 1); j < ntot; j++) {
      if (dist[i] > dist[j]) {
        // Swap distance, neighbor index, and octant simultaneously
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
  // 2. Remove the first element (closest neighbor, often the point itself)
  quita(ntot, neig, dist, oct);
  // 3. Perform octant-based neighbor selection
  ifrac = octavino(nvec, ntot, oct, dist, neig);

  return (ifrac);
}

/**
 * @brief Sorts neighbors by distance and performs octant-based selection (without quita).
 *
 * This function is similar to `busca` but **does not call `quita`**. This means
 * the closest point (often the point itself being estimated in the full dataset, *not* cross-validation)
 * is **not** removed from the sorted list before octant selection. Used in interpolation (`main_faik`).
 *
 * @param nvec The desired number of neighbors for Kriging.
 * @param ntot The total number of points in the dataset.
 * @param oct Array of octant numbers (output/input).
 * @param dist Array of distances (output/input).
 * @param neig Array of neighbor indices (output/input).
 * @return int The result of `octavino` (0 for success, 1 for failure).
 */
int busca2(int nvec, int ntot, int *oct, float *dist, int *neig) {
  int i, j, k, ifrac = 0;
  float w;
  // Initialize neighbor indices
  for (i = 0; i < ntot; i++)
    neig[i] = i;
  // 1. Sort neighbors by distance
  for (i = 0; i < (ntot - 1); i++)
    for (j = (i + 1); j < ntot; j++) {
      if (dist[i] > dist[j]) {
        // Swap distance, neighbor index, and octant
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
  /*quita(ntot, neig, dist, oct); // <-- Commented out/omitted: Keeps the closest point*/
  // 2. Perform octant-based neighbor selection
  ifrac = octavino(nvec, ntot, oct, dist, neig);

  return (ifrac);
}

/**
 * @brief Solves a system of linear equations Ax = Y using the Gauss-Jordan method.
 *
 * The matrix A and the vector Y are modified in place. The solution (x) overwrites Y.
 * It uses partial pivoting (searching for the maximum pivot in the current column)
 * to improve numerical stability. Handles multiple right-hand side vectors (ns > 1).
 *
 * @param a Pointer to the flattened square matrix A (size ne * ne).
 * @param y Pointer to the right-hand side vector(s) Y (size ne * ns).
 * @param ne The number of equations (dimension of matrix A).
 * @param ns The number of right-hand side vectors (number of columns in Y).
 * @return int 1 for successful solution, 0 if the matrix is singular (pivot is zero).
 */
int gaussj(float *a, float *y, int ne, int ns)
{
  float *(*ap), val, *tp, *(*yp);
  int i, j, k, piv;
  // 1. Allocate memory for array of pointers to access rows easily (ap for A, yp for Y)
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

  // 2. Main Gauss-Jordan loop
  for (i = 0; i < ne; i++)
  /* for each row */
  {
    /* search for maximum pivot */
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
      return (0); /* Singular matrix */
    };

    /* swap rows (pivoting) */
    tp = ap[i];
    ap[i] = ap[piv];
    ap[piv] = tp;
    tp = yp[i];
    yp[i] = yp[piv];
    yp[piv] = tp;
    val = ap[i][i];
    /* normalize equation i (divide by pivot) */
    for (j = i; j < ne; j++)
      ap[i][j] /= val;
    for (j = 0; j < ns; j++)
      yp[i][j] /= val;
    /* reduce other equations (elimination) */
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

  // 3. Copy solution from yp back to the original Y array
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

  // 4. Free allocated memory
  free(ap);
  free(yp);
  return (1);
}

/**
 * @brief Generalized Structure Function (or Variogram) model for Kriging.
 *
 * This model seems to correspond to a Generalized Covariance (GC) of order k.
 * The parameters are stored in the array zk:
 * zk[0] is not used in the h>0 branch.
 * zk[1] * h (linear term)
 * zk[2] * h^3 (cubic term)
 * zk[3] * h^5 (fifth power term)
 * zk[4] * h^2 * log(h) (logarithmic term)
 * The function returns the negative of the calculated value.
 *
 * @param h Lag distance.
 * @param zk Array of 5 structure/covariance parameters.
 * @return float Negative of the structure function value.
 */
float estruc_gki(float h, float *zk) 
{
  float est;
  if (h > 0.) {
    // est = C1*h + C3*h^3 + C5*h^5 + Clog*h^2*log(h)
    est = h * zk[1] + h * h * h * zk[2] + h * h * h * h * h * zk[3] +
          h * h * zk[4] * (float)log((double)h);
  } else {
    // Nugget effect (C0 or -C0) when h=0
    est = zk[0];
  };
  return ((-est)); // Returns the negative value
}

/**
 * @brief Simple linear structure function model (gamma(h) = h).
 *
 * Used primarily for calculating Generalized Increments of order k (GIK) weights.
 *
 * @param h Lag distance.
 * @return float The lag distance h itself.
 */
float estruc_gki_1(float h)
{

  return ((float)h);
}

/**
 * @brief Calculates a 2D polynomial monomial value based on index i.
 *
 * These monomials form the drift (trend) for Universal Kriging (Kriging with a trend).
 * - i=0: Constant (k=0, m(x,y)=1)
 * - i=1: Linear X (k=1, m(x,y)=x)
 * - i=2: Linear Y (k=1, m(x,y)=y)
 * - i=3: Quadratic X^2 (k=2, m(x,y)=x^2)
 * - i=4: Quadratic Y^2 (k=2, m(x,y)=y^2)
 * - i=5: Quadratic XY (k=2, m(x,y)=x*y)
 *
 * @param ax X coordinate.
 * @param ay Y coordinate.
 * @param i Index of the monomial to calculate (0 to 5).
 * @return float The value of the specified monomial at (ax, ay).
 */
float monomio(float ax, float ay, int i) {
  float mono;
  mono = 0.0;
  switch (i) {
  case 0: {
    mono = 1.; // Constant
    break;
  };
  case 1: {
    mono = ax; // Linear X
    break;
  };
  case 2: {
    mono = ay; // Linear Y
    break;
  };
  case 3: {
    mono = ax * ax; // Quadratic X^2
    break;
  };
  case 4: {
    mono = ay * ay; // Quadratic Y^2
    break;
  };
  case 5: {
    mono = ax * ay; // Quadratic XY
    break;
  };
  };
  return (mono);
}

/**
 * @brief Performs Universal Kriging (UK) estimation for a single point.
 *
 * Constructs and solves the Kriging system of equations for the weights (lambda)
 * and the Lagrange multiplier (mu). It uses the Generalized Structure Function
 * from `estruc_gki` and the `gaussj` solver.
 *
 * @param mat Matrix A for the Kriging system (input/output).
 * @param vec Vector b for the Kriging system (input/output, holds weights after solve).
 * @param vec1 Copy of vector b for variance calculation (input/output).
 * @param nvec Number of neighbors used.
 * @param nork Order of Kriging (0: Simple, 1: Ordinary, 2: Universal with Quadratic Drift).
 * @param dist Array of distances to selected neighbors.
 * @param neig Indices of selected neighbors.
 * @param x Array of X coordinates of all points.
 * @param y Array of Y coordinates of all points.
 * @param z Array of Z values of all points.
 * @param zestim Output: Estimated Z value.
 * @param sigma Output: Kriging variance (square root of variance).
 * @param xx X coordinate of the point to estimate.
 * @param yy Y coordinate of the point to estimate.
 * @param zk Array of 5 structure/covariance parameters.
 */
void krige(float *mat, float *vec, float *vec1, int nvec, int nork, float *dist,
           int *neig, float *x, float *y, float *z, float *zestim, float *sigma,
           float xx, float yy, float *zk)

{
  float *(*a), *b, *bb, c0;
  int i, ii, j, jj, dim, nmon[] = {1, 3, 6}; // Monomial counts for Kriging order 0, 1, 2
  // 1. Determine the dimension of the Kriging system
  dim = nvec + nmon[nork];
  // 2. Allocate pointers for matrix A and vectors b and bb
  if ((a = calloc(dim, sizeof(float *))) == NULL) {
    printf("\nERROR: Insufficient memory for equations!");
    abort();
  };
  for (i = 0; i < dim; i++)
    a[i] = mat + i * dim;
  b = vec;
  bb = vec1;

  /* Establish the Kriging equations */

  /* Right-hand side vector (b) */
  for (i = 0; i < nvec; i++) {
    b[i] = estruc_gki(dist[i], zk); // Structure function for distance to neighbor i
    bb[i] = b[i];                  // Copy for variance calculation
  };
  // Monomial terms (constraint equations)
  for (i = 0; i < nmon[nork]; i++) {
    b[nvec + i] = -monomio(xx, yy, i); // Monomials at the estimation point (xx, yy)
    bb[nvec + i] = b[nvec + i];
  };

  /* Principal Diagonal */
  c0 = estruc_gki(0., zk); // The nugget effect (C(0) or -gamma(0))
  for (i = 0; i < nvec; i++)
    a[i][i] = c0;

  /* Null box (Lagrange multiplier block) */
  for (i = 0; i < nmon[nork]; i++)
    for (j = 0; j < nmon[nork]; j++)
      a[nvec + i][nvec + j] = 0.; // Last block of the matrix is zero

  /* Symmetric part of the box (Covariance/Structure between neighbors) */
  for (i = 0; i < (nvec - 1); i++) {
    ii = neig[i];
    for (j = i + 1; j < nvec; j++) {
      jj = neig[j];
      // Calculate distance between neighbor ii and jj
      float h_sq = (x[ii] - x[jj]) * (x[ii] - x[jj]) + (y[ii] - y[jj]) * (y[ii] - y[jj]);
      a[i][j] = estruc_gki((float)sqrt(h_sq), zk);
      a[j][i] = a[i][j];
    };
  };

  /* Lateral boxes (Monomials at neighbor locations) */
  for (i = 0; i < nvec; i++) {
    ii = neig[i];
    for (j = 0; j < nmon[nork]; j++) {
      a[i][nvec + j] = -monomio(x[ii], y[ii], j); // Monomials at neighbor ii
      a[nvec + j][i] = a[i][nvec + j];            // Symmetry
    };
  };

  /* Solve the equations */
  if ((gaussj(mat, vec, dim, 1)) == 1) { // Solve the system (weights are now in b)
    *zestim = *sigma = 0.;
    // Kriging Estimate: Z* = sum(lambda_i * Z_i)
    for (i = 0; i < nvec; i++)
      *zestim += b[i] * z[neig[i]];
    // Kriging Variance: sigma^2 = sum(lambda_i * b_i) - c0
    for (i = 0; i < dim; i++)
      *sigma += b[i] * bb[i];
    if ((*sigma -= c0) > 0.) {
      *sigma = (float)sqrt(1. * (*sigma)); // Square root for standard deviation
    } else {
      *sigma = 9e-29;
      *zestim = -*sigma;
    };
  } else {
    // Failure to solve
    *sigma = 9e-29;
    *zestim = -*sigma;
  };

  free(a);
}

/**
 * @brief Performs Kriging to calculate weights for a Generalized Increment of order k (GIK).
 *
 * This function is similar to `krige` but uses the simple linear structure function
 * `estruc_gki_1` (gamma(h) = h) instead of the general one. The resulting weights
 * are used to calculate the GIK (which is essentially a generalized residual).
 *
 * @param mat Matrix A for the Kriging system.
 * @param vec Vector b for the Kriging system (holds weights after solve).
 * @param vec1 Copy of vector b for variance calculation.
 * @param nvec Number of neighbors used.
 * @param nork Order of Kriging.
 * @param dist Array of distances to selected neighbors.
 * @param neig Indices of selected neighbors.
 * @param x Array of X coordinates.
 * @param y Array of Y coordinates.
 * @param z Array of Z values.
 * @param zestim Output: Estimated Z value (using these weights).
 * @param sigma Output: Kriging variance.
 * @param xx X coordinate of the point to estimate.
 * @param yy Y coordinate of the point to estimate.
 */
void krigea_gik_1(float *mat, float *vec, float *vec1, int nvec, int nork,
                  float *dist, int *neig, float *x, float *y, float *z,
                  float *zestim, float *sigma, float xx, float yy) {
  // Structure is identical to `krige`, but uses `estruc_gki_1` (gamma(h)=h)
  float *(*a), *b, *bb, c0;
  int i, ii, j, jj, dim, nmon[] = {1, 3, 6};
  // ... (Memory allocation and system setup, identical to krige)
  dim = nvec + nmon[nork];
  if ((a = calloc(dim, sizeof(float *))) == NULL) {
    printf("\nERROR: Insufficient memory for equations!");
    abort();
  };
  for (i = 0; i < dim; i++)
    a[i] = mat + i * dim;
  b = vec;
  bb = vec1;

  /* Establish the Kriging equations (using estruc_gki_1) */

  /* Right-hand side vector (b) */
  for (i = 0; i < nvec; i++) {
    b[i] = estruc_gki_1(dist[i]); // Uses gamma(h) = h
    bb[i] = b[i];
  };
  for (i = 0; i < nmon[nork]; i++) {
    b[nvec + i] = -monomio(xx, yy, i);
    bb[nvec + i] = b[nvec + i];
  };

  /* Principal Diagonal */
  c0 = estruc_gki_1(0.); // gamma(0) = 0
  for (i = 0; i < nvec; i++)
    a[i][i] = c0;

  /* Null box */
  for (i = 0; i < nmon[nork]; i++)
    for (j = 0; j < nmon[nork]; j++)
      a[nvec + i][nvec + j] = 0.;

  /* Symmetric part of the box (Structure between neighbors) */
  for (i = 0; i < (nvec - 1); i++) {
    ii = neig[i];
    for (j = i + 1; j < nvec; j++) {
      jj = neig[j];
      float h_sq = (x[ii] - x[jj]) * (x[ii] - x[jj]) + (y[ii] - y[jj]) * (y[ii] - y[jj]);
      a[i][j] = estruc_gki_1((float)sqrt(h_sq)); // Uses gamma(h) = h
      a[j][i] = a[i][j];
    };
  };

  /* Lateral boxes */
  for (i = 0; i < nvec; i++) {
    ii = neig[i];
    for (j = 0; j < nmon[nork]; j++) {
      a[i][nvec + j] = -monomio(x[ii], y[ii], j);
      a[nvec + j][i] = a[i][nvec + j];
    };
  };

  /* Solve the equations */

  if ((gaussj(mat, vec, dim, 1)) == 1) {
    *zestim = *sigma = 0.;
    // Kriging Estimate (using the weights lambda and Z values)
    for (i = 0; i < nvec; i++)
      *zestim += b[i] * z[neig[i]];
    // Kriging Variance
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

/**
 * @brief Calculates a specific term of the generalized covariance model based on index i.
 *
 * This function returns the value of the term for a given lag distance h, used
 * in the calculation of the **Experimental Generalized Covariance** (EGC) for a GIK.
 * - i=0: h (linear)
 * - i=1: h^3 (cubic)
 * - i=2: h^5 (fifth power)
 * - i=3: h^2 * log(h) (logarithmic, typically 0 for h <= 0)
 *
 * @param h Lag distance.
 * @param i Index of the model term (0 to 3).
 * @return float The value of the model term.
 */
float covmodel(float h, int i)
{
  float mono;
  mono = 0.0;
  switch (i) {
  case 0: {
    mono = h;
    break;
  };
  case 1: {
    mono = h * h * h;
    break;
  };
  case 2: {
    mono = h * h * h * h * h;
    break;
  };
  case 3: {
    // h^2 * log(h) is 0 for h <= 0
    mono = (h <= 0.) ? 0. : (float)h * h * log((double)h);
    break;
  };
  };
  return (mono);
}

/**
 * @brief Calculates a specific term of the generalized covariance model (identical to covmodel).
 *
 * This function is a duplicate of `covmodel` and is used in the `escribe_gki` function.
 *
 * @param h Lag distance.
 * @param i Index of the model term (0 to 3).
 * @return float The value of the model term.
 */
float covmodel_1(float h, int i)
{
  float mono;
  mono = 0.0;
  switch (i) {
  case 0: {
    mono = h;
    break;
  };
  case 1: {
    mono = h * h * h;
    break;
  };
  case 2: {
    mono = h * h * h * h * h;
    break;
  };
  case 3: {
    mono = (h <= 0.) ? 0. : (float)h * h * log((double)h);
    break;
  };
  };
  return (mono);
}

/**
 * @brief Writes the Generalized Increment of order k (GIK) and its contribution to the Experimental Generalized Covariance (EGC) model.
 *
 * It takes the Kriging weights calculated by `krigea_gik_1` (stored in `ponde`) and the neighbor
 * setup to calculate and write the GIK value and the EGC terms to files.
 *
 * @param gki File pointer for GIK data output ("gki.dat").
 * @param cov File pointer for EGC contribution output ("cov.dat").
 * @param nvec Number of neighbors used.
 * @param nork Kriging order.
 * @param neig Array of neighbor indices.
 * @param x Array of X coordinates.
 * @param y Array of Y coordinates.
 * @param z Array of Z values.
 * @param ponde Kriging weights (input: lambda_i, output: normalized weights).
 * @param ix Index of the central point.
 */
void escribe_gki(FILE *gki, FILE *cov, int nvec, int nork, int *neig, float *x,
                 float *y, float *z, float *ponde, int ix) {
  int i, j, k;
  float sq, dx, dy, cx[4], cy, cw;
  // 1. Normalize the GIK weights (sum of squares is 1)
  // The weights in `ponde` (which is `vec` from `krigea_gik_1`) are modified:
  // ponde[0] is the Lagrange Multiplier (-mu), which is shifted to ponde[nvec] and negated
  // ponde[1] to ponde[nvec] are the Kriging weights (lambda_i), shifted and normalized
  sq = 1.; // Start with the central point (weight=1)
  for (i = nvec; i >= 1; i--) {
    ponde[i] = (-ponde[i - 1]); // The weights lambda_i are stored from index 1 to nvec
    sq += ponde[i] * ponde[i];
  };
  ponde[0] = 1.; // Weight of the central point (ix) is 1
  sq = (float)sqrt(1. * sq);
  for (i = 0; i <= nvec; i++)
    ponde[i] /= sq; // Normalize all weights
  // 2. Write the GIK configuration to gki.dat
  fprintf(gki, "%d %d\n%d ", nork, nvec + 1, ix);
  for (i = 0; i < nvec; i++)
    fprintf(gki, "%d ", neig[i]);
  fprintf(gki, "\n");
  for (i = 0; i <= nvec; i++)
    fprintf(gki, "%f ", ponde[i]);
  fprintf(gki, "\n");

  // 3. Calculate the GIK value (cy)
  cy = ponde[0] * z[ix]; // Weight of central point * Z_ix
  for (i = 1; i <= nvec; i++)
    cy += ponde[i] * z[neig[i - 1]]; // sum(lambda_i * Z_i)
  cy = cy * cy; // Squared GIK value

  // 4. Recalculate neighbor indices including the central point
  for (i = nvec; i > 0; i--)
    neig[i] = neig[i - 1];
  neig[0] = ix; // neig[0] is now the central point (ix)

  // 5. Calculate the EGC contributions (cx[k])
  for (k = 0; k < 4; k++)
    cx[k] = 0;
  for (i = 0; i <= nvec; i++)
    for (j = 0; j <= nvec; j++) {
      dx = x[neig[i]] - x[neig[j]];
      dy = y[neig[i]] - y[neig[j]];
      cw = (float)sqrt(1. * (dx * dx + dy * dy)); // Distance between points i and j
      for (k = 0; k < 4; k++)
        // EGC contribution: lambda_i * lambda_j * covmodel(h_ij, k)
        cx[k] += ponde[i] * ponde[j] * covmodel_1(cw, k);
    };
  // Inverse of the second EGC term (used as weight for the squared GIK value)
  cw = 1. / (cx[1] * cx[1]);
  // 6. Write the EGC terms and GIK contribution to cov.dat
  // Format: 1. cx[0] cx[1] cx[2] cx[3] cy cw (cw is the weight)
  fprintf(cov, "1. %e %e %e %e %e %e\n", cx[0], cx[1], cx[2], cx[3], cy, cw);
}

