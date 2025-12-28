# libGCK.c Geostatistics Library Documentation

The `libGCK.c` library implements a suite of functions for geostatistical estimation, primarily focusing on **Universal Kriging (UK)** using **Generalized Covariance (GC)** functions and the **Generalized Increment of order k (GIK)** method. The library is typically compiled as a C shared library for use in higher-level languages like Python.

## I. Utility and Neighbor Search Functions

These functions handle spatial organization, neighbor selection, and array manipulation essential for local Kriging estimates.

### `int octan(float ax, float ay)`

Determines the octant (0 to 7) for a given 2D vector defined by coordinates `(ax, ay)`. This is used to ensure spatial dispersion of selected neighbors.

| Parameter | Description |
| :--- | :--- |
| `ax` | X-component of the vector (e.g., $x_{estimate} - x_{neighbor}$). |
| `ay` | Y-component of the vector (e.g., $y_{estimate} - y_{neighbor}$). |
| **Returns** | The octant number (0-7). |

### `void quita(int ntot, int *neig, float *dist, int *oct)`

Shifts array elements one position to the left, effectively removing the first element. Used to discard the point being estimated (in cross-validation) from the list of neighbors.

| Parameter | Description |
| :--- | :--- |
| `ntot` | Total number of elements in the arrays. |
| `neig` | Array of neighbor indices. |
| `dist` | Array of distances. |
| `oct` | Array of octant numbers. |

### `int octavino(int nvec, int ntot, int *oct, float *dist, int *neig)`

Selects the best set of neighbors by enforcing an octant-based distribution (max 6 neighbors per octant) to ensure spatial coverage around the estimation point.

| Parameter | Description |
| :--- | :--- |
| `nvec` | The desired number of neighbors for Kriging. |
| `ntot` | The total number of neighbors available. |
| `oct` | Array of octant numbers (input/output). |
| `dist` | Array of distances (input/output). |
| `neig` | Array of neighbor indices (input/output). |
| **Returns** | `0` for successful neighbor selection, `1` if minimum spatial coverage is not met. |

### `int busca(int nvec, int ntot, int *oct, float *dist, int *neig)`

Main neighbor search routine used for **Cross-Validation (`main_crv2`)** and **GIK Calculation (`main_gik`)**. It sorts all points by distance, calls `quita` (to remove the point itself), and then calls `octavino`.

| Parameter | Description |
| :--- | :--- |
| `nvec` | The desired number of neighbors for Kriging. |
| `ntot` | The total number of points in the dataset. |
| `oct`, `dist`, `neig` | Arrays containing octant, distance, and index data. |
| **Returns** | The result of `octavino` (`0` for success, `1` for failure). |

### `int busca2(int nvec, int ntot, int *oct, float *dist, int *neig)`

Neighbor search routine used for **Interpolation (`main_faik`)**. It sorts points by distance and calls `octavino`, but **omits the call to `quita`**, meaning the nearest point (potentially the estimation point itself if it coincides with a data point) is kept.

| Parameter | Description |
| :--- | :--- |
| `nvec` | The desired number of neighbors for Kriging. |
| `ntot` | The total number of points in the dataset. |
| `oct`, `dist`, `neig` | Arrays containing octant, distance, and index data. |
| **Returns** | The result of `octavino` (`0` for success, `1` for failure). |

## II. Mathematical Kernels and Models

These functions define the core algebraic and geostatistical models used in the Kriging system.

### `int gaussj(float *a, float *y, int ne, int ns)`

Solves a system of linear equations $\mathbf{A}\mathbf{x} = \mathbf{Y}$ using the **Gauss-Jordan elimination** method with partial pivoting.

| Parameter | Description |
| :--- | :--- |
| `a` | Pointer to the flattened input matrix **A** ($\text{size } \text{ne} \times \text{ne}$). Overwritten during calculation. |
| `y` | Pointer to the right-hand side vector(s) **Y** ($\text{size } \text{ne} \times \text{ns}$). Overwritten by the solution **x**. |
| `ne` | The number of equations (dimension of matrix **A**). |
| `ns` | The number of right-hand side vectors. |
| **Returns** | `1` for successful solution, `0` if the matrix is singular. |

### `float monomio(float ax, float ay, int i)`

Calculates the value of a 2D polynomial monomial used for the **Drift (Trend)** in Universal Kriging.

| Index (`i`) | Monomial $m(x,y)$ | Kriging Order |
| :---: | :--- | :--- |
| 0 | 1 | 0 (Simple), 1 (Ordinary), 2 (Universal) |
| 1 | $x$ | 1, 2 |
| 2 | $y$ | 1, 2 |
| 3 | $x^2$ | 2 |
| 4 | $y^2$ | 2 |
| 5 | $x y$ | 2 |

| Parameter | Description |
| :--- | :--- |
| `ax` | X coordinate. |
| `ay` | Y coordinate. |
| `i` | Index of the monomial (0 to 5). |
| **Returns** | The value of the specified monomial at `(ax, ay)`. |

### `float estruc_gki(float h, float *zk)`

Calculates the Generalized Structure Function (or Variogram) $\gamma(h)$ using the fitted Generalized Covariance (GC) model parameters. Returns the negative of the calculated value, typically used as $C(h) = - \gamma(h)$.

The model form is: $\gamma(h) = C_0 + C_1|h| + C_3|h|^3 + C_5|h|^5 + C_{\log}|h|^2 \log(|h|)$.

| Parameter | Description |
| :--- | :--- |
| `h` | Lag distance. |
| `zk` | Array of 5 structure/covariance parameters ($C_0, C_1, C_3, C_5, C_{\log}$). |
| **Returns** | $-\gamma(h)$. |

### `float estruc_gki_1(float h)`

A simplified structure function model: $\gamma(h) = h$. Used specifically for calculating the Kriging system related to the Generalized Increments of order k (GIK).

| Parameter | Description |
| :--- | :--- |
| `h` | Lag distance. |
| **Returns** | The lag distance `h`. |

### `float covmodel(float h, int i)`

Calculates specific terms of the Generalized Covariance model based on index $i$. Used in fitting the GC model via the EGC (Experimental Generalized Covariance) calculation.

| Index (`i`) | Model Term $f(h)$ |
| :---: | :--- |
| 0 | $|h|$ |
| 1 | $|h|^3$ |
| 2 | $|h|^5$ |
| 3 | $|h|^2 \log(|h|)$ |

| Parameter | Description |
| :--- | :--- |
| `h` | Lag distance. |
| `i` | Index of the model term (0 to 3). |
| **Returns** | The value of the specified model term. |

## III. Kriging Solvers

These functions build and solve the Kriging system for weights ($\lambda_i$) and the Lagrange multiplier ($\mu$).

### `void krige(...)`

Performs **Universal Kriging (UK)** estimation for a single point. It constructs the Kriging matrix using `estruc_gki` (the full GC model) and solves it using `gaussj`.

| Parameter | Description |
| :--- | :--- |
| `mat` | Kriging system matrix $\mathbf{A}$ (input/output). |
| `vec` | Kriging right-hand side vector $\mathbf{b}$ (input/output, holds weights $\lambda_i$ after solve). |
| `vec1` | Copy of vector $\mathbf{b}$ for variance calculation. |
| `nvec` | Number of neighbors used. |
| `nork` | Kriging order (0: Simple, 1: Ordinary, 2: Universal). |
| `dist`, `neig` | Distances and indices of selected neighbors. |
| `x`, `y`, `z` | Arrays of coordinates and values for all data points. |
| `zestim` | Output: Estimated Z value. |
| `sigma` | Output: Kriging standard deviation ($\sqrt{\sigma^2}$). |
| `xx`, `yy` | Coordinates of the point to estimate. |
| `zk` | Array of 5 GC model parameters. |

### `void krigea_gik_1(...)`

Performs Kriging estimation specifically for the **GIK (Generalized Increment of order k)**. It uses the simplified linear structure function $\gamma(h) = h$ (`estruc_gki_1`) to calculate the weights necessary to define the GIK.

*(Parameters are identical to `krige`, except for the missing `zk` parameter as the structure function is fixed.)*

## IV. Workflow Main Functions

These functions orchestrate the main geostatistical tasks: GIK calculation, model fitting, cross-validation, and interpolation.

### `void main_gik(int ntot, float *x, float *y, float *z, int nork, int nvec)`

The primary function for the **Generalized Covariance approach**. It iterates through all data points, calculates the GIK weights using `krigea_gik_1`, writes the GIK definitions (`gki.dat`) and the contributions to the Experimental Generalized Covariance (`cov.dat`), and finally calls the model fitting and cross-validation routines.

**Workflow Steps:**
1. Calculate GIKs and EGC contributions.
2. Call `main_geko()`.
3. Call `main_crv2()`.

### `void main_geko()`

**GEneralized KOvariance Model Fitting.** Reads the EGC terms and GIK contributions from `cov.dat`, sets up and solves a **weighted least squares system** to fit the best GC model parameters (`$C_0, C_1, C_3, C_5, C_{\log}$`) based on model configurations defined in `modelos.dat`. Writes the resulting stable GC models to `modcov.dat`.

### `void main_crv2(int ntot, ...)`

**CRoss-Validation (Leave-One-Out).** Reads the fitted GC models from `modcov.dat`, and performs Leave-One-Out Cross-Validation for each model across all data points. Calculates the Mean Squared Error (MSE) and other statistics. Writes the CV results to `cv.dat`.

### `void main_faik(int ntot, ..., int bins, int hist)`

**FAIK Kriging Interpolation.** Performs Universal Kriging interpolation over a defined grid (from `xmin, ymin` to `xmax, ymax`) using a specified GC model (`zk`).

| Parameter | Description |
| :--- | :--- |
| `zk[5]` | Array of 5 GC model parameters to use for interpolation. |
| `xmin`, `ymin` | Minimum coordinates of the interpolation grid. |
| `xmax`, `ymax` | Maximum coordinates of the interpolation grid. |
| `bins` | Number of cells in the X direction (columns). |
| `hist` | Number of cells in the Y direction (rows). |
| **Output** | Writes the estimates to `faik.grd` and the Kriging standard deviations to `faike.grd` ([DSAA grid format](https://surferhelp.goldensoftware.com/topics/ascii_grid_file_format.htm)). |

## V. External Files

The library relies on several files for configuration, input, and output:

| File Name | R/W | Description |
| :--- | :--- | :--- |
| **`modelos.dat`** | Read | Configuration file defining which GC model terms ($C_0, C_1, C_3, C_5, C_{\log}$) to fit. Used by `main_geko`. |
| **`gki.dat`** | Write | Stores the definition and weights of each calculated Generalized Increment of order k (GIK). Used by `escribe_gki`. |
| **`cov.dat`** | Write/Read | Stores the contributions of each GIK to the **Experimental Generalized Covariance (EGC)**. Written by `escribe_gki`, read by `main_geko`. |
| **`modcov.dat`** | Write/Read | Contains the parameters of the stable and fitted Generalized Covariance models. Written by `main_geko`, read by `main_crv2`. |
| **`cv.dat`** | Write | Results of the Cross-Validation, including the Mean Squared Error (MSE) for each fitted model. Written by `main_crv2`. |
| **`faik.grd`** | Write | Grid file ([DSAA format](https://surferhelp.goldensoftware.com/topics/ascii_grid_file_format.htm)) containing the Kriging estimated values. Written by `main_faik`. |
| **`faike.grd`** | Write | Grid file ([DSAA format](https://surferhelp.goldensoftware.com/topics/ascii_grid_file_format.htm)) containing the Kriging standard deviation (error) values. Written by `main_faik`. |