# libGCK: Universal Kriging with Generalized Covariance Library (C)

`libGCK` is a lightweight C library implementing geostatistical routines, primarily focusing on **Universal Kriging (UK)** using the **Generalized Covariance (GC)** function approach. It utilizes the Generalized Increment of order k (GIK) method for robust covariance modeling and includes functionality for model fitting, cross-validation, and grid interpolation.

This library is designed to be compiled as a shared object (`.so` or `.dll`) for integration with other systems, such as Python via `ctypes`.

## 🌟 Features

* **Universal Kriging (UK):** Supports Kriging with a trend (up to quadratic drift).
* **Generalized Covariance (GC):** Implements the GIK method to calculate the Experimental Generalized Covariance (EGC) and fit stable GC models ($C_0, C_1|h|, C_3|h|^3, C_5|h|^5, C_{\log}|h|^2 \log|h|$).
* **Octant-Based Neighbor Search:** Includes specialized functions (`octan`, `busca`, `octavino`) to select neighbors that are spatially dispersed, improving the stability of the Kriging system.
* **Gauss-Jordan Solver:** Includes a robust `gaussj` function for solving the Kriging system of equations.
* **Cross-Validation:** Automated Leave-One-Out Cross-Validation (LOOCV) for fitted models.
* **Gridded Interpolation:** Output interpolated estimates and Kriging error in the common [DSAA grid format](https://surferhelp.goldensoftware.com/topics/ascii_grid_file_format.htm).

## 🛠️ Compilation

The library requires standard C compilation and linking against the math library (`-lm`).

```bash
# Compile as a shared library (e.g., for Python integration)
gcc -shared -o libGCK.so libGCK.c -lm
