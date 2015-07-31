#include "DosValues.h"

double* DosValues(HTightBinding *Hrs, gsl_matrix *R, int num_k_per_dim, double *Es, double num_dos) {
    int num_bands = Hrs->num_bands;
    // Use GCC nested function to make closure.
    // https://gcc.gnu.org/onlinedocs/gcc/Nested-Functions.html
    // Efn puts the eigenvalues of H(k) into energies.
    void Efn(double k[3], gsl_vector *energies) {
        // For large Wannier system, no performance impact of allocating
        // Hk and work here versus outside Efn.
        // Allocate here to make Efn thread-safe (as long as each thread
        // calls with different `energies`).
        gsl_matrix_complex *Hk = gsl_matrix_complex_calloc(num_bands, num_bands);
        gsl_eigen_herm_workspace *work = gsl_eigen_herm_alloc(num_bands);
        // Set Hk = H(k).
        HkRecip(Hrs, k, Hk);
        // Calculate eigenvalues.
        gsl_eigen_herm(Hk, energies, work);
        gsl_sort_vector(energies);
        // Clean up.
        gsl_eigen_herm_free(work);
        gsl_matrix_complex_free(Hk);
    }
    // Calculate DOS values.
    double *dos_vals = Tetra_DosList(Efn, num_k_per_dim, num_bands, R, Es, num_dos);
    return dos_vals;
}
