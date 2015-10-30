#include "PartialDosValues.h"

double** PartialDosValues(HTightBinding *Hrs, gsl_matrix *R, int num_k_per_dim, double sigma, double **Es, double num_dos) {
    int num_bands = Hrs->num_bands;
    // Use GCC nested function to make closure.
    // https://gcc.gnu.org/onlinedocs/gcc/Nested-Functions.html
    // Efn puts the eigenvalues of H(k) into energies.
    void UEfn(double k[3], gsl_vector *energies, gsl_matrix_complex *evecs) {
        // For large Wannier system, no performance impact of allocating
        // Hk and work here versus outside Efn.
        // Allocate here to make Efn thread-safe (as long as each thread
        // calls with different `energies` and `evecs`).
        gsl_matrix_complex *Hk = gsl_matrix_complex_calloc(num_bands, num_bands);
        gsl_eigen_hermv_workspace *work = gsl_eigen_hermv_alloc(num_bands);
        // Set Hk = H(k).
        HkRecip(Hrs, k, Hk);
        // Calculate eigenvalues.
        gsl_eigen_hermv(Hk, energies, evecs, work);
        // Sort the eigenvalues and corresponding eigenvectors.
        sort_evals_evecs(energies, evecs, num_bands);
        // Clean up.
        gsl_eigen_hermv_free(work);
        gsl_matrix_complex_free(Hk);
    }
    // Calculate DOS values.
    double **dos_vals = Gauss_PartialDosList(UEfn, num_k_per_dim, sigma, num_bands, R, Es, num_dos);
    return dos_vals;
}

void sort_evals_evecs(gsl_vector *energies, gsl_matrix_complex *evecs, int num_bands) {
    int band_index;

    gsl_permutation *p = gsl_permutation_calloc(num_bands);

    gsl_sort_vector_index(p, energies);
    gsl_permute_vector(p, energies);

    for (band_index = 0; band_index < num_bands; band_index++) {
        // Sort each row: equivalent to permuting eigenvectors
        // (rows hold original basis amplitude for each eigenvector).
        gsl_vector_complex_view row = gsl_matrix_complex_row(evecs, band_index);
        gsl_permute_vector_complex(p, &(row.vector));
    }

    gsl_permutation_free(p);
}
