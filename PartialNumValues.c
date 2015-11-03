#include "PartialNumValues.h"

double** PartialNumValues(HTightBinding *Hrs, gsl_matrix *R, int num_k_per_dim, double num_total_electrons, double **Es, double num_E, double *E_Fermi, double **num_states_Fermi) {
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
    // Calculate num values and E_Fermi.
    double **num_vals = partial_num_states(UEfn, num_k_per_dim, num_bands, num_total_electrons, R, Es, num_E, E_Fermi, num_states_Fermi);
    return num_vals;
}
