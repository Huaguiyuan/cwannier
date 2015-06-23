#include "BandEnergy.h"

double BandEnergy(double *E_Fermi, HTightBinding *Hrs, gsl_matrix *R, double num_electrons, int n0, double tol) {
    int num_bands = Hrs->num_bands;
    // Use GCC nested function to make closure.
    // https://gcc.gnu.org/onlinedocs/gcc/Nested-Functions.html
    // Efn puts the eigenvalues of H(k) into energies.
    void Efn(double k[3], gsl_vector *energies) {
        // For large Wannier system, no performance impact of allocating
        // Hk and work here versus outside Efn.
        // Allocate here to make Efn thread-safe (as long as called with
        // different energies).
        gsl_matrix_complex *Hk = gsl_matrix_complex_alloc(num_bands, num_bands);
        gsl_eigen_herm_workspace *work = gsl_eigen_herm_alloc(num_bands);
        // Set Hk = H(k).
        HkRecip(Hrs, k, Hk);
        // Calculate eigenvalues.
        gsl_eigen_herm(Hk, energies, work);
        // Clean up.
        gsl_eigen_herm_free(work);
        gsl_matrix_complex_free(Hk);
    }
    // Calculate band energy.
    bool use_cache = true;
    double esum = SumEnergy(E_Fermi, Efn, n0, num_bands, num_electrons, R, use_cache);
    printf("esum = %f\n", esum);
    printf("E_Fermi = %f\n", *E_Fermi);

    return esum;
}
