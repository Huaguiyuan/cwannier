#include "HTightBinding.h"
#include "BandEnergy.h"

double BandEnergy(HTightBinding *Hrs, double R[3][3], double num_electrons, int n0, double tol) {
    int num_bands = Hrs->num_bands;
    gsl_matrix_complex *Hk = gsl_matrix_complex_alloc(num_bands, num_bands);
    gsl_eigen_herm_workspace *work = gsl_eigen_herm_alloc(num_bands);
    // Use GCC nested function to make closure.
    // https://gcc.gnu.org/onlinedocs/gcc/Nested-Functions.html
    // Efn puts the eigenvalues of H(k) into energies.
    void Efn(double k[3], gsl_vector *energies) {
        // Zero out Hk.
        gsl_matrix_complex_set_zero(Hk);
        // Set Hk = H(k).
        HkRecip(Hrs, k, Hk);
        // Calculate eigenvalues.
        gsl_eigen_herm(Hk, energies, work);
    }
    // Calculate band energy.
    double esum = 0.0;
    // Clean up.
    gsl_eigen_herm_free(work);
    gsl_matrix_complex_free(Hk);

    return esum;
}
