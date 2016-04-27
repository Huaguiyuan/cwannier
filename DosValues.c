#include "DosValues.h"

double* DosValues(HTightBinding *Hrs, gsl_matrix *R, int na, int nb, int nc, double *Es, double num_dos, bool all_Es) {
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
    if (all_Es) {
        double *dos_vals = Tetra_AllDosList(Efn, na, nb, nc, num_bands, R, Es, num_dos);
        return dos_vals;
    } else {
        double *dos_vals = Tetra_DosList(Efn, na, nb, nc, num_bands, R, Es, num_dos);
        return dos_vals;
    }
}

double* DosEnergyDerivValues(HTightBinding *Hrs, gsl_matrix *R, int na, int nb, int nc, double *Es, int num_dos, double num_electrons, double *fermi, double *dos_fermi, double *dos_deriv_fermi) {
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
    double *dos_deriv_vals = Tetra_DosEnergyDerivList(Efn, na, nb, nc, num_bands, R, Es, num_dos, num_electrons, fermi, dos_fermi, dos_deriv_fermi);
    return dos_deriv_vals;
}
