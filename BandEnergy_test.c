#include <gsl/gsl_matrix.h>
#include "BandEnergy.h"

int main(int argc, char *argv[]) {
    HTightBinding *Hrs = ExtractHTightBinding("test_data/Fe_soc/Fe_soc_hr.dat");

    gsl_matrix *R = gsl_matrix_calloc(3, 3); // R -> all zeros
    // Overall scale for R doesn't matter.
    gsl_matrix_set(R, 0, 1, 1.0);
    gsl_matrix_set(R, 0, 2, 1.0);
    gsl_matrix_set(R, 1, 0, 1.0);
    gsl_matrix_set(R, 1, 2, 1.0);
    gsl_matrix_set(R, 2, 0, 1.0);
    gsl_matrix_set(R, 2, 1, 1.0);
    double num_electrons = 8.0;
    int n0 = 8;
    bool use_cache = true;

    double tol = 1e-6;
    double E_Fermi = 0.0;

    double energy = BandEnergy(&E_Fermi, Hrs, R, num_electrons, n0, tol, use_cache);

    printf("energy = %f\n", energy);
    printf("E_Fermi = %f\n", E_Fermi);

    return 0;
}
