#include <stdio.h>
#include <math.h>
#include "HTightBinding.h"

int main(int argc, char *argv[]) {
    HTightBinding *Hrs = ExtractHTightBinding("test_data/Fe_soc/Fe_soc_hr.dat");
    int num_bands = Hrs->num_bands;
    if (num_bands != 18) {
        printf("Error: Incorrect band number extracted.");
        return 1;
    }
    gsl_matrix_complex *Hk = gsl_matrix_complex_alloc(num_bands, num_bands);
    double k0[3] = {0.0, 0.0, 0.0};
    HkRecip(Hrs, k0, Hk);

    double real_expected = 28.252214;
    double imag_expected = -2.312399946e-18;
    gsl_complex Hk_00 = gsl_matrix_complex_get(Hk, 0, 0);
    double eps = 1e-12;
    if (fabs(GSL_REAL(Hk_00) - real_expected) > eps) {
        printf("Error: incorrect real part of H_{k=0}[0, 0]\n");
        return 1;
    }
    if (fabs(GSL_IMAG(Hk_00) - imag_expected) > eps) {
        printf("Error: incorrect imaginary part of H_{k=0}[0, 0]\n");
        return 1;
    }
    printf("HTightBinding_test passed.\n");

    FreeHTightBinding(Hrs);
    return 0;
}
