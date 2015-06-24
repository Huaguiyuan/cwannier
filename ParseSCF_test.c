#include "ParseSCF.h"

int main(int argc, char *argv[]) {
    double num_electrons, alat;
    gsl_matrix *R = gsl_matrix_alloc(3, 3);

    char *filePath = "test_data/Fe_sp/wannier/scf.out";
    int err = ParseSCF(filePath, &num_electrons, &alat, R);
    if (err != CWANNIER_PARSESCF_OK) {
        printf("Error code = %d returned from ParseSCF\n", err);
        return err;
    }

    printf("num_electrons = %f, alat = %f\n", num_electrons, alat);
    printf("R:\n");
    gsl_matrix_fprintf(stdout, R, "%f");
    return 0;
}
