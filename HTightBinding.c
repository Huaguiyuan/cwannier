#include "HTightBinding.h"

/*
HTightBinding* ExtractHTightBinding(char *filePath) {
    HTightBinding* Hrs = (HTightBinding *)malloc(sizeof(HTightBinding));

    FILE *fp = fopen(filePath, "r");

    char *comment = (char *)malloc(500 * sizeof(char));
    fscanf("%s\n");

    double *ras = (double *)malloc(num_rs * sizeof(double));
    double *rbs = (double *)malloc(num_rs * sizeof(double));
    double *rcs = (double *)malloc(num_rs * sizeof(double));
    double *degen = (double *)malloc(num_rs * sizeof(double));

    gsl_matrix_complex **Hrs = (gsl_matrix_complex **)malloc(num_rs * sizeof(gsl_matrix_complex*));
}

void FreeHTightBinding(HTightBinding *Hrs) {
    free(Hrs->ras);
    free(Hrs->rbs);
    free(Hrs->rcs);
    free(Hrs->degens);

    int i;
    for (i = 0; i < Hrs->num_rs; i++) {
        gsl_matrix_complex_free(Hrs->Hrs[i])
    }
    free(Hrs->Hrs);

    free(Hrs);
}
*/

// Put the value H(k) derived from Hrs into Hk.
// When HkRecip() is called, Hk should be initialized with all zeros.
void HkRecip(HTightBinding *Hrs, double k[3], gsl_matrix_complex *Hk) {
    int nr = Hrs->num_rs;
    int nb = Hrs->num_bands;
    double ra, rb, rc;
    gsl_complex weight, ikr, mul, Hr_val, Hk_val, Hk_new;
    int i, row, col;
    gsl_matrix_complex *this_Hr;

    for (i = 0; i < nr; i++) {
        ra = Hrs->ras[i];
        rb = Hrs->rbs[i];
        rc = Hrs->rcs[i];
        weight = gsl_complex_rect(1.0/Hrs->degens[i], 0.0);
        ikr = gsl_complex_rect(0.0, 2.0*M_PI*(k[0]*ra + k[1]*rb + k[2]*rc));
        mul = gsl_complex_mul(weight, gsl_complex_exp(ikr));
        this_Hr = Hrs->Hrs[i];

        for (row = 0; row < nb; row++) {
            for (col = 0; col < nb; col++) {
                Hr_val = gsl_matrix_complex_get(this_Hr, row, col);
                Hk_val = gsl_matrix_complex_get(Hk, row, col);
                Hk_new = gsl_complex_add(Hk_val, gsl_complex_mul(mul, Hr_val));
                gsl_matrix_complex_set(Hk, row, col, Hk_new);
            }
        }
    }
}

gsl_matrix_complex* HrAtR(HTightBinding *Hrs, double R[3]) {
    double eps = 1e-12;
    double ra, rb, rc;
    int nr = Hrs->num_rs;
    int i;

    for (i = 0; i < nr; i++) {
        ra = Hrs->ras[i];
        rb = Hrs->rbs[i];
        rc = Hrs->rcs[i];
        if ((fabs(ra - R[0]) < eps) && (fabs(rb - R[1]) < eps) && (fabs(rc - R[2]) < eps)) {
            return Hrs->Hrs[i];
        }
    }
    return NULL;
}
