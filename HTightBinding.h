#ifndef CWANNIER_HTIGHTBINDING_H
#define CWANNIER_HTIGHTBINDING_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "bstrlib/bstrlib.h"

typedef struct {
    int num_rs, num_bands;
    double *ras;
    double *rbs;
    double *rcs;
    double *degens;
    gsl_matrix_complex **Hrs;
} HTightBinding;

HTightBinding* ExtractHTightBinding(char *filePath);

void FreeHTightBinding(HTightBinding *Hrs);

void HkRecip(HTightBinding *Hrs, double k[3], gsl_matrix_complex *Hk);

gsl_matrix_complex* HrAtR(HTightBinding *Hrs, double R[3]);

#endif // CWANNIER_HTIGHTBINDING_H
