#ifndef CWANNIER_SPINORBIT_H
#define CWANNIER_SPINORBIT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include "HTightBinding.h"

typedef struct {
    gsl_complex val;
    int row;
    int col;
} socElem;

HTightBinding* HamiltonianWithSOC(double soc_strength, double theta, double phi, HTightBinding *Hrs_up, HTightBinding *Hrs_dn);

void addSOCElems(gsl_matrix_complex *Hr_onsite, double soc_strength, double theta, double phi);

gsl_matrix_complex* onSiteSOC_SpinZ();

void setWithHermitianConjugate(gsl_matrix_complex *M, int row, int col, gsl_complex val);

gsl_matrix_complex* onSiteSOC(double theta, double phi);

gsl_matrix_complex *RotationMatrix(double theta, double phi);
#endif // CWANNIER_SPINORBIT_H
