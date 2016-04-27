#ifndef CWANNIER_PARTIALDOS_VALUES_H
#define CWANNIER_PARTIALDOS_VALUES_H

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "HTightBinding.h"
#include "ctetra/partial.h"

double** PartialDosValues(HTightBinding *Hrs, gsl_matrix *R, int na, int nb, int nc, double sigma, double **Es, double num_dos);

void sort_evals_evecs(gsl_vector *energies, gsl_matrix_complex *evecs, int num_bands);

#endif // CWANNIER_PARTIALDOS_VALUES_H
