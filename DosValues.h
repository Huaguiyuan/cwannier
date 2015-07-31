#ifndef CWANNIER_DOS_VALUES_H
#define CWANNIER_DOS_VALUES_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "HTightBinding.h"
#include "ctetra/dos.h"

double* DosValues(HTightBinding *Hrs, gsl_matrix *R, int num_k_per_dim, double *Es, double num_dos);

#endif // CWANNIER_DOS_VALUES_H
