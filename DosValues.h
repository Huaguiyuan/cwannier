#ifndef CWANNIER_DOS_VALUES_H
#define CWANNIER_DOS_VALUES_H

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "HTightBinding.h"
#include "ctetra/dos.h"

double* DosValues(HTightBinding *Hrs, gsl_matrix *R, int num_k_per_dim, double *Es, double num_dos, bool all_Es);

double* DosEnergyDerivValues(HTightBinding *Hrs, gsl_matrix *R, int num_k_per_dim, double *Es, int num_dos, double num_electrons, double *fermi, double *dos_fermi, double *dos_deriv_fermi);

#endif // CWANNIER_DOS_VALUES_H
