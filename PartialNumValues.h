#ifndef CWANNIER_PARTIALNUMVALUES_H
#define CWANNIER_PARTIALNUMVALUES_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "HTightBinding.h"
#include "PartialDosValues.h"
#include "ctetra/sum.h"

double** PartialNumValues(HTightBinding *Hrs, gsl_matrix *R, int na, int nb, int nc, double num_total_electrons, double **Es, double num_E, double *E_Fermi, double **num_states_Fermi);

#endif //CWANNIER_PARTIALNUMVALUES_H
