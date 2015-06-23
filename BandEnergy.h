#ifndef CWANNIER_BAND_ENERGY
#define CWANNIER_BAND_ENERGY

#include <stdbool.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "HTightBinding.h"
#include "ctetra/sum.h"

double BandEnergy(double *E_Fermi, HTightBinding *Hrs, gsl_matrix *R, double num_electrons, int n0, double tol, bool use_cache);

#endif // CWANNIER_BAND_ENERGY
