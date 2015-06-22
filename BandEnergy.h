#ifndef CWANNIER_BAND_ENERGY
#define CWANNIER_BAND_ENERGY

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

// TODO - include tetra

double BandEnergy(HTightBinding *Hrs, double R[3][3], double num_electrons, int n0, double tol);

#endif // CWANNIER_BAND_ENERGY
