#ifndef CWANNIER_DOS_UTIL
#define CWANNIER_DOS_UTIL

#include <stdlib.h>
#include <stdio.h>
#include "bstrlib/bstrlib.h"
#include <gsl/gsl_matrix.h>

gsl_matrix* parse_R_from_bs(char *b1, char *b2, char *b3);
int set_R_row(gsl_matrix *R, char *b, int row);
double* linspace(double a, double b, int num);

#endif //CWANNIER_DOS_UTIL
