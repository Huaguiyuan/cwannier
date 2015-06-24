#ifndef CWANNIER_PARSESCF_H
#define CWANNIER_PARSESCF_H

#include <stdbool.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include "bstrlib/bstrlib.h"

#define CWANNIER_PARSESCF_OK 0
#define CWANNIER_PARSESCF_ERR 1

int ParseSCF(char *filePath, double *num_electrons, double *alat, gsl_matrix *R);

void get_b_row(gsl_matrix *R, int row, bstring line);

#endif // CWANNIER_PARSESCF_H
