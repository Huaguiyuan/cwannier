#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include "bstrlib/bstrlib.h"
#include "HTightBinding.h"
#include "DosValues.h"

gsl_matrix* parse_R_from_bs(char *b1, char *b2, char *b3);
int set_R_row(gsl_matrix *R, char *b, int row);
double* linspace(double a, double b, int num);
int write_dos_vals(char *outPath, int num_dos, double *Es, double *dos_vals);

int main(int argc, char *argv[]) {
    if (argc < 10) {
        printf("Usage: DosValues.out 'wannier_hr_path' 'out_path' 'start_energy' 'stop_energy' 'num_dos' 'num_k_per_dim' 'b_1_vec' 'b_2_vec' 'b_3_vec'\n");
        printf("Example: DosValues.out 'Fe_up_hr.dat' 'Fe_up_dos' '-5.0' '5.0' '500' '8' '1.0 0.0 0.0' '0.0 1.0 0.0' '0.0 0.0 1.0'\n");
        return 2;
    }
    // Parse arguments.
    char *hr_path = argv[1];
    char *out_path = argv[2];
    double start_energy = atof(argv[3]);
    double stop_energy = atof(argv[4]);
    int num_dos = atoi(argv[5]);
    int num_k_per_dim = atoi(argv[6]);
    gsl_matrix *R = parse_R_from_bs(argv[7], argv[8], argv[9]);

    // Set up data.
    HTightBinding *Hrs = ExtractHTightBinding(hr_path);
    if (Hrs == NULL) {
        printf("Error: failed to extract Hamiltonian\n");
        return 1;
    }
    double *Es = linspace(start_energy, stop_energy, num_dos);

    // Calculate DOS.
    double *dos_vals = DosValues(Hrs, R, num_k_per_dim, Es, num_dos);

    // Write out DOS.
    int write_err = write_dos_vals(out_path, num_dos, Es, dos_vals);
    if (write_err != 0) {
        printf("Error writing DOS values\n");
    }

    // Cleanup. 
    free(dos_vals);
    free(Es);
    gsl_matrix_free(R);
    return 0;
}

gsl_matrix* parse_R_from_bs(char *b1, char *b2, char *b3) {
    gsl_matrix *R = gsl_matrix_alloc(3, 3);

    int b1_err = set_R_row(R, b1, 0);
    int b2_err = set_R_row(R, b2, 1);
    int b3_err = set_R_row(R, b3, 2);
    if (b1_err != 0 || b2_err != 0 || b3_err != 0) {
        printf("ERROR: failed to parse reciprocal lattice vectors\n");
    }

    return R;
}

int set_R_row(gsl_matrix *R, char *b, int row) {
    bstring b_val_str = bfromcstr(b);
    if (b_val_str == NULL) {
        return 1;
    }
    struct bstrList *b_val_split = bsplit(b_val_str, ' ');

    int col = 0;
    int j;
    for (j = 0; j < b_val_split->qty; j++) {
        if (b_val_split->entry[j]->slen != 0) {
            char *entry_cstr = bstr2cstr(b_val_split->entry[j], 'N');
            double entry = atof(entry_cstr);
            gsl_matrix_set(R, row, col, entry);
            bcstrfree(entry_cstr);
            col++;
        }
    }
    
    bstrListDestroy(b_val_split);
    bdestroy(b_val_str);
    return 0;
}

double* linspace(double a, double b, int num) {
    double *vals = malloc(num * sizeof(double));
    double step = (b - a) / (num - 1);
    double x = a;
    int i;
    for (i = 0; i < num; i++) {
        vals[i] = x;
        x += step;
    }
    return vals;
}

int write_dos_vals(char *outPath, int num_dos, double *Es, double *dos_vals) {
    FILE *fp = fopen(outPath, "w");
    if (fp == NULL) {
        return 1;
    }

    // Write header.
    fprintf(fp, "E\tDOS\n");
    // Write DOS data.
    int i;
    for (i = 0; i < num_dos; i++) {
        fprintf(fp, "%.10f\t%.10f\n", Es[i], dos_vals[i]);
    }

    fclose(fp);
    return 0;
}
