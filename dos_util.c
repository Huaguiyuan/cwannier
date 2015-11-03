#include "dos_util.h"

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
