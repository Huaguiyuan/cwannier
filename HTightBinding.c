#include "HTightBinding.h"

// Put the value H(k) derived from Hrs into Hk.
// When HkRecip() is called, Hk should be initialized with all zeros.
void HkRecip(HTightBinding *Hrs, double k[3], gsl_matrix_complex *Hk) {
    int nr = Hrs->num_rs;
    int nb = Hrs->num_bands;
    double ra, rb, rc;
    gsl_complex weight, ikr, mul, Hr_val, Hk_val, Hk_new;
    int i, row, col;
    gsl_matrix_complex *this_Hr;

    for (i = 0; i < nr; i++) {
        ra = Hrs->ras[i];
        rb = Hrs->rbs[i];
        rc = Hrs->rcs[i];
        weight = gsl_complex_rect(1.0/Hrs->degens[i], 0.0);
        ikr = gsl_complex_rect(0.0, 2.0*M_PI*(k[0]*ra + k[1]*rb + k[2]*rc));
        mul = gsl_complex_mul(weight, gsl_complex_exp(ikr));
        this_Hr = Hrs->Hrs[i];

        for (row = 0; row < nb; row++) {
            for (col = 0; col < nb; col++) {
                Hr_val = gsl_matrix_complex_get(this_Hr, row, col);
                Hk_val = gsl_matrix_complex_get(Hk, row, col);
                Hk_new = gsl_complex_add(Hk_val, gsl_complex_mul(mul, Hr_val));
                gsl_matrix_complex_set(Hk, row, col, Hk_new);
            }
        }
    }
}

// Return the matrix H[R] stored in Hrs.
gsl_matrix_complex* HrAtR(HTightBinding *Hrs, double R[3]) {
    double eps = 1e-12;
    double ra, rb, rc;
    int nr = Hrs->num_rs;
    int i;

    for (i = 0; i < nr; i++) {
        ra = Hrs->ras[i];
        rb = Hrs->rbs[i];
        rc = Hrs->rcs[i];
        if ((fabs(ra - R[0]) < eps) && (fabs(rb - R[1]) < eps) && (fabs(rc - R[2]) < eps)) {
            return Hrs->Hrs[i];
        }
    }
    return NULL;
}

// Read the file at filePath, which is an Hr.dat file from Wannier90.
// Extract the data contained there into the returned *HTightBinding.
HTightBinding* ExtractHTightBinding(char *filePath) {
    HTightBinding* Hrs = (HTightBinding *)malloc(sizeof(HTightBinding));

    FILE *fp = fopen(filePath, "r");
    struct bStream *infile = bsopen((bNread)fread, fp);
    //infile = bread((bNread)fread, fp);

    // Need to initialize a bstring for bsreadln.
    bstring line = bfromcstr("init");
    if (line == NULL) {
        return NULL;
    }

    // Comment line.
    int err = bsreadln(line, infile, '\n');
    if (err == BSTR_ERR) {
        return NULL;
    }

    // Number of bands = num_bands.
    err = bsreadln(line, infile, '\n');
    if (err == BSTR_ERR) {
        return NULL;
    }
    btrimws(line);
    char *line_cstr = bstr2cstr(line, 'N');
    int num_bands = atoi(line_cstr);
    bcstrfree(line_cstr);
    Hrs->num_bands = num_bands;

    // Number of ri - rj values = num_rs.
    bsreadln(line, infile, '\n');
    if (err == BSTR_ERR) {
        return NULL;
    }
    btrimws(line);
    line_cstr = bstr2cstr(line, 'N');
    int num_rs = atoi(line_cstr);
    bcstrfree(line_cstr);
    Hrs->num_rs = num_rs;

    // Allocate arrays of length num_rs.
    double *ras = (double *)malloc(num_rs * sizeof(double));
    double *rbs = (double *)malloc(num_rs * sizeof(double));
    double *rcs = (double *)malloc(num_rs * sizeof(double));
    double *degens = (double *)malloc(num_rs * sizeof(double));
    gsl_matrix_complex **Hr_values = (gsl_matrix_complex **)malloc(num_rs * sizeof(gsl_matrix_complex*));
    Hrs->ras = ras;
    Hrs->rbs = rbs;
    Hrs->rcs = rcs;
    Hrs->degens = degens;
    Hrs->Hrs = Hr_values;

    // Degeneracies.
    long int num_degen_lines = lrint(ceil(((double)num_rs) / 15.0));
    long int degen_line;
    int j, r_index = 0;
    for (degen_line = 0; degen_line < num_degen_lines; degen_line++) {
        bsreadln(line, infile, '\n');
        if (err == BSTR_ERR) {
            return NULL;
        }
        btrimws(line);
        struct bstrList *line_split = bsplit(line, ' ');
        for (j = 0; j < line_split->qty; j++) {
            if (line_split->entry[j]->slen != 0) {
                char *degen_cstr = bstr2cstr(line_split->entry[j], 'N');
                double degen = atof(degen_cstr);
                bcstrfree(degen_cstr);
                degens[r_index] = degen;
                r_index++;
            }
        }
        bstrListDestroy(line_split);
    }

    // Matrix values.
    long int num_matrix_lines = num_rs * num_bands * num_bands;
    long int matrix_line;
    int last_ra = 0;
    int last_rb = 0;
    int last_rc = 0;
    int ra = 0;
    int rb = 0;
    int rc = 0;
    int row = 0;
    int col = 0;
    double realpart, imagpart;
    r_index = -1;
    for (matrix_line = 0; matrix_line < num_matrix_lines; matrix_line++) {
        bsreadln(line, infile, '\n');
        if (err == BSTR_ERR) {
            return NULL;
        }
        btrimws(line);
        struct bstrList *line_split = bsplit(line, ' ');
        int val_index = 0;
        for (j = 0; j < line_split->qty; j++) {
            if (line_split->entry[j]->slen != 0) {
                char *val_cstr = bstr2cstr(line_split->entry[j], 'N');
                if (val_index == 0) {
                    ra = atoi(val_cstr);
                } else if (val_index == 1) {
                    rb = atoi(val_cstr);
                } else if (val_index == 2) {
                    rc = atoi(val_cstr);
                } else if (val_index == 3) {
                    row = atoi(val_cstr);
                } else if (val_index == 4) {
                    col = atoi(val_cstr);
                } else if (val_index == 5) {
                    realpart = atof(val_cstr);
                } else if (val_index == 6) {
                    imagpart = atof(val_cstr);
                    // Finished collecting parts of the line; add it to Hrs.
                    // Do we have a new value of (ra, rb, rc)?
                    if (matrix_line == 0 || (ra != last_ra || rb != last_rb || rc != last_rc)) {
                        r_index++;
                        ras[r_index] = (double)ra;
                        rbs[r_index] = (double)rb;
                        rcs[r_index] = (double)rc;
                        Hr_values[r_index] = gsl_matrix_complex_alloc(num_bands, num_bands);
                    }
                    gsl_complex complex_val = gsl_complex_rect(realpart, imagpart);
                    gsl_matrix_complex_set(Hr_values[r_index], row-1, col-1, complex_val);
                    last_ra = ra;
                    last_rb = rb;
                    last_rc = rc;
                }
                bcstrfree(val_cstr);
                val_index++;
            }
        }
        bstrListDestroy(line_split);
    }

    bsclose(infile);
    return Hrs;
}

void FreeHTightBinding(HTightBinding *Hrs) {
    free(Hrs->ras);
    free(Hrs->rbs);
    free(Hrs->rcs);
    free(Hrs->degens);

    int i;
    for (i = 0; i < Hrs->num_rs; i++) {
        gsl_matrix_complex_free(Hrs->Hrs[i]);
    }
    free(Hrs->Hrs);

    free(Hrs);
}


