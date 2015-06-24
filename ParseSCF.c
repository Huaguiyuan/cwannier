#include "ParseSCF.h"

// Scan through the Quantum Espresso scf.out file at filePath and extract:
// the number of electrons; the lattice scale alat (in units of the Bohr radius);
// the reciprocal lattice vectors R[i, :] (in units of 2pi/alat).
int ParseSCF(char *filePath, double *num_electrons, double *alat, gsl_matrix *R) {
    FILE *fp = fopen(filePath, "r");
    if (fp == NULL) {
        return CWANNIER_PARSESCF_ERR;
    }
    struct bStream *infile = bsopen((bNread)fread, fp);

    // Need to initialize a bstring for bsreadln.
    bstring line = bfromcstr("init");
    if (line == NULL) {
        return CWANNIER_PARSESCF_ERR;
    }

    // Iterate over lines.
    // Look for: "lattice parameter"; "number of electrons"; "b(1)", "b(2)", "b(3)".
    bstring lpstr = bfromcstr("lattice parameter");
    bstring nelstr = bfromcstr("number of electrons");
    bstring b1 = bfromcstr("b(1)");
    bstring b2 = bfromcstr("b(2)");
    bstring b3 = bfromcstr("b(3)");
    while (true) {
        int err = bsreadln(line, infile, '\n');
        if (err == BSTR_ERR) {
            // Hit end of file.
            break;
        }
        btrimws(line);

        int index = binstr(line, 0, lpstr);
        if (index != BSTR_ERR) {
            // Found "lattice parameter"
            int substr_index = 0;
            struct bstrList *line_split = bsplit(line, ' ');
            int j;
            for (j = 0; j < line_split->qty; j++) {
                if (line_split->entry[j]->slen != 0) {
                    if (substr_index == 4) {
                        char *alat_cstr = bstr2cstr(line_split->entry[j], 'N');
                        *alat = atof(alat_cstr);
                        bcstrfree(alat_cstr);
                    }
                    substr_index++;
                }
            }
        }
        index = binstr(line, 0, nelstr);
        if (index != BSTR_ERR) {
            // Found "number of electrons"
            int substr_index = 0;
            struct bstrList *line_split = bsplit(line, ' ');
            int j;
            for (j = 0; j < line_split->qty; j++) {
                if (line_split->entry[j]->slen != 0) {
                    if (substr_index == 4) {
                        char *nel_cstr = bstr2cstr(line_split->entry[j], 'N');
                        *num_electrons = atof(nel_cstr);
                        bcstrfree(nel_cstr);
                    }
                    substr_index++;
                }
            }
        }
        index = binstr(line, 0, b1);
        if (index != BSTR_ERR) {
            // Found "b(1)"
            get_b_row(R, 0, line);
        }
        index = binstr(line, 0, b2);
        if (index != BSTR_ERR) {
            // Found "b(2)"
            get_b_row(R, 1, line);
        }
        index = binstr(line, 0, b3);
        if (index != BSTR_ERR) {
            // Found "b(3)"
            get_b_row(R, 2, line);
        }
    }
    bdestroy(lpstr);
    bdestroy(nelstr);
    bdestroy(b1);
    bdestroy(b2);
    bdestroy(b3);

    bdestroy(line);
    bsclose(infile);
    fclose(fp);
    return CWANNIER_PARSESCF_OK;
}

void get_b_row(gsl_matrix *R, int row, bstring line) {
    int substr_index = 0;
    struct bstrList *line_split = bsplit(line, ' ');
    int j;
    for (j = 0; j < line_split->qty; j++) {
        if (line_split->entry[j]->slen != 0) {
            if (substr_index == 3 || substr_index == 4 || substr_index == 5) {
                char *entry_cstr = bstr2cstr(line_split->entry[j], 'N');
                double entry = atof(entry_cstr);
                gsl_matrix_set(R, row, substr_index - 3, entry);
                bcstrfree(entry_cstr);
            }
            substr_index++;
        }
    }
}
