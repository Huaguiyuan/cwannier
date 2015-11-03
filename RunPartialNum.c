#include "dos_util.h"

#include "HTightBinding.h"
#include "PartialNumValues.h"

int write_pnum_vals(char *outPath, int num_bands, int num_E, double *Es, double **num_vals, double E_Fermi, double *num_Fermi_vals);

int main(int argc, char *argv[]) {
    if (argc < 11) {
        printf("Usage: RunPartialNum.out 'wannier_hr_path' 'out_path' 'start_energy' 'stop_energy' 'num_E' 'num_electrons' 'num_k_per_dim' 'b_1_vec' 'b_2_vec' 'b_3_vec'\n");
        printf("Example: RunPartialNum.out 'Fe_up_hr.dat' 'Fe_up_dos' '-5.0' '5.0' '500' '1.0' '8' '1.0 0.0 0.0' '0.0 1.0 0.0' '0.0 0.0 1.0'\n");
        return 2;
    }
    // Parse arguments.
    char *hr_path = argv[1];
    char *out_path = argv[2];
    double start_energy = atof(argv[3]);
    double stop_energy = atof(argv[4]);
    int num_E = atoi(argv[5]);
    double num_electrons = atof(argv[6]);
    int num_k_per_dim = atoi(argv[7]);
    gsl_matrix *R = parse_R_from_bs(argv[8], argv[9], argv[10]);

    // Set up data.
    HTightBinding *Hrs = ExtractHTightBinding(hr_path);
    if (Hrs == NULL) {
        printf("Error: failed to extract Hamiltonian\n");
        return 1;
    }
    int num_bands = Hrs->num_bands;
    double *Es;
    double E_Fermi;
    double *num_Fermi_vals;

    printf("about to calculate num\n");
    // Calculate n(E).
    double **num_vals = PartialNumValues(Hrs, R, num_k_per_dim, num_electrons, &Es, num_E, &E_Fermi, &num_Fermi_vals);
    printf("num calculation finished\n");

    // Write out n(E).
    int write_err = write_pnum_vals(out_path, num_bands, num_E, Es, num_vals, E_Fermi, num_Fermi_vals);
    if (write_err != 0) {
        printf("Error writing num values\n");
    }

    // Cleanup.
    int band_index;
    for (band_index = 0; band_index < num_bands; band_index++) {
        free(num_vals[band_index]);
    }
    free(num_vals);
    free(num_Fermi_vals);
    free(Es);
    gsl_matrix_free(R);
    return 0;
}

int write_pnum_vals(char *outPath, int num_bands, int num_E, double *Es, double **num_vals, double E_Fermi, double *num_Fermi_vals) {
    int band_index, E_index;

    FILE *fp = fopen(outPath, "w");
    if (fp == NULL) {
        return 1;
    }

    // Write Fermi header.
    fprintf(fp, "E_Fermi");
    for (band_index = 0; band_index < num_bands; band_index++) {
        fprintf(fp, "\tN_Fermi(%d)", band_index);
    }
    fprintf(fp, "\n");

    // Write num at Fermi data (tab-separated).
    fprintf(fp, "%.10f", E_Fermi);
    for (band_index = 0; band_index < num_bands; band_index++) {
        fprintf(fp, "\t%.10f", num_Fermi_vals[band_index]);
    }
    fprintf(fp, "\n");

    // Write header.
    fprintf(fp, "E");
    for (band_index = 0; band_index < num_bands; band_index++) {
        fprintf(fp, "\tN(%d)", band_index);
    }
    fprintf(fp, "\n");

    // Write num data (tab-separated).
    for (E_index = 0; E_index < num_E; E_index++) {
        fprintf(fp, "%.10f", Es[E_index]);
        for (band_index = 0; band_index < num_bands; band_index++) {
            fprintf(fp, "\t%.10f", num_vals[band_index][E_index]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    return 0;
}
