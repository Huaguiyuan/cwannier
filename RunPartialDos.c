#include "dos_util.h"

#include "HTightBinding.h"
#include "PartialDosValues.h"

int write_pdos_vals(char *outPath, int num_bands, int num_dos, double *Es, double **dos_vals);

int main(int argc, char *argv[]) {
    if (argc < 11) {
        printf("Usage: RunPartialDos.out 'wannier_hr_path' 'out_path' 'start_energy' 'stop_energy' 'num_dos' 'sigma' 'num_k_per_dim' 'b_1_vec' 'b_2_vec' 'b_3_vec'\n");
        printf("Example: RunPartialDos.out 'Fe_up_hr.dat' 'Fe_up_dos' '-5.0' '5.0' '500' '0.001' '8' '1.0 0.0 0.0' '0.0 1.0 0.0' '0.0 0.0 1.0'\n");
        return 2;
    }
    // Parse arguments.
    char *hr_path = argv[1];
    char *out_path = argv[2];
    double start_energy = atof(argv[3]);
    double stop_energy = atof(argv[4]);
    int num_dos = atoi(argv[5]);
    double sigma = atof(argv[6]);
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

    printf("about to calculate DOS\n");
    // Calculate DOS.
    double **dos_vals = PartialDosValues(Hrs, R, num_k_per_dim, sigma, &Es, num_dos);
    printf("DOS calculation finished\n");

    // Write out DOS.
    int write_err = write_pdos_vals(out_path, num_bands, num_dos, Es, dos_vals);
    if (write_err != 0) {
        printf("Error writing DOS values\n");
    }

    // Cleanup.
    int band_index;
    for (band_index = 0; band_index < num_bands; band_index++) {
        free(dos_vals[band_index]);
    }
    free(dos_vals);
    free(Es);
    gsl_matrix_free(R);
    return 0;
}

int write_pdos_vals(char *outPath, int num_bands, int num_dos, double *Es, double **dos_vals) {
    int band_index, E_index;

    FILE *fp = fopen(outPath, "w");
    if (fp == NULL) {
        return 1;
    }

    // Write header.
    fprintf(fp, "E");
    for (band_index = 0; band_index < num_bands; band_index++) {
        fprintf(fp, "\tDOS(%d)", band_index);
    }
    fprintf(fp, "\n");

    // Write DOS data (tab-separated).
    for (E_index = 0; E_index < num_dos; E_index++) {
        fprintf(fp, "%.10f", Es[E_index]);
        for (band_index = 0; band_index < num_bands; band_index++) {
            fprintf(fp, "\t%.10f", dos_vals[band_index][E_index]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    return 0;
}
