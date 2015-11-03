#include "dos_util.h"

#include "HTightBinding.h"
#include "DosValues.h"

int write_dos_vals(char *outPath, int num_dos, double *Es, double *dos_vals);

int main(int argc, char *argv[]) {
    if (argc < 10) {
        printf("Usage: RunDosValues.out 'wannier_hr_path' 'out_path' 'start_energy' 'stop_energy' 'num_dos' 'num_k_per_dim' 'b_1_vec' 'b_2_vec' 'b_3_vec'\n");
        printf("Example: RunDosValues.out 'Fe_up_hr.dat' 'Fe_up_dos' '-5.0' '5.0' '500' '8' '1.0 0.0 0.0' '0.0 1.0 0.0' '0.0 0.0 1.0'\n");
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
