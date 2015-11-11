#include "dos_util.h"

#include "HTightBinding.h"
#include "DosValues.h"

void print_usage();
int write_dos_vals(char *outPath, int num_dos, double *Es, double *dos_vals, int calc_deriv, double *dos_deriv_vals, double fermi, double dos_fermi, double dos_deriv_fermi);

void print_usage() {
    printf("Usage: RunDosValues.out 'wannier_hr_path' 'out_path' 'start_energy' 'stop_energy' 'num_dos' 'num_k_per_dim' 'b_1_vec' 'b_2_vec' 'b_3_vec' [calc_deriv] [num_electrons (required if calc_deriv == 1)]\n");
    printf("Example: RunDosValues.out 'Fe_up_hr.dat' 'Fe_up_dos' '-5.0' '5.0' '500' '8' '1.0 0.0 0.0' '0.0 1.0 0.0' '0.0 0.0 1.0'\n");
    printf("OR: RunDosValues.out 'Fe_up_hr.dat' 'Fe_up_dos' '-5.0' '5.0' '500' '8' '1.0 0.0 0.0' '0.0 1.0 0.0' '0.0 0.0 1.0' '1' '1.0'\n");
}

int main(int argc, char *argv[]) {
    if (argc < 10) {
        print_usage();
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
    int calc_deriv = 0;
    double num_electrons = 0.0;
    if (argc > 10) {
        calc_deriv = atoi(argv[10]);
        if (calc_deriv && (argc < 12)) {
            print_usage();
            return 2;
        }
        if (calc_deriv) {
            num_electrons = atof(argv[11]);
        }
    }

    // Set up data.
    HTightBinding *Hrs = ExtractHTightBinding(hr_path);
    if (Hrs == NULL) {
        printf("Error: failed to extract Hamiltonian\n");
        return 1;
    }
    double *Es = linspace(start_energy, stop_energy, num_dos);

    // Calculate DOS.
    double *dos_vals = NULL;
    double *dos_deriv_vals = NULL;
    double fermi, dos_fermi, dos_deriv_fermi;

    if (calc_deriv) {
        bool all_Es = true;
        dos_vals = DosValues(Hrs, R, num_k_per_dim, Es, num_dos, all_Es);
        dos_deriv_vals = DosEnergyDerivValues(Hrs, R, num_k_per_dim, Es, num_dos, num_electrons, &fermi, &dos_fermi, &dos_deriv_fermi);
    } else {
        bool all_Es = false;
        dos_vals = DosValues(Hrs, R, num_k_per_dim, Es, num_dos, all_Es);
    }

    // Write out DOS.
    int write_err = write_dos_vals(out_path, num_dos, Es, dos_vals, calc_deriv, dos_deriv_vals, fermi, dos_fermi, dos_deriv_fermi);
    if (write_err != 0) {
        printf("Error writing DOS values\n");
    }

    // Cleanup. 
    if (calc_deriv) {
        free(dos_deriv_vals);
    }
    free(dos_vals);
    free(Es);
    gsl_matrix_free(R);
    return 0;
}

int write_dos_vals(char *outPath, int num_dos, double *Es, double *dos_vals, int calc_deriv, double *dos_deriv_vals, double fermi, double dos_fermi, double dos_deriv_fermi) {
    FILE *fp = fopen(outPath, "w");
    if (fp == NULL) {
        return 1;
    }

    // Write Fermi energy header and data, if we have it.
    if (calc_deriv) {
        fprintf(fp, "E_Fermi\tDOS_Fermi\t(d/dE)DOS_Fermi\n");
        fprintf(fp, "%.10f\t%.10f\t%.10f\n", fermi, dos_fermi, dos_deriv_fermi);
    }

    // Write header.
    fprintf(fp, "E\tDOS");
    if (calc_deriv) {
        fprintf(fp, "\t(d/dE)DOS\n");
    } else {
        fprintf(fp, "\n");
    }
    // Write DOS data.
    int i;
    for (i = 0; i < num_dos; i++) {
        if (calc_deriv) {
            fprintf(fp, "%.10f\t%.10f\t%.10f\n", Es[i], dos_vals[i], dos_deriv_vals[i]);
        } else {
            fprintf(fp, "%.10f\t%.10f\n", Es[i], dos_vals[i]);
        }
    }

    fclose(fp);
    return 0;
}
