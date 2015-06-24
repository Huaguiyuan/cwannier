#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include "bstrlib/bstrlib.h"
#include "paths.h"
#include "ParseSCF.h"
#include "HTightBinding.h"
#include "BandEnergy.h"
#include "SpinOrbit.h"

int main(int argc, char *argv[]) {
    if (argc < 8) {
        printf("SOC-induced anisotropy calculation -- invoke with:\n");
        printf("anisotropy.out (system_name) (n) (soc_strength) (theta1) (phi1) (theta2) (phi2)\n");
        return 1;
    }
    char *system_name = argv[1];
    int n = atoi(argv[2]);
    double soc_strength = atof(argv[3]);
    double theta1 = atof(argv[4]);
    double phi1 = atof(argv[5]);
    double theta2 = atof(argv[6]);
    double phi2 = atof(argv[7]);
    
    char *scf_path = cwannier_data_path(system_name, "wannier", "\0", "scf.out");
    double num_electrons, alat;
    gsl_matrix *R = gsl_matrix_alloc(3, 3);

    int err = ParseSCF(scf_path, &num_electrons, &alat, R);
    if (err != CWANNIER_PARSESCF_OK) {
        printf("Error code = %d returned from ParseSCF\n", err);
        return err;
    }
    bcstrfree(scf_path);

    char *hr_up_path = cwannier_data_path(system_name, "wannier", system_name, "_up_hr.dat");
    HTightBinding *Hrs_up = ExtractHTightBinding(hr_up_path);
    bcstrfree(hr_up_path);

    char *hr_dn_path = cwannier_data_path(system_name, "wannier", system_name, "_dn_hr.dat");
    HTightBinding *Hrs_dn = ExtractHTightBinding(hr_dn_path);
    bcstrfree(hr_dn_path);

    printf("Hamiltonian loaded.\n");

    HTightBinding *Hrs_soc_1 = HamiltonianWithSOC(soc_strength, theta1, phi1, Hrs_up, Hrs_dn);
    HTightBinding *Hrs_soc_2 = HamiltonianWithSOC(soc_strength, theta2, phi2, Hrs_up, Hrs_dn);

    bool use_cache = true;
    double E_Fermi_1 = 0.0;
    double E_Fermi_2 = 0.0;

    double energy1 = BandEnergy(&E_Fermi_1, Hrs_soc_1, R, num_electrons, n, use_cache);
    printf("Got energy1 = %f\n", energy1);

    double energy2 = BandEnergy(&E_Fermi_2, Hrs_soc_2, R, num_electrons, n, use_cache);
    printf("Got energy2 = %f\n", energy2);

    printf("energy1 - energy2 = %e\n", energy1 - energy2);

    FreeHTightBinding(Hrs_up);
    FreeHTightBinding(Hrs_dn);
    FreeHTightBinding(Hrs_soc_1);
    FreeHTightBinding(Hrs_soc_2);
    return 0;
}
