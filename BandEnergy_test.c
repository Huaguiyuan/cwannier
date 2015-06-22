#include "BandEnergy.h"

int main(int argc, char *argv[]) {
    HTightBinding *Hrs = ExtractHTightBinding("test_data/Fe_soc/Fe_soc_hr.dat");

    double R[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double num_electrons = 1.0;
    int n0 = 8;
    double tol = 1e-6;
    double E_Fermi = 0.0;
    double energy = BandEnergy(&E_Fermi, Hrs, R, num_electrons, n0, tol);

    printf("energy = %f\n", energy);
    printf("E_Fermi = %f\n", E_Fermi);

    return 0;
}
