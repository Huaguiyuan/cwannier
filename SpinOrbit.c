#include "SpinOrbit.h"

HTightBinding* HamiltonianWithSOC(double soc_strength, double theta, double phi, HTightBinding *Hrs_up, HTightBinding *Hrs_dn) {
    // bands per spin
    int N = Hrs_up->num_bands;

    int num_rs = Hrs_up->num_rs;

    double R0[3] = {0.0, 0.0, 0.0};
    double degen = 0.0;
    gsl_matrix_complex *Hr_up_onsite = HrAtR(Hrs_up, R0, &degen);
    if (Hr_up_onsite == NULL) {
        printf("Error: failed to retrieve on-site Hamiltonian.");
        exit(EXIT_FAILURE);
    }
    gsl_matrix_complex *Hr_dn_onsite = HrAtR(Hrs_dn, R0, &degen);
    if (Hr_dn_onsite == NULL) {
        printf("Error: failed to retrieve on-site Hamiltonian.");
        exit(EXIT_FAILURE);
    }

    // Initialize on-site part of spin-orbit coupled H with up-up and down-down blocks
    // given by Hr_up_onsite and Hr_dn_onsite; set up-down and down-up blocks to 0.
    gsl_matrix_complex *Hr_soc_onsite = gsl_matrix_complex_calloc(2*N, 2*N);
    int row, col;
    for (row = 0; row < N; row++) {
        for (col = 0; col < N; col++) {
            gsl_matrix_complex_set(Hr_soc_onsite, row, col, gsl_matrix_complex_get(Hr_up_onsite, row, col));
            gsl_matrix_complex_set(Hr_soc_onsite, row+N, col+N, gsl_matrix_complex_get(Hr_dn_onsite, row, col));
        }
    }
    addSOCElems(Hr_soc_onsite, soc_strength, theta, phi);

    // Construct HTightBinding with SOC contribution on the on-site term.
    HTightBinding *Hrs_soc = (HTightBinding*)malloc(sizeof(HTightBinding));
    Hrs_soc->num_bands = 2*N;
    Hrs_soc->num_rs = num_rs;

    Hrs_soc->ras = (double*)malloc(num_rs * sizeof(double));
    Hrs_soc->rbs = (double*)malloc(num_rs * sizeof(double));
    Hrs_soc->rcs = (double*)malloc(num_rs * sizeof(double));
    Hrs_soc->degens = (double*)malloc(num_rs * sizeof(double));
    Hrs_soc->Hrs = (gsl_matrix_complex**)malloc(num_rs * sizeof(gsl_matrix_complex*));

    int i;
    double eps = 1e-12;
    for (i = 0; i < num_rs; i++) {
        double ra = Hrs_up->ras[i];
        double rb = Hrs_up->rbs[i];
        double rc = Hrs_up->rcs[i];
        double degen = Hrs_up->degens[i];
        Hrs_soc->ras[i] = ra;
        Hrs_soc->rbs[i] = rb;
        Hrs_soc->rcs[i] = rc;
        Hrs_soc->degens[i] = degen;

        // On-site part?
        if (fabs(ra) < eps && fabs(rb) < eps && fabs(rc) < eps) {
            Hrs_soc->Hrs[i] = Hr_soc_onsite;
        } else {
            // Not on-site: copy up-up block from Hrs_up and down-down block from Hrs_dn;
            // set up-down and down-up blocks to 0.
            Hrs_soc->Hrs[i] = gsl_matrix_complex_calloc(2*N, 2*N);
            for (row = 0; row < N; row++) {
                for (col = 0; col < N; col++) {
                    gsl_matrix_complex_set(Hrs_soc->Hrs[i], row, col, gsl_matrix_complex_get(Hrs_up->Hrs[i], row, col));
                    gsl_matrix_complex_set(Hrs_soc->Hrs[i], row+N, col+N, gsl_matrix_complex_get(Hrs_dn->Hrs[i], row, col));
                }
            }
        }
    }
    return Hrs_soc;
}

void addSOCElems(gsl_matrix_complex *Hr_onsite, double soc_strength, double theta, double phi) {
    // number of bands per spin
    int N = Hr_onsite->size1 / 2;
    // number of atoms
    int N_atom = N / 9;

    // Spin-orbit coupling part of the Hamiltonian.
    // The value returned from onSiteSOC is independent of the SOC strength:
    // need to multiply to get the value we will use.
    gsl_matrix_complex *H_soc = onSiteSOC(theta, phi);
    gsl_matrix_complex_scale(H_soc, gsl_complex_rect(soc_strength, 0.0));

    int i, s, sp, row, col;
    // Iterate over atoms.
    for (i = 0; i < N_atom; i++) {
        // Iterate over spins.
        for (s = 0; s < 2; s++) {
            for (sp = 0; sp < 2; sp++) {
                int offset_s = 9*i + 9*N_atom*s;
                int offset_sp = 9*i + 9*N_atom*sp;
                // Set corresponding blocks in Hr_onsite and H_soc:
                // H_onsite[offset_s:offset_s+9, offset_sp:offset_sp+9] += H_soc[s*9:s*9+9, sp*9:sp*9+9]
                for (row = 0; row < 9; row++) {
                    for (col = 0; col < 9; col++) {
                        gsl_complex soc_val = gsl_matrix_complex_get(H_soc, s*9 + row, sp*9 + col);
                        gsl_complex orig_val = gsl_matrix_complex_get(Hr_onsite, offset_s + row, offset_sp + col);
                        gsl_complex new_val = gsl_complex_add(orig_val, soc_val);
                        gsl_matrix_complex_set(Hr_onsite, offset_s + row, offset_sp + col, new_val);
                    }
                }
            }
        }
    }
    gsl_matrix_complex_free(H_soc);
}

gsl_matrix_complex* onSiteSOC_SpinZ() {
    // Spin offsets = (0, 0)
    socElem up_up[3] = {{gsl_complex_rect(0.0, -0.5), 2, 3}, {gsl_complex_rect(0.0, 0.5), 6, 5},
            {gsl_complex_rect(0.0, 1.0), 8, 7}};

    // Spin offsets = (9, 9)
    socElem dn_dn[3] = {{gsl_complex_rect(0.0, -0.5), 6, 5}, {gsl_complex_rect(0.0, 0.5), 2, 3},
            {gsl_complex_rect(0.0, -1.0), 8, 7}};

    // Spin offsets = (0, 9)
    socElem up_dn[16] = {{gsl_complex_rect(-0.5, 0.0), 1, 2}, {gsl_complex_rect(-0.5, 0.0), 6, 8},
            {gsl_complex_rect(-0.5, 0.0), 5, 7}, {gsl_complex_rect(0.0, -0.5), 3, 1},
            {gsl_complex_rect(0.0, -0.5), 8, 5}, {gsl_complex_rect(0.0, -0.5), 6, 7},
            {gsl_complex_rect(0.5, 0.0), 2, 1}, {gsl_complex_rect(0.5, 0.0), 8, 6},
            {gsl_complex_rect(0.0, 0.5), 7, 6}, {gsl_complex_rect(0.0, 0.5), 1, 3},
            {gsl_complex_rect(0.0, 0.5), 5, 8}, {gsl_complex_rect(0.5, 0.0), 7, 5},
            {gsl_complex_rect(0.0, -sqrt(3)/2.0), 6, 4}, {gsl_complex_rect(-sqrt(3)/2.0, 0.0), 4, 5},
            {gsl_complex_rect(0.0, sqrt(3)/2.0), 4, 6}, {gsl_complex_rect(sqrt(3)/2.0, 0.0), 5, 4}};

    // Initialize socMatrix to 0, then set specified elements.
    gsl_matrix_complex *socMatrix = gsl_matrix_complex_calloc(18, 18);
    int i, row, col;
    socElem elem;
    for (i = 0; i < 3; i++) {
        elem = up_up[i];
        row = elem.row;
        col = elem.col;
        setWithHermitianConjugate(socMatrix, row, col, elem.val);
    }
    for (i = 0; i < 3; i++) {
        elem = dn_dn[i];
        row = elem.row + 9;
        col = elem.col + 9;
        setWithHermitianConjugate(socMatrix, row, col, elem.val);
    }
    for (i = 0; i < 16; i++) {
        elem = up_dn[i];
        row = elem.row;
        col = elem.col + 9;
        setWithHermitianConjugate(socMatrix, row, col, elem.val);
    }
    return socMatrix;
}

void setWithHermitianConjugate(gsl_matrix_complex *M, int row, int col, gsl_complex val) {
    gsl_matrix_complex_set(M, row, col, val);
    gsl_matrix_complex_set(M, col, row, gsl_complex_conjugate(val));
}

gsl_matrix_complex* onSiteSOC(double theta, double phi) {
    gsl_matrix_complex *H_soc_z = onSiteSOC_SpinZ();
    gsl_matrix_complex *R = RotationMatrix(theta, phi);

    // Initialize H_soc to 0, then set elements corresponding to those in H_soc_z.
    gsl_matrix_complex *H_soc = gsl_matrix_complex_calloc(18, 18);
    int l, lp, s, sp;
    // Iterate over angular momentum states.
    for (l = 0; l < 9; l++) {
        for (lp = 0; lp < 9; lp++) {
            // Iterate over spins.
            for (s = 0; s < 2; s++) {
                for (sp = 0; sp < 2; sp++) {
                    // up, up
                    gsl_complex Rval1 = gsl_complex_mul(gsl_complex_conjugate(gsl_matrix_complex_get(R, 0, s)), gsl_matrix_complex_get(R, 0, sp));
                    gsl_complex c1 = gsl_complex_mul(Rval1, gsl_matrix_complex_get(H_soc_z, l, lp));
                    // down, up
                    gsl_complex Rval2 = gsl_complex_mul(gsl_complex_conjugate(gsl_matrix_complex_get(R, 1, s)), gsl_matrix_complex_get(R, 0, sp));
                    gsl_complex c2 = gsl_complex_mul(Rval2, gsl_matrix_complex_get(H_soc_z, l+9, lp));
                    // up, down
                    gsl_complex Rval3 = gsl_complex_mul(gsl_complex_conjugate(gsl_matrix_complex_get(R, 0, s)), gsl_matrix_complex_get(R, 1, sp));
                    gsl_complex c3 = gsl_complex_mul(Rval3, gsl_matrix_complex_get(H_soc_z, l, lp+9));
                    // down, down
                    gsl_complex Rval4 = gsl_complex_mul(gsl_complex_conjugate(gsl_matrix_complex_get(R, 1, s)), gsl_matrix_complex_get(R, 1, sp));
                    gsl_complex c4 = gsl_complex_mul(Rval4, gsl_matrix_complex_get(H_soc_z, l+9, lp+9));

                    gsl_complex total = gsl_complex_add(c1, gsl_complex_add(c2, gsl_complex_add(c3, c4)));
                    gsl_matrix_complex_set(H_soc, l + s*9, lp + sp*9, total);
                }
            }
        }
    }
    gsl_matrix_complex_free(R);
    gsl_matrix_complex_free(H_soc_z);
    return H_soc;
}

gsl_matrix_complex* RotationMatrix(double theta, double phi) {
    double ct = cos(theta / 2.0);
    double st = sin(theta / 2.0);
    gsl_complex ip2 = gsl_complex_rect(0.0, phi / 2.0);
    gsl_complex mip2 = gsl_complex_rect(0.0, -phi / 2.0);
    gsl_matrix_complex *R = gsl_matrix_complex_calloc(2, 2);

    gsl_matrix_complex_set(R, 0, 0, gsl_complex_mul_real(gsl_complex_exp(mip2), ct));
    gsl_matrix_complex_set(R, 0, 1, gsl_complex_mul_real(gsl_complex_exp(mip2), -st));
    gsl_matrix_complex_set(R, 1, 0, gsl_complex_mul_real(gsl_complex_exp(ip2), st));
    gsl_matrix_complex_set(R, 1, 1, gsl_complex_mul_real(gsl_complex_exp(ip2), ct));

    return R;
}
