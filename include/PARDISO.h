#ifndef __PARDISO_H_
#define __PARDISO_H_


void pardiso(int *ia, int *ja, double *a, int nnz, int n, double *b, double *x){
    // ============================= CALL PARDISO ============================ //
    std::cout << "Setting up PARDISO parameters." << std::endl;

    //int nrhs = 1;          /* Number of right hand sides. */

    // must choose -2 for iterative solver
    int      mtype = -2;        /* Symmetric positive definite matrix */

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */
    void    *pt[64];

    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      maxfct, mnum, phase, error, msglvl, solver;

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    int      k;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */

    /* -------------------------------------------------------------------- */
    /* ..  Setup Pardiso control parameters.                                */
    /* -------------------------------------------------------------------- */
    std::cout << "Calling PARDISO init." << std::endl;

    error  = 0;
    solver = 0;

    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

    std::cout << "error : " << error << std::endl;

    if (error != 0)
    {
        if (error == -10 )
            printf("No license file found \n");
        if (error == -11 )
            printf("License is expired \n");
        if (error == -12 )
            printf("Wrong username or hostname \n");
        return 1;
    }

    std::cout << "[PARDISO]: License check was successful ... " << std::endl;

    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }

    iparm[2]  = num_procs;

    maxfct = 1;         /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */

    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */


    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */

    for (int i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (int i = 0; i < nnz; i++) {
        ja[i] += 1;
    }

    std::cout << "after 1 based conversion" << std::endl;


    /* -------------------------------------------------------------------- */
    /*  .. pardiso_chk_matrix(...)                                          */
    /*     Checks the consistency of the given matrix.                      */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */

    pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
    if (error != 0) {
        printf("\nERROR in consistency of matrix: %d", error);
        exit(1);
    }

    /* -------------------------------------------------------------------- */
    /* ..  pardiso_chkvec(...)                                              */
    /*     Checks the given vectors for infinite and NaN values             */
    /*     Input parameters (see PARDISO user manual for a description):    */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */

    pardiso_chkvec (&n, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR  in right hand side: %d", error);
        exit(1);
    }

    /* -------------------------------------------------------------------- */
    /* .. pardiso_printstats(...)                                           */
    /*    prints information on the matrix to STDOUT.                       */
    /*    Use this functionality only for debugging purposes                */
    /* -------------------------------------------------------------------- */

    pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR right hand side: %d", error);
        exit(1);
    }

    /* -------------------------------------------------------------------- */
    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    /*     all memory that is necessary for the factorization.              */
    /* -------------------------------------------------------------------- */

    std::cout << "After PARDISO checks." << std::endl;

    iparm[19-1] = -1; // in order to compute Gflops
    printf("\nGFlops factorisation : %i", iparm[19-1]);

    // start timer phase 1
    double timespent_p11 = -omp_get_wtime();

    phase = 11;

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, &idum, &nrhs,
                iparm, &msglvl, &ddum, &ddum, &error, dparm);

    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }

    // get time phase 1
    timespent_p11 += omp_get_wtime();
    // printf("\nTime spent on Phase 1 : %f", time_spent_p11);

    /* -------------------------------------------------------------------- */
    /* ..  Numerical factorization.                                         */
    /* -------------------------------------------------------------------- */

    // start timer phase 2
    double timespent_p22 = -omp_get_wtime();

    phase = 22;
    iparm[32] = 1; /* compute determinant */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs,
                iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }

    // get time phase 2
    timespent_p22 += omp_get_wtime();

    // printf("\nFactorization completed ...\n ");
    printf("\nFactorization completed .. \n");

    double log_det = dparm[32];
    printf("\nPardiso   log(det) = %f ", log_det);

    int gflops_fact = iparm[19-1];
    int mem_fact_solve = iparm[17-1];

    printf("\nGFlops factorisation : %i", iparm[19-1]);
    printf("\nMem fact + solve     : %i", mem_fact_solve);

    /* -------------------------------------------------------------------- */
    /* ..  Back substitution and iterative refinement.                      */
    /* -------------------------------------------------------------------- */

    // start timer phase 3
    //double timespent_p33 = 0;
    double timespent_p33 = -omp_get_wtime();

    phase = 33;

    iparm[7] = 0;       /* Max numbers of iterative refinement steps. */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs,
                iparm, &msglvl, b, x, &error,  dparm);

    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }

    // get time phase 3
    timespent_p33 += omp_get_wtime();

    printf("\nSolve completed ... ");

    /* -------------------------------------------------------------------- */
    /* ... Inverse factorization.                                           */
    /* -------------------------------------------------------------------- */

    double timespent_sel_inv = 0;

    if (solver == 0)
    {
    printf("\nCompute Diagonal Elements of the inverse of A ... \n");
    timespent_sel_inv = -omp_get_wtime();

    phase = -22;
    //iparm[35]  = 1; /*  no not overwrite internal factor L // crashes for larger matrices if uncommented */
    pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
            iparm, &msglvl, b, x, &error,  dparm);

    // get time to compute selected inverse
    timespent_sel_inv += omp_get_wtime();

    }

    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix back to 0-based C-notation.                       */
    /* -------------------------------------------------------------------- */
    for (int i = 0; i < n+1; i++) {
        ia[i] -= 1;
    }
    for (int i = 0; i < nnz; i++) {
        ja[i] -= 1;
    }


    /* -------------------------------------------------------------------- */
    /* ..  Print statistics                                                 */
    /* -------------------------------------------------------------------- */

    printf("\nTime spent on phase 1 : %f s", timespent_p11);
    printf("\nTime spent on phase 2 : %f s", timespent_p22);
    printf("\nTime spent on phase 3 : %f s", timespent_p33);
    printf("\nTime spent on sel inv : %f s\n", timespent_sel_inv);
}

#endif // __PARDISO_H_
