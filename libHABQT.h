/*

cdecl calling convention is used.

This file documents the ABI of shared library. See README.md for details.

*/


extern void* tomography_init(int QBNUM, int PNUM, int MHMCITER, int POVMITER,
    int VERB);
/*
Initialise storage for tomography data and return a pointer to it
(MEMP metavar)

You can initialise storage as many times as necessary, it is persistent and
will exist untill freed by tomography_free. You can have several instances
active at the same time.

Arguments:
  QBNUM    -- Number of quantum bits under tomography
  PNUM     -- Number of particles to use for approximating the distribution
              over states (per rank)
  MHMCITER -- Number of Metropolis–Hastings steps to perform when resampling
              (after adjusting for acceptance rate)
  POVMITER -- Number of optimisation steps to perform when searching for
              optimal measurment
  VERB     -- Verbosity level of output to stdout from 0 (no output) to
              1 (full output)
*/


extern void tomography(void *MEMP, double *SVR, double *SVI, double *DMR,
    double *DMI, double *POVMR, double *POVMI);
/*
Run a tomography update after obtaining the measurement and return resulting
mean Bayesian estimate and optimal separable POVM to perform next.

Arguments:
  MEMP     -- Pointer to tomography state storage obtained from tomography_init
  SVR      -- Pointer to the first element of double array storing real
              part of state vector. Array must hold 2 ^ QBNUM elements.
              Contents will not be modified by the fuction.
  SVI      -- Same as SVR, but for the imaginary part.
  DMR      -- Pointer to the first element of double array into which the real
              part of the mean estimate density matrix will be written.
              Must be able to hold 2 ^ (2 * QBNUM) double elements.
  DMI      -- Same as DMR, but for the imaginary part.
  POVMR    -- Pointer to the first element of doulbe array into which the real
              part of POVM elements will be written. Must be able to hold
              2 ^ (2 * QBNUM) elements. POVM elements are stored sequentially
              as scaled pure state vectors.
  POVMI    -- Same as POVMR, but for imaginary parts.
*/


extern void tomography_free(void *MEMP);
/*
Free the internal storage. Behaviour on trying to use the pointer with
tomography function afterwards is undefined, unless it’s been reinitialised
with tomography_init.

Arguments:
  MEMP     -- Pointer to tomography state storage obtained from tomography_init
*/


/*
Example of usage:
  const int nq = 1; // Single qubit.
  const int dim = pow(2,nq); // Dimension of space.
  double dmr[dim][dim], dmi[dim][dim]; // DMR, DMI allocation.
  double povmr[dim][dim], povmi[dim][dim]; //POVMR, POVMI allocation.
  double svr[dim], svi[dim]; // SVR, SVI allocation.
  memset(svr, 0, dim*sizeof(double));
  memset(svi, 0, dim*sizeof(double));
  svr[0] = 1; // Measurement result is a projection on first basis vector.
  void *storage = tomography_init(nq,1000,50,0); // Storage allocation.
  tomography(storage,&svr[0],&svi[0],&dmr[0][0],&dmi[0][0],&povmr[0][0],
      &povmi[0][0]);
  //... Read and process values in dmr and dmi, povmr and povmi, perform a
  //mesurement (optionally using the POVM returned by the function) and write
  //results to svr and svi, then call the tomography fuction again:
  tomography(storage,&svr[0],&svi[0],&dmr[0][0],&dmi[0][0],&povmr[0][0],
      &povmi[0][0]);
  //... Process the new estimate in dmr and dmi and repeat calls to tomography
  //as necessary.  When finished, free the internal storage used by tomography.
  tomography_free(storage);
*/
