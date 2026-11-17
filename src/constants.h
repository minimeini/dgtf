#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H

// Armadillo configuration macros (ARMA_*) only take effect if defined before
// <RcppArmadillo.h> is first included, e.g. through PKG_CPPFLAGS in
// src/Makevars. Every translation unit includes it before this header, so
// they are not set here.

// Use OpenMP only when the compiler enables it (-fopenmp defines _OPENMP).
// src/Makevars does not enable it, so the parallel sections run serially
// and <omp.h> is not needed.
#if defined(_OPENMP) && !defined(DGTF_USE_OPENMP)
    #define DGTF_USE_OPENMP
#endif

// #if !defined(DGTF_TIMING_HVA)
//     #define DGTF_TIMING_HVA
// #endif

// #if !defined(DGTF_TIMING_SMC)
//     #define DGTF_TIMING_SMC
// #endif


// #if !defined(DGTF_DO_BOUND_CHECK)
//     #define DGTF_DO_BOUND_CHECK
// #endif

// #if !defined(DGTF_DETAILED_OUTPUT)
//     #define DGTF_DETAILED_OUTPUT
// #endif


inline constexpr unsigned int NUM_THREADS = 8;

inline constexpr double EPS = 2.220446e-16;
inline constexpr double EPS8 = 1.e-8;
inline constexpr double UPBND = 700.;
inline constexpr double LOG2PI = 1.837877;

inline constexpr double covid_m = 4.7;
inline constexpr double covid_s = 2.9;
inline constexpr unsigned int MIN_LAG = 2;

inline constexpr double LN_MU = 1.386262;
inline constexpr double LN_SD2 = 0.3226017;
inline constexpr double NB_KAPPA = 0.395;
inline constexpr double NB_R = 6;

inline constexpr double NB_LAMBDA = 0.;
inline constexpr double NB_DELTA = 30;

inline constexpr bool VERBOSE = false;
inline constexpr unsigned int MAX_ITER = 100;
#endif
