/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2013.
 */

#include <fftw3.h>
#include "FFTW.h"

namespace edu {
namespace mit {
namespace fftw {

    FFTW_FFTWPlan FFTWWrapper::fftwPlanDft1d(int N, fftw_complex* input, fftw_complex* output, bool forward) {
        return fftw_plan_dft_1d(N, input, output, forward?-1:1, FFTW_ESTIMATE);
    }

    FFTW_FFTWPlan FFTWWrapper::fftwPlanDft1d(int N, double* input, fftw_complex* output) {
        return fftw_plan_dft_r2c_1d(N, input, output, FFTW_ESTIMATE);
    }

    FFTW_FFTWPlan FFTWWrapper::fftwPlanDft1d(int N, fftw_complex* input, double* output) {
        return fftw_plan_dft_c2r_1d(N, input, output, FFTW_ESTIMATE);
    }

    FFTW_FFTWPlan FFTWWrapper::fftwPlanDft1d(int N, int howmany, fftw_complex* input, fftw_complex* output, bool forward) {
        return fftw_plan_many_dft(1, &N, howmany,
                              input, &N,
                              1, N,
                              output, &N,
                              1, N,
                              forward?-1:1, FFTW_ESTIMATE);
    }

    FFTW_FFTWPlan FFTWWrapper::fftwPlanDft1d(int N, int howmany, fftw_complex* input, double* output) {
        return fftw_plan_many_dft_c2r(1, &N, howmany,
                              input, &N,
                              1, N,
                              output, &N,
                              1, N,
                              FFTW_ESTIMATE);
    }

    FFTW_FFTWPlan FFTWWrapper::fftwPlanDft1d(int N, int howmany, double* input, fftw_complex* output) {
        return fftw_plan_many_dft_r2c(1, &N, howmany,
                              input, &N,
                              1, N,
                              output, &N,
                              1, N,
                              FFTW_ESTIMATE);
    }
 
    FFTW_FFTWPlan FFTWWrapper::fftwPlanDft3d(int N1, int N2, int N3, fftw_complex* input, fftw_complex* output, bool forward) {
        return fftw_plan_dft_3d(N1, N2, N3, input, output, forward?-1:1, FFTW_ESTIMATE);
    }

    FFTW_FFTWPlan FFTWWrapper::fftwPlanDft3d(int N1, int N2, int N3, double* input, fftw_complex* output) {
        return fftw_plan_dft_r2c_3d(N1, N2, N3, input, output, FFTW_ESTIMATE);
    }

    FFTW_FFTWPlan FFTWWrapper::fftwPlanDft3d(int N1, int N2, int N3, fftw_complex* input, double* output) {
        return fftw_plan_dft_c2r_3d(N1, N2, N3, input, output, FFTW_ESTIMATE);
    }

    void FFTWWrapper::fftwExecute(FFTW_FFTWPlan p) {
        fftw_execute(p);
    }

    void FFTWWrapper::fftwDestroyPlan(FFTW_FFTWPlan p) {
        fftw_destroy_plan(p);
    }

    // TODO threaded FFTW.  This is not needed ATM.
    /*
    int FFTWWrapper::fftwInitThreads() {
        return fftw_init_threads();
    }

    void FFTWWrapper::fftwPlanWithNThreads(int nThreads) {
        fftw_plan_with_nthreads(nThreads);
    }

    void FFTWWrapper::fftwCleanupThreads() {
        fftw_cleanup_threads();
    }
    */
}  // namespace fftw
}  // namespace mit
}  // namespace edu

