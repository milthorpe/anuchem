#include "FFTW.h"
#include <fftw3.h>

namespace edu {
    namespace mit {
        namespace fftw {
            int FFTWWrapper::fftwInitThreads() {
                return fftw_init_threads();
            }

            void FFTWWrapper::fftwPlanWithNThreads(int nThreads) {
                fftw_plan_with_nthreads(nThreads);
            }

            FFTW_FFTWPlan FFTWWrapper::fftwPlanDft1d(int N, fftw_complex* input, fftw_complex* output, bool forward) {
                return fftw_plan_dft_1d(N, input, output, forward?-1:1, FFTW_ESTIMATE);
            }
            
            FFTW_FFTWPlan FFTWWrapper::fftwPlanDft3d(int N1, int N2, int N3, fftw_complex* input, fftw_complex* output, bool forward) {
                return fftw_plan_dft_3d(N1, N2, N3, input, output, forward?-1:1, FFTW_ESTIMATE);
            }

            void FFTWWrapper::fftwExecute(FFTW_FFTWPlan p) {
                fftw_execute(p);
            }

            void FFTWWrapper::fftwDestroyPlan(FFTW_FFTWPlan p) {
                fftw_destroy_plan(p);
            }

            void FFTWWrapper::fftwCleanupThreads() {
                fftw_cleanup_threads();
            }
        }
    }
}


