#ifndef EDU_MIT_FFTW_FFTW_H
#define EDU_MIT_FFTW_FFTW_H

#include <fftw3.h>

namespace edu {
    namespace mit {
        namespace fftw {
            typedef fftw_plan FFTW_FFTWPlan;

            class FFTWWrapper {
                public:
                    static int fftwInitThreads();
                    static void fftwPlanWithNThreads(int nthreads);
                    static FFTW_FFTWPlan fftwPlanDft1d(int, fftw_complex*, fftw_complex*, bool);
                    static FFTW_FFTWPlan fftwPlanDft1d(int, double*, fftw_complex*);
                    static FFTW_FFTWPlan fftwPlanDft1d(int, fftw_complex*, double*);
                    static FFTW_FFTWPlan fftwPlanDft1d(int N, int howmany, fftw_complex* input, fftw_complex* output, bool forward);
                    static FFTW_FFTWPlan fftwPlanDft1d(int N, int howmany, double* input, fftw_complex* output);
                    static FFTW_FFTWPlan fftwPlanDft1d(int N, int howmany, fftw_complex* input, double* output);
                    static FFTW_FFTWPlan fftwPlanDft3d(int N1, int N2, int N3, fftw_complex* input, fftw_complex* output, bool forward);
                    static FFTW_FFTWPlan fftwPlanDft3d(int N1, int N2, int N3, fftw_complex* input, double* output);
                    static FFTW_FFTWPlan fftwPlanDft3d(int N1, int N2, int N3, double* input, fftw_complex* output);
                    static void fftwExecute(FFTW_FFTWPlan);
                    static void fftwDestroyPlan(FFTW_FFTWPlan);
                    static void fftwCleanupThreads();
            };
        }
    }
}
#endif


