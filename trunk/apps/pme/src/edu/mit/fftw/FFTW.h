#ifndef EDU_MIT_FFTW_FFTW_H
#define EDU_MIT_FFTW_FFTW_H

#include <fftw3.h>

namespace edu {
    namespace mit {
        namespace fftw {
            typedef fftw_plan FFTW_FFTWPlan;

            class FFTWWrapper {
                public:
                   static FFTW_FFTWPlan fftwPlanDft1d(int, fftw_complex*, fftw_complex*, bool);
                   static FFTW_FFTWPlan fftwPlanDft3d(int N1, int N2, int N3, fftw_complex* input, fftw_complex* output, bool forward);
                   static void fftwExecute(FFTW_FFTWPlan);
                   static void fftwDestroyPlan(FFTW_FFTWPlan);
            };
        }
    }
}
#endif


