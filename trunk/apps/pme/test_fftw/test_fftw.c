#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
int main() {
    fftw_complex *in, *out, *roundtrip;
    fftw_plan p;
    int N = 50;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    for (int i=0; i<N; i++) {
        in[i] = 1.0 / (i+1) * (1 + I);
        printf("%g + %gi ", creal(in[i]), cimag(in[i]));
    }
    fftw_execute(p);
    printf("result for %d:\n", N);
    for (int i=0; i<N; i++) {
        printf("%g + %gi ", creal(out[i]), cimag(out[i]));
    }
    fftw_destroy_plan(p);

    roundtrip = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, out, roundtrip, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    printf("roundtrip for %d:\n", N);
    for (int i=0; i<N; i++) {
        printf("%g + %gi ", creal(roundtrip[i]) / N, cimag(roundtrip[i]) / N);
    }
    fftw_destroy_plan(p);

    fftw_free(in); fftw_free(out); fftw_free(roundtrip);
}

