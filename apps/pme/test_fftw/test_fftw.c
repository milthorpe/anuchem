#include <stdio.h>
#include <complex.h>
#include <fftw3.h>

using namespace std;

int main() {
    fftw_complex *in, *out, *roundtrip;
    fftw_plan p;
    int N = 50;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    complex<double> onePlusi = complex<double>(1.0, 1.0);
    for (int i=0; i<N; i++) {
        complex<double> val = 1.0 / (i+1) * onePlusi;
        in[i][0] = real(val);
        in[i][1] = imag(val);
        printf("%g + %gi ", in[i][0], in[i][1]);
    }
    fftw_execute(p);
    printf("result for %d:\n", N);
    for (int i=0; i<N; i++) {
        printf("%g + %gi ", in[i][0], in[i][1]);
    }
    fftw_destroy_plan(p);

    roundtrip = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, out, roundtrip, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    printf("roundtrip for %d:\n", N);
    for (int i=0; i<N; i++) {
        printf("%g + %gi ", in[i][0] / N, in[i][1] / N);
    }
    fftw_destroy_plan(p);

    fftw_free(in); fftw_free(out); fftw_free(roundtrip);
}

