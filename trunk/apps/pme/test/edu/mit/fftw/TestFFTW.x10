package edu.mit.fftw;

public class TestFFTW {
    public static def main(args : Rail[String]!) {
        val N : Int = 50;
        val input = Rail.make[Complex](N, (i : Int) => 1.0 / (i+1.0) * (1.0 + Complex.I));
        for (var i : int = 0; i < input.length; i++) {
            Console.OUT.println(input(i));
        }
        val output = Rail.make[Complex](N);
        val plan : FFTW.FFTWPlan = FFTW.fftwPlan1D(N, input, output, false);
        FFTW.fftwExecute(plan);
        for (var i : int = 0; i < output.length; i++) {
            Console.OUT.println("output(" + i + ") =" + output(i));
        }
        FFTW.fftwDestroyPlan(plan);
        val roundtrip = Rail.make[Complex](N);
        val plan2 : FFTW.FFTWPlan = FFTW.fftwPlan1D(N, output, roundtrip, true);
        FFTW.fftwExecute(plan2);
        for (var i : int = 0; i < roundtrip.length; i++) {
            Console.OUT.println("roundtrip(" + i + ") =" + roundtrip(i) / N);
        }
    }
}
