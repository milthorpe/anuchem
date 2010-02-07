package edu.mit.fftw;

import x10.compiler.Native;
import x10.compiler.NativeRep;

@NativeRep("c++", "edu:mit:fftw:FFTWWrapper", "edu:mit:fftw:FFTWWrapper", null)
public class FFTW {
    @Native("c++", "FFTWWrapper::fftwPlanDft1d(#1, reinterpret_cast<fftw_complex*>(#2._val->_data), reinterpret_cast<fftw_complex*>(#3._val->_data), #4)")
    public native static def fftwPlan1D(n : Int, input : Rail[Complex], output : Rail[Complex], forward : Boolean) : FFTWPlan;

    @Native("c++", "FFTWWrapper::fftwExecute(#1)")
    public native static def fftwExecute(p : FFTWPlan) : Void;

    @Native("c++", "FFTWWrapper::fftwDestroyPlan(#1)")
    public native static def fftwDestroyPlan(p : FFTWPlan) : Void;

    @NativeRep("c++", "FFTW_FFTWPlan", "FFTW_FFTWPlan", null)
    public struct FFTWPlan { };
}
