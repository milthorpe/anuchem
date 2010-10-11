package edu.mit.fftw;

import x10.compiler.Native;
import x10.compiler.NativeRep;
import x10.compiler.NativeCPPCompilationUnit;

@NativeCPPCompilationUnit("FFTW.cc")
@NativeRep("c++", "::edu::mit::fftw::FFTWWrapper", "::edu::mit::fftw::FFTWWrapper", null)
public class FFTW {
    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft1d(#1, reinterpret_cast<fftw_complex*>(#2._val->raw().data), reinterpret_cast<fftw_complex*>(#3._val->raw().data), #4)")
    @Native("java", "null")
    public native static def fftwPlan1d(n : Int, input : Array[Complex]{rail}, output : Array[Complex]{rail}, forward : Boolean) : FFTWPlan;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft3d(#1, #2, #3, reinterpret_cast<fftw_complex*>(#4._val->raw().data), reinterpret_cast<fftw_complex*>(#5._val->raw().data), #6)")
    @Native("java", "null")
    public native static def fftwPlan3d(n1 : Int, n2 : Int, n3 : Int, input : DistArray[Complex], output : DistArray[Complex], forward : Boolean) : FFTWPlan;
    
    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwExecute(#1)")
    @Native("java", "")
    public native static def fftwExecute(p : FFTWPlan) : Void;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwDestroyPlan(#1)")
    @Native("java", "")
    public native static def fftwDestroyPlan(p : FFTWPlan) : Void;

    // TODO threaded FFTW.  This is not needed ATM.
    /*
    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwInitThreads()")
    @Native("java", "0")
    public native static def fftwInitThreads() : Int;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanWithNThreads(#1)")
    @Native("java", "")
    public native static def fftwPlanWithNThreads(nThreads : Int) : Void;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwCleanupThreads()")
    @Native("java", "")
    public native static def fftwCleanupThreads() : Void;
    */

    @NativeRep("c++", "::edu::mit::fftw::FFTW_FFTWPlan", "::edu::mit::fftw::FFTW_FFTWPlan", null)
    @NativeRep("java", "edu.mit.fftw.FFTWPlan", "edu.mit.fftw.FFTWPlan", null)
    public static struct FFTWPlan { };
}
