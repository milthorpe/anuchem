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
package edu.mit.fftw;

import x10.compiler.Native;
import x10.compiler.NativeRep;
import x10.compiler.NativeCPPCompilationUnit;
import x10.compiler.NativeCPPInclude;
import x10.regionarray.DistArray;

@NativeCPPCompilationUnit("FFTW.cc")
@NativeRep("c++", "::edu::mit::fftw::FFTWWrapper", "::edu::mit::fftw::FFTWWrapper", null)
public class FFTW {
    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft1d(#1, reinterpret_cast<fftw_complex*>(#2->raw), reinterpret_cast<fftw_complex*>(#3->raw), #4)")
    @Native("java", "null")
    public native static def fftwPlan1d(n:Int, input:Rail[Complex], output:Rail[Complex], forward:Boolean):FFTWPlan;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft1d(#1, reinterpret_cast<double*>(#2->raw), reinterpret_cast<fftw_complex*>(#3->raw))")
    @Native("java", "null")
    public native static def fftwPlan1d(n:Int, input:Rail[Double], output:Rail[Complex]):FFTWPlan;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft1d(#1, reinterpret_cast<fftw_complex*>(#2->raw), reinterpret_cast<double*>(#3->raw))")
    @Native("java", "null")
    public native static def fftwPlan1d(n:Int, input:Rail[Complex], output:Rail[Double]):FFTWPlan;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft1d(#1, #2, reinterpret_cast<fftw_complex*>(#3->raw()->raw), reinterpret_cast<fftw_complex*>(#4->raw()->raw), #5)")
    @Native("java", "null")
    public native static def fftwPlan1d(n:Int, howMany:Int, input:DistArray[Complex], output:DistArray[Complex], forward:Boolean):FFTWPlan;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft1d(#1, #2, reinterpret_cast<double*>(#3->raw()->raw), reinterpret_cast<fftw_complex*>(#4->raw()->raw))")
    @Native("java", "null")
    public native static def fftwPlan1d(n:Int, howMany:Int, input:DistArray[Double], output:DistArray[Complex]):FFTWPlan;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft1d(#1, #2, reinterpret_cast<fftw_complex*>(#3->raw()->raw), reinterpret_cast<double*>(#4->raw()->raw))")
    @Native("java", "null")
    public native static def fftwPlan1d(n:Int, howMany:Int, input:DistArray[Complex], output:DistArray[Double]):FFTWPlan;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft3d(#1, #2, #3, reinterpret_cast<fftw_complex*>(#4->raw()->raw), reinterpret_cast<fftw_complex*>(#5->raw()->raw), #6)")
    @Native("java", "null")
    public native static def fftwPlan3d(n1:Int, n2:Int, n3:Int, input:DistArray[Complex], output:DistArray[Complex], forward:Boolean):FFTWPlan;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft3d(#1, #2, #3, reinterpret_cast<double*>(#4->raw()->raw), reinterpret_cast<fftw_complex*>(#5->raw()->raw))")
    @Native("java", "null")
    public native static def fftwPlan3d(n1:Int, n2:Int, n3:Int, input:DistArray[Double], output:DistArray[Complex]):FFTWPlan;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanDft3d(#1, #2, #3, reinterpret_cast<fftw_complex*>(#4->raw()->raw), reinterpret_cast<double*>(#5->raw()->raw))")
    @Native("java", "null")
    public native static def fftwPlan3d(n1:Int, n2:Int, n3:Int, input:DistArray[Complex], output:DistArray[Double]):FFTWPlan;
    
    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwExecute(#1)")
    @Native("java", "")
    public native static def fftwExecute(p:FFTWPlan):void;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwDestroyPlan(#1)")
    @Native("java", "")
    public native static def fftwDestroyPlan(p:FFTWPlan):void;

    // TODO threaded FFTW.  This is not needed ATM.
    /*
    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwInitThreads()")
    @Native("java", "0")
    public native static def fftwInitThreads():Int;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwPlanWithNThreads(#1)")
    @Native("java", "")
    public native static def fftwPlanWithNThreads(nThreads:Int):void;

    @Native("c++", "::edu::mit::fftw::FFTWWrapper::fftwCleanupThreads()")
    @Native("java", "")
    public native static def fftwCleanupThreads():void;
    */

    @NativeRep("c++", "::edu::mit::fftw::FFTW_FFTWPlan", "::edu::mit::fftw::FFTW_FFTWPlan", null)
    @NativeRep("java", "edu.mit.fftw.FFTWPlan", "edu.mit.fftw.FFTWPlan", null)
    public static struct FFTWPlan { };
}
