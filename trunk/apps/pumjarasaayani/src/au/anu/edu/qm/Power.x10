package au.anu.edu.qm;

public class Power { 
    var l:Int, m:Int, n:Int;

    public def this() { l = m = n = 0; }

    public def this(l:Int, m:Int, n:Int) { 
        this.l = l; this.m = m; this.n = n;
    } 

    public def getL() : Int = this.l;
    public def getM() : Int = this.m;
    public def getN() : Int = this.n;
}

