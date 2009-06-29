package x10.lang;

/**
 * This class represents a complex number (a + b*i).
 * @author milthorpe
 */
public value Complex {
	public val real : Double;
	public val imaginary : Double;
	
	public static val ZERO : Complex = new Complex(0.0, 0.0);
    public static val ONE : Complex = new Complex(1.0, 0.0);
    public static val I : Complex = new Complex(0.0, 1.0);
	public static val INF : Complex = new Complex(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
    public static val NaN : Complex = new Complex(Double.NaN, Double.NaN);

	public def this(real : Double, imaginary : Double) {
		this.real = real;
		this.imaginary = imaginary;
	}
	
    /**
     * @return the sum of this complex number and the given complex number
     */
	public def add(a : Complex) : Complex {
		return new Complex(real + a.real, imaginary + a.imaginary);
	}

    /**
     * @return the sum of this complex number and the given Double
     */
	public def add(a : Double) : Complex {
		return new Complex(real + a, imaginary);
	}

    /**
     * @return the difference between this complex number and the given complex number
     */
    public def subtract(a : Complex) : Complex {
        if (isNaN() || a.isNaN()) {
            return NaN;
        }
        
        return new Complex(real - a.real, imaginary - a.imaginary);
    }

    /**
     * @return the difference between this complex number and the given Double
     */
	public def subtract(a : Double) : Complex {
        if (isNaN()) {
            return NaN;
        }
		return new Complex(real - a, imaginary);
	}
	
    /**
     * @return the product of this complex number and the given complex number
     */
	public def multiply(a : Complex) : Complex {
		return new Complex(real * a.real - imaginary * a.imaginary, 
                             real * a.imaginary + imaginary * a.real);
	}

    /**
     * @return the product of this complex number and the given Double
     */
	public def multiply(a : Double) : Complex {
		return new Complex(real * a, imaginary * a);
	}

    /**
     * Gets the quotient of this complex number and the given complex number.
	 * Uses Smith's algorithm <a href="http://doi.acm.org/10.1145/368637.368661"/>
	 * TODO: consider using Priest's algorithm <a href="http://doi.acm.org/10.1145/1039813.1039814"/>
     * @return the quotient of this complex number and the given complex number
     */
    public def divide(w : Complex) : Complex {
        if (isNaN() || w.isNaN()) {
            return NaN;
        }

        val c : Double = w.real;
        val d : Double = w.imaginary;
        if (c == 0.0 && d == 0.0) {
            return NaN;
        }
        
        if (w.isInfinite() && !isInfinite()) {
            return ZERO;
        }
		 	 
        if (Math.abs(d) <= Math.abs(c)) {
			if (c == 0.0) {
                return new Complex(imaginary/d, -real/c);
            }
            val r : Double =  d / c;
            val denominator : Double = c + d * r;
            return new Complex((real + imaginary * r) / denominator,
                (imaginary - real * r) / denominator);
        } else {
			if (d == 0.0) {
                return new Complex(real/c, imaginary/c);
            }
            val r : Double = c / d;
            val denominator : Double = c * r + d;
            return new Complex((real * r + imaginary) / denominator,
                (imaginary * r - real) / denominator);
        }
    }

    /**
     * Gets the quotient of this complex number and the given Double.
     */
    public def divide(a : Double) : Complex {
        return new Complex(real / a, imaginary / a);
    }
	
	public def equals(a : Complex) : boolean {
        if (a.isNaN()) return isNaN();
        else return real == a.real && imaginary == a.imaginary;
	}

    /** 
     * @return the conjugate of this complex number
     */
    public def conjugate() : Complex {
        if (isNaN()) {
            return NaN;
        }   
        return new Complex(real, -imaginary);
    }

    /**
     * @return the negation of this complex number
     */
    public def negate() : Complex {
        if (isNaN()) {
            return NaN;
        }
        
        return new Complex(-real, -imaginary);
    }

    /**
     * Return the absolute value of this complex number.
     * <p>
     * Returns <code>NaN</code> if either real or imaginary part is
     * <code>NaN</code> and <code>Double.POSITIVE_INFINITY</code> if
     * neither part is <code>NaN</code>, but at least one part takes an infinite
     * value.
     *
     * @return the absolute value
     */
    public def abs() : Double {
        if (isNaN()) {
            return Double.NaN;
        }
        
        if (isInfinite()) {
            return Double.POSITIVE_INFINITY;
        }

        if (imaginary == 0.0) {
            return Math.abs(real);
        } else if (real == 0.0) {
            return Math.abs(imaginary);
        } else {
            // use hypot to avoid unnecessary under/overflow
            return Math.hypot(real, imaginary);
        }
    }

    /**
     * @return true if either part of this complex number is NaN
     */
    public def isNaN() : boolean {
        return real.isNaN() || imaginary.isNaN();        
    }

    /**
     * @return true if either part of this complex number is infinite
     * and neither part is <code>NaN</code>
     */
    public def isInfinite() : boolean {
        return !isNaN() && 
        (real.isInfinite() || imaginary.isInfinite());
    }

	public def toString() : String {
		return (real + " + " + imaginary + "i");
	}
}
