public class Complex {
	public val imaginary : double;
    public val real : double;
	
	public static val ZERO : Complex = new Complex(0.0, 0.0);
    public static val ONE : Complex = new Complex(1.0, 0.0);
    public static val I : Complex = new Complex(0.0, 1.0);
    public static val NaN : Complex = new Complex(Double.NaN, Double.NaN);

	public def this(imaginary : double, real : double) {
		this.imaginary = imaginary; this.real = real;
	}
	
    /**
     * @return the sum of this complex number and the given complex number
     */
	public def add(a : Complex) : Complex {
		return new Complex(imaginary + a.imaginary, real + a.real);
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
     * @return the product of this complex number and the given complex number
     */
	public def multiply(a : Complex) : Complex {
		return new Complex(real * a.real - imaginary * a.imaginary, 
                             real * a.imaginary + imaginary * a.real);
	}

    /**
     * Gets the quotient of this complex number and the given complex number.
     * <p>
     * Uses 
     * <a href="http://doi.acm.org/10.1145/1039813.1039814">
     * </p>
     * @return the quotient of this complex number and the given complex number
     */
    public def divide(a : Complex) : Complex {
        if (isNaN() || a.isNaN()) {
            return NaN;
        }

        val c : double = a.real;
        val d : double = a.imaginary;
        if (c == 0.0 && d == 0.0) {
            return NaN;
        }
        
        if (a.isInfinite() && !isInfinite()) {
            return ZERO;
        }

        if (Math.abs(c) < Math.abs(d)) {
            if (d == 0.0) {
                return new Complex(real/c, imaginary/c);
            }
            val q : double = c / d;
            val denominator : double = c * q + d;
            return new Complex((real * q + imaginary) / denominator,
                (imaginary * q - real) / denominator);
        } else {
            if (c == 0.0) {
                return new Complex(imaginary/d, -real/c);
            }
            val q : double =  d / c;
            val denominator : double = d * q + c;
            return new Complex((imaginary * q + real) / denominator,
                (imaginary - real * q) / denominator);
        }
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
    public def abs() : double {
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
}
