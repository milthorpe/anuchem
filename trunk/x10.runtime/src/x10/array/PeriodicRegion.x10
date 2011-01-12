/*
 *  This file is part of the X10 project (http://x10-lang.org).
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Australian National University 2010.
 */

package x10.array;

import x10.array.Region.Scanner;
import x10.compiler.Inline;
import x10.compiler.TempNoInline_1;

/**
 * A periodic region decorates a standard X10 region by implementing 
 * periodic boundary conditions, in which elements at each edge of
 * the region are considered to be neighbours, and indexes that fall
 * outside the "home" region in any dimension are wrapped around modulo 
 * the size of the region in that dimension.
 */
public final class PeriodicRegion extends Region{rect} {
    public val baseRegion : RectRegion{self.rank==this.rank};

    public def this(base : RectRegion) : PeriodicRegion(base.rank) {
        super(base.rank, base.rect, base.zeroBased);
        this.baseRegion = base;
    }

    private @Inline def getPeriodicIndex(index : Int, dim : Int) : Int {
        val regionMin = baseRegion.min(dim);
        val regionMax = baseRegion.max(dim);
        var actualIndex : Int = index;
        while (actualIndex < regionMin) actualIndex += regionMax-regionMin+1;
        while (actualIndex > regionMax) actualIndex -= regionMax-regionMin+1; 
        return actualIndex;
    }

    public @Inline def indexOf(pt:Point):Int {
        val actualPt = Point.make(new Array[Int](rank, (i : Int) => getPeriodicIndex(pt(i), i)));
	    return baseRegion.indexOf(actualPt);
    }

    public @Inline def indexOf(i0:int)  = baseRegion.indexOf(getPeriodicIndex(i0, 0));
    public @Inline def indexOf(i0:int, i1:int) = baseRegion.indexOf(getPeriodicIndex(i0, 0), getPeriodicIndex(i1, 1));
    public @Inline def indexOf(i0:int, i1:int, i2:int) = baseRegion.indexOf(getPeriodicIndex(i0, 0), getPeriodicIndex(i1, 1), getPeriodicIndex(i2, 2));
    public @Inline def indexOf(i0:int, i1:int, i2:int, i3:int) = baseRegion.indexOf(getPeriodicIndex(i0, 0), getPeriodicIndex(i1, 1), getPeriodicIndex(i2, 2), getPeriodicIndex(i3, 3));

    public @Inline def size(): int = baseRegion.size();
    public @Inline def isConvex(): boolean = baseRegion.isConvex();
    public @Inline def isEmpty(): boolean = baseRegion.isEmpty();
    protected def computeBoundingBox(): Region(rank) = baseRegion.computeBoundingBox();
    public @TempNoInline_1 def min():(int)=>int = (i:int)=> baseRegion.min(i);
    public @TempNoInline_1 def max():(int)=>int = (i:int)=> baseRegion.max(i);
    public def intersection(that: Region(rank)): Region(rank) = new PeriodicRegion(baseRegion.intersection(that) as RectRegion);
    public def product(that: Region): Region = new PeriodicRegion(baseRegion.product(that) as RectRegion);
    public def translate(v: Point(rank)): Region(rank) = new PeriodicRegion(baseRegion.translate(v) as RectRegion);
    public def projection(axis: int): Region(1) = new PeriodicRegion(baseRegion.projection(axis) as RectRegion);
    public def eliminate(axis: int): Region /*(rank-1)*/ = new PeriodicRegion(baseRegion.eliminate(axis) as RectRegion);
    public def iterator(): Iterator[Point(rank)] = baseRegion.iterator();
    public def scanners(): Iterator[Scanner] = baseRegion.scanners();
    public def contains(that: Region(rank)): boolean = baseRegion.contains(that);
    public def contains(p:Point):boolean = baseRegion.contains(p);

    public def toString():String {
        return "Periodic: " + baseRegion.toString();
    }
}

